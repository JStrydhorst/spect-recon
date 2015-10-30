#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <ctime>

#include "fft.h"
#include "dicom.h"
#include "proj_table.h"
#include "recon4d.h"


#ifdef _WIN32
	#define __SL "\\"
#elif defined __unix
	#define __SL "/"
#endif

//using namespace std;

typedef unsigned char uchar;

void resort(int* list, int n);
float fit_cosine(float *input, int N, float *coeffs);

/* Gaussian Rotator Parameters */
#define ROT_SUBPIX	65
#define ROT_MAT_SIZE	3 // must be an odd number (3 or 5)
#define ROT_FWHM	1.0f

/* Coordinate system
Internally, x is increasing left-to-right, y is increasing bottom-to-top, and
z is directed into the bore (left-hand coordinate system). 
Rotations are specified in the couterclockwise direction, with 0 defined as the positive x-axis.

	y
	^
	|
	|
	|
	|
	0-----------> x

The DICOM standard expects the top row first, the order of the y-vector must be flipped when saving DICOM files

	0-----------> x
	|
	|
	|
	|
	v
	y


/********************************************************************************************
/ Attenuation map
/********************************************************************************************/
class AttenMap
{
public:
	AttenMap();
	AttenMap(const char* filename, int new_num_proj, int new_nz, int new_nxy);	// create attenuation map from attenmap
	~AttenMap();

	friend class Reconstruction;
private:
	int num_proj;
	int nxy, nz;
	float**** atten_map;
};

AttenMap::AttenMap()
{
	num_proj = 0;
	nz = 0;
	nxy = 0;
	atten_map = NULL;
}

AttenMap::AttenMap(const char *filename, int new_num_proj, int new_nz, int new_nxy)
{	
	int angle, i, j;
	ifstream fin;

	num_proj = new_num_proj;
	nz = new_nz;
	nxy = new_nxy;

	try
	{
		atten_map = new float***[num_proj];
		for(angle=0;angle<num_proj;angle++)
			atten_map[angle] = new float**[nz];
		for(angle=0;angle<num_proj;angle++)
			for(i=0;i<nz;i++)
				atten_map[angle][i] = new float*[nxy];
		for(angle=0;angle<num_proj;angle++)
			for(i=0;i<nz;i++)
				for(j=0;j<nxy;j++)
					atten_map[angle][i][j] = new float[nxy];
	}
	catch(bad_alloc& ba)
	{
		cout << "bad_alloc caught in AttenMap::AttenMAP  --  " << ba.what() << endl;
		exit(-1);
	}

	fin.open(filename,ifstream::binary);
	for(angle=0;angle<num_proj;angle++)
		for(i=0;i<nz;i++)
			for(j=0;j<nxy;j++)
				fin.read(reinterpret_cast<char*>(atten_map[angle][i][j]), nxy*sizeof(float));
	fin.close();
}

AttenMap::~AttenMap()
{
	int angle,i,j;
	for(angle = 0;angle<num_proj;angle++)
		for(i=0;i<nz;i++)
			for(j=0;j<nxy;j++)
				delete [] atten_map[angle][i][j];
	for(angle = 0;angle<num_proj;angle++)
		for(i=0;i<nz;i++)
			delete [] atten_map[angle][i];
	for(angle = 0;angle<num_proj;angle++)
		delete [] atten_map[angle];
	delete [] atten_map;
}

/********************************************************************************************
/ Projection
/********************************************************************************************/
Projection::Projection()	// default constructon, empty projection, is this ever called???
{
	int i;

	n_timeslots = 1;
	ny = 0;
	nz = 0;
	res = 0;
	num_proj = 0;
	angles = NULL;
	z_shift = NULL;
	detector = NULL;
	angular_view = NULL;
	timeslot = NULL;
	frame_time = NULL;
	_data = NULL;
	data = NULL;

	for(i=0;i<4;i++)
	{
		AxialDetectorOffset[i] = 0.0;
		TransaxialDetectorOffset[i] = 0.0;

		AxialApertureOffset[i] = 0.0;
		TransaxialApertureOffset[i] = 0.0;
		RadiusOfRotation[i] = 0.0;
	}

	aperture[0] = 0;
}

// creates blank projection with specified parameters, creates a uniformly spaced helical scan, single time frame
// (including both aperture_name and n_pinholes overspecifies the geometry)
Projection::Projection(int projections, int size_y, int size_z, float new_res,
					   char* aperture_name, const char* filename)
{
	int i,j;
	ifstream fin;

	n_timeslots = 1;
	ny = size_y;
	nz = size_z;
	res = new_res;
	num_proj = projections;
	if(aperture_name)
		strcpy(aperture, aperture_name);

	// initialize angle array to evenly spaced projections starting at 90 degrees
	angles = new float[num_proj];
	for(i=0;i<num_proj;i++)
		angles[i] = i * (360.0f/num_proj);// - 90.0f;	

	// initialize z_shift array to 0 (circular scan)
	z_shift = new float[num_proj];
	for(i=0;i<num_proj;i++)
		z_shift[i] = 0.0f;

	// initialize z_shift array to 0 (circular scan)
	detector = new unsigned short[num_proj];
	for(i=0;i<num_proj;i++)
		detector[i] = 1;

	for(i=0;i<4;i++)
	{
		AxialDetectorOffset[i] = 0.0;
		TransaxialDetectorOffset[i] = 0.0;
		AxialApertureOffset[i] = 0.0;
		TransaxialApertureOffset[i] = 0.0;
		RadiusOfRotation[i] = 45.0;
	}
	
	angular_view = new unsigned short[num_proj];
	for(i=0;i<num_proj;i++)
		angular_view[i] = i;

	timeslot = new unsigned short[num_proj];
	for(i=0;i<num_proj;i++)
		timeslot[i] = 1;

	frame_time = new float[num_proj];
	for(i=0;i<num_proj;i++)
		frame_time[i] = 0.0f;

	try
	{
		// allocate memory for data
		_data = new float[num_proj*nz*ny];
		data = new float**[num_proj];
		for(i=0;i<num_proj;i++)
			data[i] = new float*[nz];
		for(i=0;i<num_proj;i++)
			for(j=0;j<nz;j++)
				data[i][j] = _data + i*nz*ny + j*ny;
	}
	catch(bad_alloc& ba)
	{
		cout << "bad_alloc caught in Projection::Projection(int, int, int, etc... )  --  " << ba.what() << endl;
		exit(-1);
	}

	if(filename)
	{
		fin.open(filename,ios::binary);
		for(i=0;i<num_proj;i++)
			for(j=0;j<nz;j++)
				fin.read(reinterpret_cast<char*>(data[i][j]), ny*sizeof(float));
		fin.close();
	}
	else
		Clear();

}

Projection::Projection(const char* DCMfilename)
{
	int i,j,k;
	int len;
	unsigned short temp_us;
	char temp_str[64];
	unsigned short* temp_data;
	char *temp_ch, *p_ch;
	float max_z_shift;

	float start_angle[4];	// number of detectors is hard coded
	float angular_step;

	int hour, min;

	RootDicomObj* DCM = new RootDicomObj(DCMfilename);
	DicomObj* tempDCM;

	DCM->GetValue(0x0028,0x0010,&temp_us,sizeof(temp_us));
	nz = temp_us;
	DCM->GetValue(0x0028,0x0011,&temp_us,sizeof(temp_us));
	ny = temp_us;

	if(DCM->GetValue(0x0054,0x0071,&temp_us,sizeof(temp_us)) != 0xFFFFFFFF)
		n_timeslots = temp_us;
	else
		n_timeslots = 1;

	len = DCM->GetValue(0x0028,0x0030,temp_str,64);
	temp_str[len] = 0;
	res = float(atof(temp_str));
	res = float(floor(res*10 + 0.5)/10); // rounds to a precision of 0.1 mm

	len = DCM->GetValue(0x0028,0x0008,temp_str,64);
	temp_str[len] = 0;
	num_proj = atoi(temp_str);

	tempDCM = DCM->GetSQObject(0x0054,0x0022,0);
	if(tempDCM)
	{
		len = tempDCM->GetValue(0x0009,0x1032,aperture,15);
		aperture[len] = 0;
	}
	else
		sprintf(aperture,"APT2");

	if(!strncmp(aperture,"APT3",5))	// bug in some data sets reports aperture APT2 as APT3
	{
		cout << "Fixing aperture name!" << endl;
		sprintf(aperture,"APT2");
	}

	try
	{
		// allocate memory for data
		_data = new float[num_proj*nz*ny];
		data = new float**[num_proj];
		for(i=0;i<num_proj;i++)
			data[i] = new float*[nz];
		for(i=0;i<num_proj;i++)
			for(j=0;j<nz;j++)
				data[i][j] = _data + i*nz*ny + j*ny;
	}
	catch(bad_alloc& ba)
	{
		cout << "bad_alloc caught in Projection::Projection(const char*)  --  " << ba.what() << endl;
		exit(-1);
	}

	// read in projection data
	len = num_proj * nz * ny;
	temp_data = new unsigned short[len];
	DCM->GetValue(0x7FE0,0x0010,temp_data,len*sizeof(unsigned short));
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
				data[i][j][k] = temp_data[i*nz*ny+j*ny+k];
	delete [] temp_data;

	// scale projection data
	float rescale, offset;
	memset(temp_str,0,sizeof(temp_str));
	if(DCM->GetValue(0x0028,0x1053,temp_str,sizeof(temp_str)) != -1L)
	{
		rescale = float(atof(temp_str));

		memset(temp_str,0,sizeof(temp_str));
		DCM->GetValue(0x0028,0x1052,temp_str,sizeof(temp_str));
		offset = float(atof(temp_str));

		for(i=0;i<num_proj;i++)
			for(j=0;j<nz;j++)
				for(k=0;k<ny;k++)
					data[i][j][k] = offset + rescale * data[i][j][k];
	}

	// read detector number
	detector = new unsigned short[num_proj];
	DCM->GetValue(0x0054,0x0020,detector,num_proj*sizeof(unsigned short));

	// read start angles and offsets
	for(i=0;i<4;i++)
	{
		tempDCM = DCM->GetSQObject(0x0054,0x0022, i);
		len = tempDCM->GetValue(0x0054,0x0200,temp_str,64);
		temp_str[len] = 0;
		start_angle[i] = float(270 - atoi(temp_str));

		tempDCM->GetValue(0x0009,0x107d,&TransaxialApertureOffset[i],sizeof(double));
		tempDCM->GetValue(0x0009,0x107e,&AxialApertureOffset[i],sizeof(double));

		tempDCM->GetValue(0x0009,0x107f,&TransaxialDetectorOffset[i],sizeof(double));
		tempDCM->GetValue(0x0009,0x1080,&AxialDetectorOffset[i],sizeof(double));

		tempDCM->GetValue(0x0009,0x1081,&RadiusOfRotation[i],sizeof(double));
	}

	// read increment
	tempDCM = DCM->GetSQObject(0x0054,0x0052);
	len = tempDCM->GetValue(0x0018,0x1144,temp_str,64);
	temp_str[len] = 0;
	angular_step = float(atof(temp_str));

	angles = new float[num_proj];
	angular_view = new unsigned short[num_proj];
	DCM->GetValue(0x0054,0x0090,angular_view,num_proj*sizeof(unsigned short));
	for(i=0;i<num_proj;i++)
	{
		angles[i] =  start_angle[detector[i]-1] + angular_step*(angular_view[i] - 1);
		// following is completely unncessary, but it makes the numbers fall into [0, 360)
		while (angles[i] < 0.0)
			angles[i] += 360.0f;
		while(angles[i] >= 360.0f)
			angles[i] -= 360.0f;
	}

	// set up z_shift for helical scan
	z_shift = new float[num_proj];
	tempDCM = DCM->GetSQObject(0x0054,0x0052);
	len = tempDCM->GetLength(0x0009,0x1024);
	temp_ch = new char[len];
	len = tempDCM->GetValue(0x0009,0x1024,temp_ch,len);
	p_ch = temp_ch;

	max_z_shift = 0.0f;
	for(i=0;i<num_proj;i++)
	{
		z_shift[i] = float(atof(p_ch));
		p_ch = strchr(p_ch,'\\') + 1;

		if(i)
		{
			if (z_shift[i] > max_z_shift)
				max_z_shift = z_shift[i];
		}
		else
			max_z_shift = z_shift[i];
	}
	for(i=0;i<num_proj;i++)
		z_shift[i] = max_z_shift - z_shift[i];	// shift to get a zero-based z_shift
	delete [] temp_ch;

	// initialize time_slot vector
	timeslot = new unsigned short[num_proj];
	if (n_timeslots > 1)
		DCM->GetValue(0x0054,0x0070,timeslot,num_proj*sizeof(unsigned short));
	else
		for(i=0;i<num_proj;i++)
			timeslot[i] = 1;

	// frame time vector
	frame_time = new float[num_proj];
	tempDCM = DCM->GetSQObject(0x0054,0x0052);
	len =  tempDCM->GetLength(0x0009,0x1027);
	temp_ch = new char[len];
	len = tempDCM->GetValue(0x0009,0x1027,temp_ch,len);
	p_ch = temp_ch;
	for(i=0;i<num_proj;i++)
	{
		strncpy(temp_str,p_ch+8,2);
		temp_str[2] = '\0';
		hour = atoi(temp_str);
		strncpy(temp_str,p_ch+10,2);
		temp_str[2] = '\0';
		min = atoi(temp_str);
		frame_time[i] = float(hour*60 + min);
		p_ch = strchr(p_ch,'\\') + 1;

		if(i>0)
			frame_time[i] -= frame_time[0];
	}
	frame_time[0] = 0.0f;
	delete [] temp_ch;

	delete DCM;
}

// creates a copy of the projection
Projection::Projection(const Projection &param)
{
	int i,j;

	n_timeslots = param.n_timeslots;
	ny = param.ny;
	nz = param.nz;
	res = param.res;
	num_proj = param.num_proj;
	angles = new float[num_proj];
	z_shift = new float[num_proj];
	detector = new unsigned short[num_proj];
	angular_view = new unsigned short[num_proj];
	timeslot = new unsigned short[num_proj];
	frame_time = new float[num_proj];

	for(i=0; i<num_proj;i++)
	{
		angles[i] = param.angles[i];
		z_shift[i] = param.z_shift[i];
		detector[i] = param.detector[i];
		angular_view[i] = param.angular_view[i];
		timeslot[i] = param.timeslot[i];
		frame_time[i] = param.frame_time[i];
	}

	for(i=0;i<4;i++)
	{
		AxialDetectorOffset[i] = param.AxialDetectorOffset[i];
		TransaxialDetectorOffset[i] = param.TransaxialDetectorOffset[i];
		AxialApertureOffset[i] = param.AxialApertureOffset[i];
		TransaxialApertureOffset[i] = param.TransaxialApertureOffset[i];
		RadiusOfRotation[i] = param.RadiusOfRotation[i];
	}

	_data = new float[num_proj*nz*ny];
	data = new float**[num_proj];
	for(i=0;i<num_proj;i++)
		data[i] = new float*[nz];
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			data[i][j] = _data + i*nz*ny + j*ny;

	memcpy(aperture,param.aperture,16);

	// copy data
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			memcpy(data[i][j],param.data[i][j],ny*sizeof(float));

}

Projection::~Projection()
{
	int i;

	for(i=0;i<num_proj;i++)
		delete [] data[i];
	delete [] data;
	delete [] _data;

	delete [] angles;
	delete [] z_shift;
	delete [] detector;
	delete [] angular_view;
	delete [] timeslot;
	delete [] frame_time;
}

void Projection::ProjectionInfo()
{
	cout << "******Projection Info************" << endl;
	cout << "Number of projections: " << num_proj << endl;
	cout << nz << " x " << ny << " (" << res << " mm)" << endl;
	cout << "Aperture: " << aperture << endl;
	cout << "Scan Length: " << GetScanLength() << " mm" << endl;
	cout << "Scan Duration: " << GetScanDuration() << " min" << endl;
	if(n_timeslots>1)
		cout << "Gates: " << GetNumTimeSlots() << endl;
	cout << "*********************************" << endl << endl;
}

// only works when initial projection set has 1 timeslot
void Projection::AddTimeslots(int n_new_timeslots)
{
	int i, j, det, ts, pos;
	int n_positions = num_proj/4;
	int new_offset;
	int old_offset;

	int n_new_proj = num_proj * n_new_timeslots;

	float* new_angles = new float[n_new_proj];
	float *new_z_shift = new float[n_new_proj];		// z_shift of projection
	unsigned short *new_detector = new unsigned short[n_new_proj];		// array of detector heads
	unsigned short *new_angular_view = new unsigned short[n_new_proj];
	float *new_frame_time = new float[n_new_proj];	// in minutes

	delete [] timeslot;
	timeslot = new unsigned short[n_new_proj];
	
	for(det=0;det<4;det++)
		for(ts=0;ts<n_new_timeslots;ts++)
			for(pos=0;pos<n_positions;pos++)
			{
				new_offset = det*n_positions*n_new_timeslots + ts*n_positions + pos;
				old_offset = det*n_positions + pos;

				new_angles[new_offset]=angles[old_offset];
				new_z_shift[new_offset]=z_shift[old_offset];
				new_detector[new_offset]=detector[old_offset];
				new_angular_view[new_offset]=angular_view[old_offset];
				new_frame_time[new_offset]=frame_time[old_offset];

				timeslot[new_offset] = ts+1;
			}
	delete [] angles;
	delete [] z_shift;
	delete [] detector;
	delete [] angular_view;
	delete [] frame_time;

	angles = new_angles;
	z_shift = new_z_shift;
	detector = new_detector;
	angular_view = new_angular_view;
	frame_time = new_frame_time;
	n_timeslots = n_new_timeslots;

	// clear the old memory
	for(i=0;i<num_proj;i++)
		delete [] data[i];
	delete [] data;
	delete [] _data;

	// update the number of projections
	num_proj = n_new_proj;

	// reallocate memory
	try
	{
		// allocate memory for data
		_data = new float[num_proj*nz*ny];
		data = new float**[num_proj];
		for(i=0;i<num_proj;i++)
			data[i] = new float*[nz];
		for(i=0;i<num_proj;i++)
			for(j=0;j<nz;j++)
			{
				data[i][j] = _data + i*nz*ny + j*ny;
				memset(data[i][j],0,ny*sizeof(float));
			}
	}
	catch(bad_alloc& ba)
	{
		cout << "bad_alloc caught in Projection::Projection(int, int, int, etc... )  --  " << ba.what() << endl;
		exit(-1);
	}

}


// Changes the timeslot setting, necessary because the same temporary projection is reused in the reconstruction algorithm
void Projection::SetTimeslot(int timeslot)
{
	for(int i=0;i<num_proj;i++)
		this->timeslot[i] = timeslot;
}

void Projection::Clear()
{
	int i, j, k;
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
				data[i][j][k] = 0.0;
}

// divides the values stored in the supplied parameter
// by the values of the current object and saves the
// results in the current object
void Projection::Ratio(const Projection* param)
{
	int i, j, k;
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
			{				
				if(data[i][j][k] > 0.01f)
					data[i][j][k] = param->data[i][j][k]/data[i][j][k];
				else
					data[i][j][k] = 1.0f;
			}
}


void Projection::Sum(const Projection* param, float frac)
{
	int i, j, k;
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
					data[i][j][k] += frac*param->data[i][j][k];
}


// intrinsic resolution filtering
// applies a gaussian filter to all of the projections
void Projection::Filter(float fwhm)
{
	int i,j,k;
	float **filter;
	float **data_i;	// imaginary data;
	float sigma, a;

	// create filter
	filter = new float*[nz];
	for(i=0;i<nz;i++)
		filter[i] = new float[ny];

	sigma = fwhm/2.35482f;	// 2*sqrt(2*ln(2))
	a = float(M_PI * M_PI * 2 * sigma * sigma / (res * res * ny * ny));
	for(i=0;i<=nz/2;i++)
	{
		for(j=0;j<=ny/2;j++)
			filter[i][j] = exp(-a*(i*i+j*j));
		for(;j<ny;j++)
			filter[i][j] = filter[i][ny-j];
	}
	for(;i<nz;i++)
		for(j=0;j<ny;j++)
			filter[i][j] = filter[nz-i][j];

	// set up array for imaginary data
	data_i = new float*[nz];
	for(i=0;i<nz;i++)
		data_i[i] = new float[ny];

	for(i=0;i<num_proj;i++)
	{
		// clear imaginary data set
		for(j=0;j<nz;j++)
			memset(data_i[j],0,ny*sizeof(float));
		// fft
		fft2(data[i],data_i,ny,nz,1);

		// apply filter
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
			{
				data[i][j][k] *= filter[j][k];
				data_i[j][k] *= filter[j][k];
			}

		// inverse fft
		fft2(data[i],data_i,ny,nz,-1);
	}


	// free memory
	for(i=0;i<nz;i++)
	{
		delete [] filter[i];
		delete [] data_i[i];
	}
	delete filter;
	delete data_i;

	// ensure non-negativity
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
				if(data[i][j][k] < 0)
					data[i][j][k] = 0;
}


void Projection::ApplyOffsets(int dir)
{
	int i,j,k;

	int j_src, k_src;
	float j_fr, k_fr;

	float j_shift, k_shift;

	float **temp;

	temp = new float*[nz];
	for(i=0;i<nz;i++)
		temp[i] = new float[ny];
	

	for(i=0;i<num_proj;i++)
	{
		if(dir)
		{
			j_shift = float(AxialApertureOffset[detector[i]-1] - AxialDetectorOffset[detector[i]-1]) / res;
			k_shift = float(TransaxialApertureOffset[detector[i]-1] - TransaxialDetectorOffset[detector[i]-1]) / res;
		}
		else
		{
			j_shift = float(AxialDetectorOffset[detector[i]-1] - AxialApertureOffset[detector[i]-1]) / res;
			k_shift = float(TransaxialDetectorOffset[detector[i]-1] - TransaxialApertureOffset[detector[i]-1]) / res;
		}

		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
			{
				temp[j][k] = data[i][j][k];
				data[i][j][k] = 0.0f;
			}

		for(j=0;j<nz;j++)
		{
			j_fr = j + j_shift;

			j_src = int(j_fr);
			if(j_src<0)
				continue;
			if(j_src>(nz-2))
				continue;
			j_fr -= j_src;

			for(k=0;k<ny;k++)
			{
				k_fr = k + k_shift;

				k_src = int(k_fr);
				if(k_src<0)
					continue;
				if(k_src>(ny-2))
					continue;
				k_fr -= k_src;

				data[i][j][k] = temp[j_src][k_src] * (1-j_fr) * (1-k_fr) +
								temp[j_src+1][k_src] * j_fr * (1-k_fr) +
								temp[j_src][k_src+1] * (1-j_fr) * k_fr +
								temp[j_src+1][k_src+1] * j_fr * k_fr;

			}
		}
	}

	for(i=0;i<nz;i++)
		delete [] temp[i];
	delete [] temp;

}

Projection* Projection::CreateNorm(float halflife)
{
	int i, j, k;

	Projection* newProj = new Projection(*this);

	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
				newProj->data[i][j][k] = 1.0f;

	if(halflife != 0.0f)
		newProj->ApplyDecay(halflife);

	return newProj;
}

/*
void Projection::ApplyCos3Scaling()
{
	int i, j, k;

	float y_ctr, z_ctr;
	float y_pos, z_pos;
	float temp;

	z_ctr = (nz - 1.0f)/2;
	y_ctr = (ny - 1.0f)/2;

	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
		{
			z_pos = (j-z_ctr) * res;
			for(k=0;k<ny;k++)
			{
				y_pos = (k-y_ctr) * res;
				temp = sqrt(176.0f * 176.0f + z_pos * z_pos + y_pos * y_pos)/176.0f;
				temp = temp*temp*temp;
				data[i][j][k] *= temp;
			}
		}
	
}
*/

void Projection::ApplyDecay(float halflife)
{
	int i,j,k;
	float decay;

	for(i=0;i<num_proj;i++)
	{
		decay = exp(-log(2.0f)*frame_time[i]/halflife);

		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
				data[i][j][k] *= decay;
	}

}


void Projection::Scale(float f)
{
	int i,j,k;

	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
				data[i][j][k] *= f;

}

/* need to find a library with a poisson deviate generator
void Projection::Poisson(long seed)
{
	int i,j,k;

	long idum=seed;

	std::default_random_engine generator;
	std::poisson_distribution<int> poidev;

	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
				data[i][j][k] = poidev(data[i][j][k]);

}
*/

// saves the data in a binary file...
void Projection::SaveBIN(const char* filename)
{
	int i,j;

	ofstream fout;
	fout.open(filename,fstream::binary);
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			fout.write(reinterpret_cast<char*>(data[i][j]),ny*sizeof(float));
	fout.close();

}

void Projection::LoadBIN(const char* filename)
{
	int i,j;

	ifstream fout;
	fout.open(filename,fstream::binary);
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			fout.read(reinterpret_cast<char*>(data[i][j]),ny*sizeof(float));
	fout.close();
}

Projection** Projection::SplitTimeSlots()
{
	int i,j,k;
	int current_timeslot;	// timeslot of the current frame (starts at 0)
	int dest_frame;			// frame in destination projection

	int *dest_frame_counter;
	dest_frame_counter = new int[n_timeslots];	// 

	Projection **proj_array;

	// create array of projections
	proj_array = new Projection*[n_timeslots];
	for (i=0;i<n_timeslots;i++)
	{
		proj_array[i] = new Projection(num_proj/n_timeslots, ny, nz, res, aperture);
		dest_frame_counter[i] = 0;

		for(j=0;j<4;j++)
		{
			proj_array[i]->AxialDetectorOffset[j] = AxialDetectorOffset[j];
			proj_array[i]->TransaxialDetectorOffset[j] = TransaxialDetectorOffset[j];

			proj_array[i]->AxialApertureOffset[j] = AxialApertureOffset[j];
			proj_array[i]->TransaxialApertureOffset[j] = TransaxialApertureOffset[j];

			proj_array[i]->RadiusOfRotation[j] = RadiusOfRotation[j];
		}
	}

	// copy data
	for(i=0;i<num_proj;i++)
	{
		current_timeslot = timeslot[i]-1;	// timeslots are numbered 1, 2, etc
		dest_frame = dest_frame_counter[current_timeslot];

		proj_array[current_timeslot]->timeslot[dest_frame] = timeslot[i];

		proj_array[current_timeslot]->angles[dest_frame] = angles[i];
		proj_array[current_timeslot]->z_shift[dest_frame] = z_shift[i];
		proj_array[current_timeslot]->detector[dest_frame] = detector[i];
		proj_array[current_timeslot]->angular_view[dest_frame] = angular_view[i];
		proj_array[current_timeslot]->frame_time[dest_frame] = frame_time[i];

		// copy projection data
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
				proj_array[current_timeslot]->data[dest_frame][j][k] = data[i][j][k];

		dest_frame_counter[current_timeslot]++;
	}

	delete [] dest_frame_counter;
	return proj_array;
}

Projection** Projection::CreateOrderedSubsets(int n_sets)
{
	int i,j,k, m, n;
	Projection** proj_array;
	int* list;

	if(n_timeslots > 1)
	{
		cout << "Projections must be sorted by time slot before creating subsets" << endl;
		return NULL;
	}

	if(num_proj%n_sets)
	{
		cout << "Total projections must divide evenly..." << endl;
		return NULL;
	}

	proj_array = new Projection*[n_sets];

	list = new int[n_sets];
	for(i=0;i<n_sets;i++)
		list[i] = i;
	resort(list,n_sets);
	// copy data into new subsets;
	for(i=0;i<n_sets;i++)
	{
		proj_array[i] = new Projection(num_proj/n_sets, ny, nz, res, aperture);

		// copy detector information
		for(k=0;k<4;k++)
		{
			proj_array[i]->AxialDetectorOffset[k] = AxialDetectorOffset[k];
			proj_array[i]->TransaxialDetectorOffset[k] = TransaxialDetectorOffset[k];
			proj_array[i]->AxialApertureOffset[k] = AxialApertureOffset[k];
			proj_array[i]->TransaxialApertureOffset[k] = TransaxialApertureOffset[k];
			proj_array[i]->RadiusOfRotation[k] = RadiusOfRotation[k];
		}

		// copy frame information and data
		for(j=0;j<num_proj/n_sets;j++)
		{
			// copy angles
			proj_array[i]->angles[j] = angles[list[i]+n_sets*j];
			proj_array[i]->z_shift[j] = z_shift[list[i]+n_sets*j];
			proj_array[i]->detector[j] = detector[list[i]+n_sets*j];
			proj_array[i]->timeslot[j] = timeslot[list[i]+n_sets*j];
			proj_array[i]->angular_view[j] = angular_view[list[i]+n_sets*j];


			// and copy data
			for(m=0;m<nz;m++)
				for(n=0;n<ny;n++)
					proj_array[i]->data[j][m][n] = data[list[i]+n_sets*j][m][n];
		}
	}

	delete [] list;

	return proj_array;
}

float Projection::GetScanLength()
{
	int i;
	float length = 0.0f;
	for(i=0;i<num_proj;i++)
		if(z_shift[i]>length)
			length = z_shift[i];
	return length;
}

float Projection::GetScanDuration()
{
	int i;
	float duration = 0.0f;
	for(i=0;i<num_proj;i++)
		if(frame_time[i]>duration)
			duration = frame_time[i];
	return duration;
}

// get data from header of 
void Projection::SaveDICOM(const char* filename, bool scale)
{
	int i,j,k;
//	float scale;

	RootDicomObj* proj = new RootDicomObj();
	DicomObj* DO;
	DataElement* DE, *tempDE;

	ofstream fout;
	ifstream header_info;
	
	char temp_str[64];
	char* long_str;
	unsigned short temp_us;
	unsigned short* temp_data;

	time_t _Time;
	tm* timeinfo;
	time(&_Time);
	srand((unsigned int)_Time);
	timeinfo = localtime(&_Time);

	header_info.open("c:\\SPECT\\temp\\proj_header_defaults.txt");
	proj->ReadFromFile(header_info);
	header_info.close();

	// ImageType
	if(n_timeslots==1)
		sprintf(temp_str,"ORIGINAL\\PRIMARY\\TOMO\\EMISSION");
	else
		sprintf(temp_str,"ORIGINAL\\PRIMARY\\GATED TOMO\\EMISSION");
	DE = new DataElement(0x0008,0x0008,"CS",strlen(temp_str),temp_str);
	proj->SetElement(DE);

	// SOP Instance UID*
	sprintf(temp_str,"1.3.6.1.4.1.40450.1.1.%d",rand());
	DE = new DataElement(0x0008,0x0018,"UI",strlen(temp_str),temp_str);
	proj->SetElement(DE);

	// Study Instance UID*
	sprintf(temp_str,"1.3.6.1.4.1.0.%4d%02d%02d.%02d%02d%02d.%d",
		timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, rand());
	DE = new DataElement(0x0020,0x000d,"UI",strlen(temp_str),temp_str);
	proj->SetElement(DE);

	// Series Instance UID*
	sprintf(temp_str,"1.3.6.1.4.1.0.%4d%02d%02d.%02d%02d%02d.%d",
		timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, rand());
	DE = new DataElement(0x0020,0x000e,"UI",strlen(temp_str),temp_str);
	proj->SetElement(DE);

	// Patient ID*
	sprintf(temp_str,"%d",rand());
	DE = new DataElement(0x0010,0x0020,"LO",strlen(temp_str),temp_str);
	proj->SetElement(DE);

	// Number of Frames
	sprintf(temp_str,"%d", num_proj);
	DE = new DataElement(0x0028,0x0008,"IS",strlen(temp_str),temp_str);
	proj->SetElement(DE);

	// Rows & Columns
	temp_us = nz;
	DE = new DataElement(0x0028,0x0010,"US",sizeof(temp_us),&temp_us);	// Rows
	proj->SetElement(DE);
	temp_us = ny;
	DE = new DataElement(0x0028,0x0011,"US",sizeof(temp_us),&temp_us);	// Columns
	proj->SetElement(DE);

	// Pixel Spacing
	sprintf(temp_str,"%f\\%f",res,res);
	DE = new DataElement(0x0028,0x0030,"DS",strlen(temp_str),temp_str);
	proj->SetElement(DE);


	// Frame Increment Pointer
	if(n_timeslots > 1)
	{
		unsigned short temp_at[] = {0x0054,0x0010,0x0054,0x0020,0x0054,0x0050,0x0054,0x0060,0x0054,0x0070,0x0054,0x0090};
		DE = new DataElement(0x0028,0x0009,"AT",12*sizeof(unsigned short),temp_at);
	} else {
		unsigned short temp_at[] = {0x0054,0x0010,0x0054,0x0020,0x0054,0x0050,0x0054,0x0090};
		DE = new DataElement(0x0028,0x0009,"AT",8*sizeof(unsigned short),temp_at);
	}	
	proj->SetElement(DE);

	
	// Energy window description (0054, 0010)
	temp_data = new unsigned short[num_proj];
	for(i=0;i<num_proj;i++)
		temp_data[i] = 1;		// only one energy window
	DE = new DataElement(0x0054,0x0010,"US",num_proj*sizeof(unsigned short),temp_data);
	proj->SetElement(DE);

	temp_us = 1;
	DE = new DataElement(0x0054,0x0011,"US",sizeof(temp_us),&temp_us);
	proj->SetElement(DE);
	// 0054,0012 containing the isotope information is stored in default header text file

	// Detector information (0054,0020)
	for(i=0;i<num_proj;i++)
		temp_data[i] = detector[i];
	DE = new DataElement(0x0054,0x0020,"US",num_proj*sizeof(unsigned short),temp_data);
	proj->SetElement(DE);

	temp_us = 4;
	DE = new DataElement(0x0054,0x0021,"US",sizeof(temp_us),&temp_us);
	proj->SetElement(DE);

	DE = new DataElement(0x0054,0x0022,"SQ",0,NULL);
	proj->SetElement(DE);
	for(i=0;i<4;i++)
	{
		// Create Object
		DO = new DicomObj();
		DE->SetObject(DO);

		tempDE = new DataElement(0x0009,0x1031,"DS",2,"0");
		DO->SetElement(tempDE);
		tempDE = new DataElement(0x0009,0x1032,"SH",4,"APT2");
		DO->SetElement(tempDE);
		tempDE = new DataElement(0x0009,0x107d,"FD",sizeof(double),TransaxialApertureOffset+i);
		DO->SetElement(tempDE);
		tempDE = new DataElement(0x0009,0x107e,"FD",sizeof(double),AxialApertureOffset+i);
		DO->SetElement(tempDE);
		tempDE = new DataElement(0x0009,0x107f,"FD",sizeof(double),TransaxialDetectorOffset+i);
		DO->SetElement(tempDE);
		tempDE = new DataElement(0x0009,0x1080,"FD",sizeof(double),AxialApertureOffset+i);
		DO->SetElement(tempDE);
		tempDE = new DataElement(0x0009,0x1081,"FD",sizeof(double),RadiusOfRotation+i);
		DO->SetElement(tempDE);

		j=0;
		while(detector[j]-1 != i) // searching for frame with detector==(i+1) (and angular_view==1))
			j++;
		sprintf(temp_str,"%f",270-angles[j]);
		tempDE = new DataElement(0x0054,0x0200,"DS",strlen(temp_str),temp_str);
		DO->SetElement(tempDE);
	}

	// Rotation information (0054,0050)
	for(i=0;i<num_proj;i++)
		temp_data[i] = 1;
	DE = new DataElement(0x0054,0x0050,"US",num_proj*sizeof(unsigned short),temp_data);
	proj->SetElement(DE);

	temp_us = 1;
	DE = new DataElement(0x0054,0x0051,"US",sizeof(temp_us),&temp_us);
	proj->SetElement(DE);

	DE = new DataElement(0x0054,0x0052,"SQ",0,NULL);
	proj->SetElement(DE);
	DO = new DicomObj();
	DE->SetObject(DO);

	long_str = new char[20*num_proj];
	j=0;
	for(i=0;i<num_proj;i++)
		j += sprintf(long_str+j,"%3.6f\\",z_shift[0]-z_shift[i]);
	long_str[j-1]=0;
	DE = new DataElement(0x0009,0x1024,"DS",strlen(long_str),long_str);
	DO->SetElement(DE);

	j=0;
	for(i=0;i<num_proj;i++)
		j += sprintf(long_str+j,"%04d%02d%02d%02d%02d%02d\\",
			timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
	long_str[j-1]=0;
	DE = new DataElement(0x0009,0x1027,"DT",strlen(long_str),long_str);
	DO->SetElement(DE);
	delete [] long_str;

	sprintf(temp_str,"%1.1f",angles[1]-angles[0]);
	DE = new DataElement(0x0018,0x1144,"DS",strlen(temp_str),temp_str);
	DO->SetElement(DE);

	// RR Interval Vector (0x0054,0x0060)
	if(n_timeslots > 1)
	{
		for(i=0;i<num_proj;i++)
			temp_data[i] = 1;
		DE = new DataElement(0x0054,0x0060,"US",num_proj*sizeof(unsigned short),temp_data);
		proj->SetElement(DE);
	
		temp_us = 1;
		DE = new DataElement(0x0054,0x0061,"US",sizeof(temp_us),&temp_us);
		proj->SetElement(DE);
	}


	// Timeslot vector (0x0054,0x0070)
	if(n_timeslots > 1)
	{
		DE = new DataElement(0x0054,0x0070,"US",num_proj*sizeof(unsigned short),timeslot);
		proj->SetElement(DE);

		temp_us = n_timeslots;
		DE = new DataElement(0x0054,0x0071,"US",sizeof(temp_us),&temp_us);
		proj->SetElement(DE);
	}
	

	// Angular View Vector (0054,0090)
	DE = new DataElement(0x0054,0x0090,"US",num_proj*sizeof(unsigned short),angular_view);
	proj->SetElement(DE);


	// data
	float scale_factor = 1.0f;
	if(scale)
	{
		scale_factor = 0.0f;
		for(i=0;i<num_proj;i++)
			for(j=0;j<nz;j++)
				for(k=0;k<ny;k++)
					if (data[i][j][k] > scale_factor)
						scale_factor = data[i][j][k];

		scale_factor = 32000.0f/scale_factor;

		sprintf(temp_str,"0.0");
		DE = new DataElement(0x0028,0x1052,"DS",strlen(temp_str),temp_str);	// Rescale Intercept
		proj->SetElement(DE);
		sprintf(temp_str,"%1.6f",1/scale_factor);
		DE = new DataElement(0x0028,0x1053,"DS",strlen(temp_str),temp_str); // Rescale Slope
		proj->SetElement(DE);

	}

	delete [] temp_data;
	temp_data = new unsigned short[num_proj*nz*ny];	
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			for(k=0;k<ny;k++)
				temp_data[i*nz*ny+j*ny+k] = (unsigned short)(scale_factor*data[i][j][k]);
	DE = new DataElement(0x7fe0,0x0010,"OW",num_proj*nz*ny*sizeof(unsigned short),temp_data);
	proj->SetElement(DE);
	delete [] temp_data;

	fout.open(filename,ofstream::binary);
	proj->Write(fout);
	fout.close();

	delete proj;

}


/********************************************************************************************
/ Reconstruction
/********************************************************************************************/

float** Reconstruction::sp_matrix = NULL;
bool Reconstruction::rotation_matrix_initialized = false;
ProjTable** Reconstruction::proj_table = NULL;
int Reconstruction::n_proj_tables = 0;



// constructor
Reconstruction::Reconstruction(int size_xy, int size_z, int timeslots, float new_res, char** filename)
{
	int t,i,j;
	ifstream f_in;
	ifstream::pos_type size;
	int slices, first_slice;

	nxy = size_xy;
	nz = size_z;
	res = new_res;

	n_timeslots = timeslots;

	unsigned int rad_sq = (nxy>>1) - 2; // 2 pixel margin allows 5 pixel gaussian matrix (usually use 3 anyways)
	rad_sq *= rad_sq;
	float xy_offset;

	try
	{
		_data = new float[n_timeslots*nz*nxy*nxy]; // linear data storage;
		data = new float***[n_timeslots];
		for(t=0;t<n_timeslots;t++)
			data[t] = new float**[nz];
		for(t=0;t<n_timeslots;t++)
			for(i=0;i<nz;i++)
				data[t][i] = new float*[nxy];
		for(t=0;t<n_timeslots;t++)
			for(i=0;i<nz;i++)
				for(j=0;j<nxy;j++)
					data[t][i][j] = _data + (t*nz*nxy*nxy) + (nxy*nxy*i) + nxy*j;
	
		_data2 = new float[nxy*nxy*nxy];
		data2 = new float**[nxy];
		for(i=0;i<nxy;i++)
			data2[i] = new float*[nxy];
		for(i=0;i<nxy;i++)
			for(j=0;j<nxy;j++)
				data2[i][j] = _data2 + nxy*nxy*i + nxy*j;
	}
	catch(bad_alloc& ba)
	{
		cout << "Reconstruction(int,int,float,char*):bad_alloc caught: " << ba.what() << endl;
	}

	// create rotation mask
	rot_mask = new float*[nxy];
	for(i=0;i<nxy;i++)
		rot_mask[i] = new float[nxy];
	
	xy_offset = (nxy-1)/2.0f;
	for(i=0;i<nxy;i++)
		for(j=0;j<nxy;j++)
			if( ((i-xy_offset)*(i-xy_offset)+(j-xy_offset)*(j-xy_offset)) < rad_sq)
				rot_mask[i][j] = 1;
			else
				rot_mask[i][j] = 0;

	for(t=0;t<n_timeslots;t++)
	{
		Clear(t);
		if(filename)
		{	
			f_in.open(filename[t],ios::in|ios::binary|ios::ate);
			if(f_in.is_open())
			{
				size = f_in.tellg();
				f_in.seekg(0,ios::beg);

				if( (size % (nxy*nxy*sizeof(float))) != 0)
				{
					cout << "Data file doesn't contain an integral number of slices." << endl;
				}
				else
				{
					slices = int(size)/(nxy*nxy*sizeof(float));

					if(slices>nz)
						cout << "Data file exceeds size of reconstruction." << endl;
					else
					{
						first_slice = int((nz-slices)/2);
						for(i=0;i<slices;i++)
							for(j=0;j<nxy;j++) 
								f_in.read(reinterpret_cast<char*>(data[t][first_slice + i][j]),nxy*sizeof(float));
					}
				}
				f_in.close();
			}
		}
	}

}


Reconstruction::Reconstruction(const char* filename) // only loads a single timeframe
{
	RootDicomObj *DCM = new RootDicomObj(filename);
	unsigned short* temp_data;
	int len;
	int i,j,k;
	unsigned short file_slices, file_nxy;
	float file_res;
	char temp_ch[16];

	DCM->GetValue(0x0028,0x0008,&file_slices,sizeof(unsigned short));
	DCM->GetValue(0x0028,0x0010,&file_nxy,sizeof(unsigned short));
	DCM->GetValue(0x0028,0x0030,&temp_ch,16);

	file_res = float(atof(temp_ch));

	Reconstruction(file_nxy,file_slices,1,file_res);

	len = file_slices * file_nxy * file_nxy;
	temp_data = new unsigned short[len];
	DCM->GetValue(0x7FE0,0x0010,temp_data,len*sizeof(unsigned short));
	for(i=0;i<nz;i++)
		for(j=0;j<nxy;j++)
			for(k=0;k<nxy;k++)
				data[0][i][j][k] = temp_data[i*nz*nxy+j*nxy+k];
	delete [] temp_data;


}


Reconstruction::~Reconstruction()
{
	int i, t;

	for(t=0;t<n_timeslots; t++)
		for(i=0;i<nz;i++)
			delete [] data[t][i];
	for(t=0;t<n_timeslots; t++)
		delete data[t];
	for(i=0;i<nxy;i++)
		delete [] data2[i];

	delete [] data;
	delete [] data2;

	delete [] _data;
	delete [] _data2;

	for(i=0;i<nxy;i++)
		delete [] rot_mask[i];
	delete [] rot_mask;

}

void Reconstruction::ReconstructionInfo()
{
	cout << "**********Reconstruction Info***********" << endl;
	cout << nz << "x" << nxy << "x" << nxy << " (" << res << " mm)" << endl;
	cout << "Time slots: " << n_timeslots << endl;
	cout << "****************************************" << endl << endl;

}

// Rotates and copies the data starting at z_offset into data2
// rotation angle is positive for counterclockwise rotation, looking into the bore
// z_offset, xshift, and yshift are specified in mm
// xshift and yshift are applied to data after rotation
// data2 is cleared first
void Reconstruction::RotateOut(float angle, int timeslot, float z_offset, float xshift, float yshift)
{
	int i,j,k,m,n;	// counters (i,j,k for the voxel indices)

	int x, y;		
	int n_grid = (ROT_MAT_SIZE-1)/2;
	int j_int, k_int, j_fr, k_fr;
	float x_center, y_center;		// midpoint of (0,0) pixel in world coordinates 
	x_center = (nxy-1)/2.0f;
	y_center = (nxy-1)/2.0f;

	int i_min, i_max;

	if(!rotation_matrix_initialized)
		InitRotationMatrix();

	int i_offset =  int(floor(z_offset/res));
	float i_fr = (z_offset/res) - i_offset;

	//cout << "(" << z_offset << ")";

	if (i_offset < 0)
	{
		//cout << endl << "warning: clipping in z-dim" << " (" << z_offset << ")" << endl;
		i_min = -i_offset;
	}
	else
		i_min = 0;

	if(i_offset + nxy >= nz)
	{
		//cout << endl << "warning: clipping in z-dim" << " (" << z_offset << ")" << endl;
		i_max = nz - i_offset - 1;
	}
	else
		i_max = nxy;

	float scale;

	// precalculate the sine and cosine of the angle
	float sin_th = sin(angle*atan(1.0f)/45.0f);
	float cos_th = cos(angle*atan(1.0f)/45.0f);

	// clear data2
	Clear2();

	for(j=0;j<nxy;j++)
		for(k=0;k<nxy;k++)
			if(rot_mask[j][k])
			{
				x = int(ROT_SUBPIX * ( (k-x_center)*cos_th - (j-y_center)*sin_th + x_center + xshift/res) );
				y = int(ROT_SUBPIX * ( (k-x_center)*sin_th + (j-y_center)*cos_th + y_center + yshift/res) );

				k_int = x/ROT_SUBPIX;	// integral location
				if(k_int<n_grid)
					continue;
				if(k_int>(nxy-1-n_grid))
					continue;

				j_int = y/ROT_SUBPIX;
				if(j_int<n_grid)
					continue;
				if(j_int>(nxy-1-n_grid))
					continue;

				k_fr = x%ROT_SUBPIX;	// fractional location
				j_fr = y%ROT_SUBPIX;

				for(m=-n_grid;m<=n_grid;m++)
					for(n=-n_grid;n<=n_grid;n++)
					{
						scale = sp_matrix[j_fr][m] * sp_matrix[k_fr][n];
						for(i=i_min;i<i_max;i++)
							data2[i][j_int+m][k_int+n] += scale * ((1-i_fr)*data[timeslot][i_offset+i][j][k] + i_fr*data[timeslot][i_offset+i+1][j][k]);
					}
			}
}

// Rotates and adds the data in data2 to data starting at z_offset
// rotation angle is positive for counterclockwise rotation, looking into the bore
// z_offset, xshift, and yshift are specified in mm
// xshift and yshift are applied to data2 before rotation
void Reconstruction::RotateIn(float angle, int timeslot, float z_offset, float xshift, float yshift)
{
	int i,j,k,m,n;	// counters (i,j,k for the voxel indices)

	int x, y;		
	int n_grid = (ROT_MAT_SIZE-1)/2;
	int j_int, k_int, j_fr, k_fr;
	float x_offset, y_offset;		// midpoint of (0,0) pixel in world coordinates 
	x_offset = (nxy-1)/2.0f - xshift/res;
	y_offset = (nxy-1)/2.0f - yshift/res;

	int i_min, i_max;

	int i_offset = int(floor(z_offset/res));
	float i_fr = (z_offset/res) - i_offset;

	//cout << "(" << z_offset << ")";

	if (i_offset < 0)
	{
		//cout << endl << "warning: clipping in z-dim" << " (" << z_offset << ")" << endl;
		i_min = -i_offset;
	}
	else
		i_min = 0;

	if(i_offset + nxy >= nz)
	{
		//cout << endl << "warning: clipping in z-dim" << " (" << z_offset << ")" << endl;
		i_max = nz - i_offset - 1;
	}
	else
		i_max = nxy;

	float scale;

	if(!rotation_matrix_initialized)
		InitRotationMatrix();

	// precalculate the sine and cosine of the angle
	float sin_th = sin(angle*atan(1.0f)/45.0f);
	float cos_th = cos(angle*atan(1.0f)/45.0f);

	for(j=0;j<nxy;j++)
		for(k=0;k<nxy;k++)
			if(rot_mask[j][k])
			{
				x = int(ROT_SUBPIX * ( (k-x_offset)*cos_th - (j-y_offset)*sin_th + ((nxy-1)/2.0f) ) );
				y = int(ROT_SUBPIX * ( (k-x_offset)*sin_th + (j-y_offset)*cos_th + ((nxy-1)/2.0f) ) );

				k_int = x/ROT_SUBPIX;	// integral location
				if(k_int<n_grid)
					continue;
				if(k_int>(nxy-1-n_grid))
					continue;

				j_int = y/ROT_SUBPIX;
				if(j_int<n_grid)
					continue;
				if(j_int>(nxy-1-n_grid))
					continue;

				k_fr = x%ROT_SUBPIX;	// fractional location
				j_fr = y%ROT_SUBPIX;

				for(m=-n_grid;m<=n_grid;m++)
					for(n=-n_grid;n<=n_grid;n++)
					{
						scale = sp_matrix[j_fr][m] * sp_matrix[k_fr][n];
						for(i=i_min;i<i_max;i++)
						{
							data[timeslot][i_offset+i][j_int+m][k_int+n] += scale * (1-i_fr) * data2[i][j][k];
							data[timeslot][i_offset+i+1][j_int+m][k_int+n] += scale * i_fr * data2[i][j][k];
						}
					}
			}
	}



void Reconstruction::Update(const Reconstruction* scale, const Reconstruction* norm, const Reconstruction* MAPPartial, int timeslot, float threshold)
{
	int i, j, k;

	if(MAPPartial)
	{
		for(i=0;i<nz;i++)
			for(j=0;j<nxy;j++)
				for(k=0;k<nxy;k++)
					if((norm->data[0][i][j][k] + MAPPartial->data[0][i][j][k]) > threshold)	// don't update pixels which are only marginally 'sampled'
							data[timeslot][i][j][k] = data[timeslot][i][j][k] * (scale->data[0][i][j][k]) / (norm->data[0][i][j][k] + MAPPartial->data[0][i][j][k]);
	}
	else
	{
		for(i=0;i<nz;i++)
			for(j=0;j<nxy;j++)
				for(k=0;k<nxy;k++)
					if(norm->data[0][i][j][k] > threshold)
						data[timeslot][i][j][k] = data[timeslot][i][j][k] * (scale->data[0][i][j][k]) / (norm->data[0][i][j][k]);
	}
}

void Reconstruction::UpdateMAPSpatialPartial(const Reconstruction* param, int timeslot, float beta)
{
	int i, j, k;
	int ii, jj, kk;

	float U_plus, U_minus;
//	float w;
	float diff;

	for(i=1;i<nz-1;i++)
		for(j=1;j<nxy-1;j++)
			for(k=1;k<nxy-1;k++)
			{
				U_plus = 0.0f;
				U_minus = 0.0f;
				for(ii=-1;ii<=1;ii++)
					for(jj=-1;jj<=1;jj++)
						for(kk=-1;kk<=1;kk++)
						{
							diff = param->data[timeslot][i][j][k]-param->data[timeslot][i+ii][j+jj][k+kk];
							U_plus += (diff+0.5f) * (diff+0.5f);
							U_minus += (diff-0.5f) * (diff-0.5f);
						}
				data[0][i][j][k] += beta * 0.5f * (U_plus-U_minus);
			}
}

/*
void Reconstruction::UpdateMAPSpatialPartial(const Reconstruction* param, int timeslot, float beta)
{
	int i, j, k;
	int ii, jj, kk;

	float U_plus, U_minus;
	float w;
	float diff;
	float delta = 0;	// 

	for(i=1;i<nz-1;i++)
		for(j=1;j<nxy-1;j++)
			for(k=1;k<nxy-1;k++)
			{
				U_plus = 0.0f;
				U_minus = 0.0f;
				for(ii=-1;ii<=1;ii++)
					for(jj=-1;jj<=1;jj++)
						for(kk=-1;kk<=1;kk++)
						{
							w = sqrt(float(ii*ii + jj*jj + kk*kk));
							if (w==0.0f)
								continue;
							w = 1/w;
							diff = param->data[timeslot][i][j][k]-param->data[timeslot][i+ii][j+jj][k+kk];
							U_plus += w * phi(diff+0.5f, delta);
							U_minus += w * phi(diff-0.5f, delta);
						}
				data[0][i][j][k] += beta*(U_plus-U_minus);
			}
}

float Reconstruction::phi(float u, float delta)
{
	float temp;

	temp = (u*u)/(delta*delta+u*u);
	return temp;

}


void Reconstruction::UpdateMAPCTPartial(const Reconstruction* param, const Reconstruction* CT, int timeslot, float beta, float delta)
{
	int i, j, k;
	int ii, jj, kk;

	float U_plus, U_minus;
	float w;
	float CT_diff;	//difference in CT values

	for(i=1;i<nz-1;i++)
		for(j=1;j<nxy-1;j++)
			for(k=1;k<nxy-1;k++)
			{
				U_plus = 0.0f;
				U_minus = 0.0f;
				for(ii=-1;ii<=1;ii++)
					for(jj=-1;jj<=1;jj++)
						for(kk=-1;kk<=1;kk++)
						{
							CT_diff = 20*(CT->data[0][i][j][k] - CT->data[0][i+ii][j+jj][k+kk]);
							w = sqrt(float(ii*ii + jj*jj + kk*kk));
							if (w==0.0f)
								continue;
							w = 1/(w*(CT_diff*CT_diff+1));
							U_plus += w * phi(param->data[timeslot][i][j][k]+1.0f-param->data[timeslot][i+ii][j+jj][k+kk], delta);
							U_minus += w * phi(param->data[timeslot][i][j][k]-1.0f-param->data[timeslot][i+ii][j+jj][k+kk], delta);
						}
				data[0][i][j][k] += beta*(U_plus-U_minus)/2;
			}
}
*/



void Reconstruction::UpdateMAPTemporalPartial(const Reconstruction* param, int timeslot, float beta)
{
/*	int i, j, k,t;

//	float U_plus, U_minus;

	int n_t = param->n_timeslots;

	float *values = new float[n_t];
//	float coeffs[3];

	for(i=1;i<nz;i++)
		for(j=1;j<nxy;j++)
			for(k=1;k<nxy;k++)
			{

			}

	delete [] values;
	*/
}

/*
void Reconstruction::UpdateGatedMAPPartial(const Reconstruction* param, int timeslot, float beta, float delta)
{
	int i, j, k;
	int ii, jj, kk, tt;

	int t;

	float U_plus, U_minus, diff;
	float w;

	for(i=1;i<nz-1;i++)
		for(j=1;j<nxy-1;j++)
			for(k=1;k<nxy-1;k++)
			{
				U_plus = 0.0f;
				U_minus = 0.0f;
				for(ii=-1;ii<=1;ii++)
					for(jj=-1;jj<=1;jj++)
						for(kk=-1;kk<=1;kk++)
							for(tt=-1;tt<=1;tt+=2)	// compare to adjacent time frames only
							{
								w = sqrt(float(ii*ii + jj*jj + kk*kk + tt*tt));
								if (w==0.0f)
									continue;
								w = 1/w;
								t = (timeslot+tt + param->n_timeslots) % (param->n_timeslots);
								//cout << t;
								diff = param->data[timeslot][i][j][k]+1.0f-param->data[t][i+ii][j+jj][k+kk];
								U_plus += w * phi(diff+0.5f, delta);
								U_minus += w * phi(diff-0.5f, delta);
							}
				data[0][i][j][k] += beta*(U_plus-U_minus);
			}
}*/

void Reconstruction::Clear(int timeslot)
{
	int i, j, k;

	for(i=0;i<nz;i++)
		for(j=0;j<nxy;j++)
			for(k=0;k<nxy;k++)
				data[timeslot][i][j][k] = 0.0;
}

void Reconstruction::Clear2()
{
	int i, j, k;

	for(i=0;i<nxy;i++)
		for(j=0;j<nxy;j++)
			for(k=0;k<nxy;k++)
				data2[i][j][k] = 0.0;
}

void Reconstruction::CopyGate(int dest_timeslot, int src_timeslot)
{
	int i,j;

	if (dest_timeslot==src_timeslot)
		return;
	if (dest_timeslot>=n_timeslots)
		return;
	if (src_timeslot>=n_timeslots)
		return;

	for(i=0;i<nz;i++)
		for(j=0;j<nxy;j++)
			memcpy(data[dest_timeslot][i][j],data[src_timeslot][i][j],nxy*sizeof(float));

}

void Reconstruction::Filter(float sigma)
{
	int i,j,k,t;
	int len_xy, len_z;	// length of padded data
	float ***data_r, ***data_i;

	float *filt_xy, *filt_z;

	float a;

	// pad the data to the nearest power of 2 for the FFT to work
	len_xy=1;
	len_z=1;
	while(len_xy<nxy)
		len_xy = len_xy<<1;
	while(len_z<nz)
		len_z = len_z<<1;

	data_r = new float**[len_z];
	data_i = new float**[len_z];
	for(i=0;i<len_z;i++)
	{
		data_r[i] = new float*[len_xy];
		data_i[i] = new float*[len_xy];
	}
	for(i=0;i<len_z;i++)
		for(j=0;j<len_xy;j++)
		{
			data_r[i][j] = new float[len_xy];
			data_i[i][j] = new float[len_xy];
		}

	// create the filters in fourier space
	filt_xy = new float[len_xy];
	filt_z = new float[len_z];

	a = float(2*M_PI*M_PI*sigma*sigma/(len_xy*len_xy*res*res));
	for(i=0;i<=len_xy/2;i++)
		filt_xy[i] = exp(-a*(i*i));
	for(;i<len_xy;i++)
		filt_xy[i] = filt_xy[len_xy - i];

	a = float(2*M_PI*M_PI*sigma*sigma/(len_z*len_z*res*res));
	for(i=0;i<=len_z/2;i++)
		filt_z[i] = exp(-a*(i*i));
	for(;i<len_z;i++)
		filt_z[i] = filt_z[len_z - i];

	/*
	// TEST CODE
	ofstream f;
	f.open("c:\\SPECT\\filter.txt");
	f << "XY Filter (" << len_xy << ")" << endl;
	for(i=0;i<len_xy;i++)
		f << filt_xy[i] << ", ";
	f << endl << endl;
	f << "Z Filter (" << len_z << ")" << endl;
	for(i=0;i<len_xy;i++)
		f << filt_xy[i] << ", ";
	f << endl;
	f.close();
	// END TEST CODE
	*/
	
	for(t=0;t<n_timeslots;t++)
	{
		// put the data in the complex data buffer
		// buffers are preset to zero so data will be padded with zeros
		for(i=0;i<len_z;i++)
			for(j=0;j<len_xy;j++)
			{
				memset(data_r[i][j],0,len_xy*sizeof(float));
				memset(data_i[i][j],0,len_xy*sizeof(float));
			}

		for(i=0;i<nz;i++)
			for(j=0;j<nxy;j++)
				memcpy(data_r[i][j],data[t][i][j],nxy*sizeof(float));
		
		// fourier transform the data
		fft3(data_r,data_i,len_xy,len_xy,len_z,1);

		// apply the filter in fourier space
		for (i=0;i<len_z;i++)
			for(j=0;j<len_xy; j++)
				for(k=0;k<len_xy;k++)
				{
					a = filt_z[i]*filt_xy[j]*filt_xy[k];
					data_r[i][j][k] *= a;
					data_i[i][j][k] *= a;
				}

		// inverse fourier transform
		fft3(data_r,data_i,len_xy,len_xy,len_z,-1);

		// put data back where it belongs
		for(i=0;i<nz;i++)
			for(j=0;j<nxy;j++)
				memcpy(data[t][i][j],data_r[i][j],nxy*sizeof(float));

		// clean up stray negatives caused by the filter
		for(i=0;i<nz;i++)
			for(j=0;j<nxy;j++)
				for(k=0;k<nxy;k++)
					if(data[t][i][j][k]<0)
						data[t][i][j][k] = -data[t][i][j][k];
	}

	// free allocated memory
	for(i=0;i<len_z;i++)
		for(j=0;j<len_xy;j++)
		{
			delete [] data_r[i][j];
			delete [] data_i[i][j];
		}
	for(i=0;i<len_z;i++)
	{
		delete [] data_r[i];
		delete [] data_i[i];
	}
	delete [] data_r;
	delete [] data_i;

	delete [] filt_xy;
	delete [] filt_z;

}


//void Reconstruction::CleanUpStrayPixels(Reconstruction* param, float limit)
void Reconstruction::CleanUpStrayPixels(float limit)
{
	int i,j,k,t;

//	if(limit > 0.0f)
//	{
		for(t=0;t<n_timeslots;t++)
			for(i=0;i<nz;i++)
				for(j=0;j<nxy;j++)
					for(k=0;k<nxy;k++)
						if(data[t][i][j][k] > limit)
							data[t][i][j][k] = limit;
						//if(param->data[0][i][j][k] == 0.0f)
						//	data[t][i][j][k] = 0.0f;
//	}
	
	/*else
	{
		for(t=0;t<n_timeslots;t++)
			for(i=0;i<nz;i++)
				for(j=0;j<nxy;j++)
					for(k=0;k<nxy;k++)
						if(param->data[0][i][j][k] == 0.0f)
					}
							data[t][i][j][k] = 0.0f;
	}*/

}

void Reconstruction::Mask(int radius)
{
	int i,j,k,t;

	int rad_sq = radius*radius;
	float xy_offset = (nxy-1)/2.0f;

	for(j=0;j<nxy;j++)
		for(k=0;k<nxy;k++)
			if( ((j-xy_offset)*(j-xy_offset)+(k-xy_offset)*(k-xy_offset)) > rad_sq)
			{
				for(t=0;t<n_timeslots;t++)
					for(i=0;i<nz;i++)
						data[t][i][j][k] = 0.0f;
			}
}

/*
void Reconstruction::ClearLowPixels(float limit)
{
	int i,j,k,t;

	for(t=0;t<n_timeslots;t++)
		for(i=0;i<nz;i++)
			for(j=0;j<nxy;j++)
				for(k=0;k<nxy;k++)
					if(data[t][i][j][k] < limit)
						data[t][i][j][k] = 0.0f;
} */

void Reconstruction::SetPlane(int k)
{
	int i,j;
	for(i=0;i<nz;i++)
		for(j=0;j<nxy;j++)
			data[0][i][j][k] = 1.0;
}

void Reconstruction::CylindricalPhantom(float radius, int first_slice, int last_slice, float offset)
{
	int i,j,k;
	float r_sq;

	float ctr = (nxy-1)/2.0f;
	
	Clear(0);

	if (first_slice<0)
		first_slice = 0;
	if (last_slice>(nz-1))
		last_slice = nz-1;

	float r = radius*radius;
		
	for(j=0;j<nxy;j++)
		for(k=0;k<nxy;k++)
		{
			r_sq = ((j*1.0f-ctr-offset)*(j*1.0f-ctr-offset)+(k*1.0f-ctr)*(k*1.0f-ctr));
			if( r_sq < r )
				for(i=first_slice; i<last_slice; i++)
					data[0][i][j][k] = 1.0f;
		}

}

void Reconstruction::Crosshairs()
{
	int i, j;

	for(i=0;i<nz;i++)
	{
		for(j=0;j<nxy;j++)
		{
			data[0][i][j][nxy/2] = 1.0f;
			data[0][i][j][nxy/2-1] = 1.0f;
		}

		for(j=0;j<nxy;j++)
		{
			data[0][i][nxy/2][j] = 1.0f;
			data[0][i][nxy/2-1][j] = 1.0f;
		}
	}


}

void Reconstruction::HotRodPhantom()
{
	float xp,yp, x, y;	// center of rod
	int i, j, k;
	int rod, row, m;
	float sin_th, cos_th;
	float dist;

	int sub_j, sub_k;
	float sum;

	float j_ctr = 63.5;
	float k_ctr = 63.5;

	float theta[] = {30, 90, 150, 210, 270, 330};
	float d[] = {1.2f, 1.6f, 2.4f, 3.2f, 4.0f, 4.8f};
	int n_rows[] = {8, 6, 4, 3, 2, 2};
	float offset;

	for(i=0;i<nz;i++)
		for(j=0;j<nxy;j++)
			for(k=0;k<nxy;k++)
				data[0][i][j][k] = 0.0f;

	for(rod=0; rod<6; rod++)
		for(row=0;row < n_rows[rod];row++)
			for(m = 0; m<row+1; m++)
			{
				// locate unrotated center
				offset = 4.8f + sqrt(3.0f)*d[rod]/2;
				xp = offset + sqrt(3.0f) * row * d[rod];
				yp = (2*m - row)*d[rod];

				// rotate point
				sin_th = sin(theta[rod] * atan(1.0f)/45.0f);
				cos_th = cos(theta[rod] * atan(1.0f)/45.0f);

				x = xp * cos_th - yp * sin_th;
				y = xp * sin_th + yp * cos_th;
				
				// put a point at each rod center
				for(j=0;j<nxy;j++)
					for(k=0;k<nxy;k++)
					{
						sum = 0.0f;
						for(sub_j = 0; sub_j<10; sub_j++)
							for(sub_k = 0; sub_k<10; sub_k++)
							{
								dist = (x - (k-k_ctr)*0.6f - (sub_k-4.5f)*0.06f) * (x - (k-k_ctr)*0.6f - (sub_k-4.5f)*0.06f) + ( y-(j-j_ctr)*0.6f - (sub_j-4.5f)*0.06f) * (y - (j-j_ctr)*0.6f - (sub_j-4.5f)*0.06f); 
								if (dist <= (d[rod]*d[rod])/4)
									sum += 0.01f;
							}
				
						for(i=0;i<nz;i++)
							data[0][i][j][k] += sum;
					}
			}

}


void Reconstruction::ColdRodPhantom()
{
	int i,j,k;
	HotRodPhantom();

	for(i=0;i<nz;i++)
		for(j=0;j<nxy;j++)
			for(k=0;k<nxy;k++)
				if(rot_mask[j][k])
					data[0][i][j][k] = 1 - data[0][i][j][k];

}


void Reconstruction::TestRotation()
{
	// Create Phantom
	Clear(0);
	Crosshairs();
	// save original data
	SaveBIN("c:\\SPECT\\temp\\original.bin");

	// RotateOut, with offsets
	RotateOut(30.0f,0,15,5,10);

	Clear(0);
	RotateIn(0.0f,0,15,0,0);
	SaveBIN("c:\\SPECT\\temp\\rotated.bin");

	// RotateIn, with opposite offsets
	Clear(0);
	RotateIn(0.0f,0,15,0,0);
	SaveBIN("c:\\SPECT\\temp\\unrotated.bin");

	Clear(0);
	RotateIn(0.0f,0,15,0,0);
	SaveBIN("c:\\SPECT\\temp\\unrotated2.bin");
}


float Reconstruction::MSE(Reconstruction* param)
{
	int i,j,k;
	float sum;

	float ratio;

	ratio = Activity()/param->Activity();

	sum = 0.0;

	for(i=0;i<nz;i++)
		for(j=0;j<nxy;j++)
			for(k=0;k<nxy;k++)
				sum += ( data[0][i][j][k] - (ratio*param->data[0][i][j][k]) )*( data[0][i][j][k] - (ratio*param->data[0][i][j][k]) );

	sum /= (nxy*nxy*nz);

	return sqrt(sum);
}

float Reconstruction::Activity()
{
	int i,j,k;
	float sum;

	sum = 0.0;

	for(i=0;i<nz;i++)
		for(j=0;j<nxy;j++)
			for(k=0;k<nxy;k++)
				sum += data[0][i][j][k];

	return sum;
}

/*
float Reconstruction::VOIActivity(int timeslot)
{

	int i,j,k;
	float sum = 0.0f;

	for(i=82;i<116;i++)
		for(j=65;j<91;j++)
			for(k=45;k<78;k++)
				sum+= data[timeslot][i][j][k];
	return sum;

}
*/


float Reconstruction::VOIActivity(int timeslot)
{

	int i,j,k;
	float sum = 0.0f;

	for(i=-14;i<15;i++)
		for(j=-14;j<15;j++)
			for(k=-14;k<15;k++)
				sum+= data[timeslot][nz/2+i][nxy/2+j][nxy/2+k];
	return sum;

}


void Reconstruction::SaveBIN(const char* filename, int mask)
{
	int i,j,k,t;
	char _filename[260];
	char* _ptr;
	const char* _ptr2;
	ofstream fout;

	float r;

	// create mask in the first slice of data2
	for(j=0;j<nxy;j++)
		for(k=0;k<nxy;k++)
		{
			r = float((j-(nxy-1)/2)*(j-(nxy-1)/2)+(k-(nxy-1)/2)*(k-(nxy-1)/2));
			if(r <= (mask*mask))
				data2[0][j][k] = 1;
			else
				data2[0][j][k] = 0;
		}

	for(t=0;t<n_timeslots;t++)
	{
		strcpy(_filename,filename);
		
		if(n_timeslots>1)
		{
			_ptr = strrchr(_filename,'.');
			_ptr2 = strrchr(filename,'.');
			if(_ptr)
				sprintf(_ptr,"_%d%s",t+1,_ptr2);
			else
				sprintf(_filename,"%s_%d",_filename,t+1);
		}
		

		fout.open(_filename,ofstream::binary);
		if(fout.is_open())
		{
			for(i=0;i<nz;i++)
			{
				// copy frame to second slice of data2
				for(j=0;j<nxy;j++)
					for(k=0;k<nxy;k++)
						data2[1][j][k] = data2[0][j][k] * data[t][i][j][k];

				// write frame
				for(j=0;j<nxy;j++)
					fout.write(reinterpret_cast<char*>(data2[1][j]),nxy*sizeof(float));
			}
			fout.close();
		}
	}

}

void Reconstruction::LoadBIN(const char* filename, bool duplicate_gate0)
{
	int i,j,k,t;
	ifstream fin;

	//TODO: check for file existence and size compatibility

	fin.open(filename,ofstream::binary);
	for(i=0;i<nz;i++)
		for(j=0;j<nxy;j++)
			fin.read(reinterpret_cast<char*>(data[0][i][j]),nxy*sizeof(float));

	// copy data into all timeslots
	if(duplicate_gate0)
		for(t=1;t<n_timeslots;t++)
			CopyGate(t,0);
	else
		for(t=1;t<n_timeslots;t++)
			for(i=0;i<nz;i++)
				for(j=0;j<nxy;j++)
					fin.read(reinterpret_cast<char*>(data[t][i][j]),nxy*sizeof(float));
	fin.close();

	for(t=1;t<n_timeslots;t++)
		for(i=0;i<nz;i++)
			for(j=0;j<nxy;j++)
				for(k=0;k<nxy;k++)
					data[t][i][j][k] *= rot_mask[j][k];

}


int Reconstruction::LoadProjTables(const Projection *proj, const char* dir)
{
	char filename[260];
	int n_pinholes;
	int ph;
	ifstream fin;

	if(proj_table)
		return 0;

	// load proj_tables (currently only works for APT2)
	if(!strncmp(proj->aperture,"APT2",5))
		n_pinholes = 9;
	else
		return -1;	// error

	proj_table = new ProjTable*[n_pinholes];
	for(ph = 0; ph < n_pinholes; ph++)
	{
		sprintf(filename, "%s" __SL "%s_%d_%d_%3.1f_%d_%d_%3.1f_%d.pta",dir,proj->aperture,nxy,nxy,res,proj->ny,proj->nz,proj->res,ph);
		//cout << "Loading " << filename << endl;
		try
		{
			proj_table[n_proj_tables] = new ProjTable(filename);
		}
		catch (bad_alloc &ba)
		{
			cout << "bad_alloc caught: " << ba.what() << endl;
			return -1;
		}
		if(proj_table[n_proj_tables])
			n_proj_tables++;
		else
			return -1;
			
	}

	return 0;

}

void Reconstruction::UnloadProjTables()
{
	while(n_proj_tables)
	{
		delete proj_table[n_proj_tables-1];
		n_proj_tables--;
	}

	delete [] proj_table;

}

bool Reconstruction::ProjTableLoaded()
{
	if(proj_table)
		return true;
	else
		return false;
}

void Reconstruction::Project(Projection* proj, const char* mmp_atten)
{
	int i, j, k;
	int angle;
	float x_shift, y_shift, z_shift;

	AttenMap* atten;

	char filename[260];

	if(!ProjTableLoaded())
		return;

	cout << "Projecting ";

	// clear the projections
	proj->Clear();

	for(angle = 0;angle<proj->num_proj;angle++)	// calling this variable angle is confusing, it's the index of the projection, not the actual projection angle
	{
		//cout << "(" << proj->angles[angle] << ")";
		cout << ".";

		x_shift = (45.0f - float(proj->RadiusOfRotation[proj->detector[angle]-1]));
		y_shift = float(proj->TransaxialApertureOffset[proj->detector[angle]-1]);
		z_shift = float(proj->AxialApertureOffset[proj->detector[angle]-1]);

		//cout << proj->detector[angle] << " : (" << x_shift << ", " << y_shift << ", " << z_shift << ")" << endl;

		// load the attenuation map for this angle
		if(mmp_atten)
		{
			sprintf(filename, "%s" __SL "att_map_%d_%d.bin",mmp_atten,proj->detector[angle],proj->angular_view[angle]);
			atten = new AttenMap(filename,n_proj_tables,nxy,nxy);

		}

		for(int ph=0; ph < n_proj_tables; ph++)
		{

			// rotate the projection data
			RotateOut(-proj->angles[angle], proj->timeslot[angle]-1, proj->z_shift[angle] + z_shift,x_shift, y_shift);

			if(mmp_atten)
			{
				// apply attenuation to rotated data
				for(i=0;i<nxy;i++)
					for(j=0;j<nxy;j++)
						for(k=0;k<nxy;k++)
							data2[i][j][k] *= atten->atten_map[ph][i][j][k];
			}

			// projection matrix
			for(unsigned int i=0;i<proj_table[ph]->n;i++)
				proj->_data[(angle<<16) + (proj_table[ph]->pi[i]<<8) + proj_table[ph]->pj[i]] += proj_table[ph]->weight[i] * data2[proj_table[ph]->ri[i]][proj_table[ph]->rj[i]][proj_table[ph]->rk[i]];
//				proj->_data[(angle<<16) + (proj_table[ph]->pi[i]<<8) + proj_table[ph]->pj[i]] += proj_table[ph]->weight[i] * _data2[ (proj_table[ph]->ri[i] << 14) + (proj_table[ph]->rj[i] << 7) + (proj_table[ph]->rk[i]) ];
				// this hardcodes a 256 x 256 projection matrix
		}
		if(mmp_atten)
			delete atten;

	}

	cout << endl;
	
}

/******************************************
/ Backprojects the projection data into the reconstruction
/ Only backprojects into timeslot 0
/******************************************/
void Reconstruction::Backproject(const Projection* proj, const char* mmp_atten)
{
	int i,j,k;
	int angle;

	uchar src_i, src_j, dest_i, dest_j, dest_k;

	float x_shift, y_shift, z_shift;

	char filename[260];

	AttenMap* atten;

	if(!ProjTableLoaded())
		return;

	cout << "Backprojecting ";
	Clear(0);
	for(angle=0; angle<proj->num_proj; angle++)
	{
		//cout << "(" << proj->angles[angle] << ")";
		cout << ".";

		x_shift = -(45.0f - float(proj->RadiusOfRotation[proj->detector[angle]-1]));
		y_shift = -float(proj->TransaxialApertureOffset[proj->detector[angle]-1]);
		z_shift = float(proj->AxialApertureOffset[proj->detector[angle]-1]);

		// load attenuation map
		if(mmp_atten)
		{
			sprintf(filename, "%s" __SL "att_map_%d_%d.bin",mmp_atten,proj->detector[angle],proj->angular_view[angle]);
			try
			{
				atten = new AttenMap(filename,n_proj_tables,nxy,nxy);
			}
			catch(bad_alloc &ba)
			{
				cout << "bad_alloc caught: " << ba.what() << endl;
				return;
			}

		}

		for(int ph=0;ph<n_proj_tables;ph++)
		{
			Clear2();

			// project
			for(unsigned int s=0;s<proj_table[ph]->n;s++)
			{
				src_i = proj_table[ph]->pi[s];
				src_j = proj_table[ph]->pj[s];
				dest_i = proj_table[ph]->ri[s];
				dest_j = proj_table[ph]->rj[s];
				dest_k = proj_table[ph]->rk[s];
//				_data2[(dest_i << 14) + (dest_j << 7) + dest_k] += proj_table[ph]->weight[s] * proj->_data[(angle << 16) + (src_i << 8) + src_j];	// hardcoded dimensions
				data2[dest_i][dest_j][dest_k] += proj_table[ph]->weight[s] * proj->_data[(angle << 16) + (src_i << 8) + src_j];	// hardcoded dimensions

			}
					
			if(mmp_atten)
			{
				// apply attenuation
				for(i=0;i<nxy;i++)
					for(j=0;j<nxy;j++)
						for(k=0;k<nxy;k++)
							data2[i][j][k] *= atten->atten_map[ph][i][j][k];
			}

			RotateIn(proj->angles[angle], 0, proj->z_shift[angle] + z_shift, x_shift, y_shift);
		}

		if(mmp_atten)
			delete atten;
	}

	cout << endl;

}

// the following routine intializes the subpixel matrix
// for Gaussian rotation and must be called once
// before attempting any rotation or (back)projections
void Reconstruction::InitRotationMatrix()
{
	int i, j;
	int n_grid = (ROT_MAT_SIZE-1)/2;
	float z, sum;
	float sigma = ROT_FWHM/(2.0f*sqrt(2.0f*log(2.0f)));
	
	// allocate memory 
	// table is formatted as sp_table[subpixel offset][relative pixel]
	sp_matrix = new float*[ROT_SUBPIX];
	for(i=0;i<ROT_SUBPIX;i++)
		sp_matrix[i] = new float[n_grid*2+1] + n_grid; // allows us to index the pixels as -2, -1, 0, 1, 2
	// assign the values
	for(i=0;i<ROT_SUBPIX;i++)
	{
		sum = 0.0;
		for(j=-n_grid;j<=n_grid;j++)
		{
			if(j<n_grid)
			{
				z = (1 + j - i/(1.0f*ROT_SUBPIX))/sigma;
				sp_matrix[i][j] = 0.5f*erfc(-z/sqrt(2.0f)) - sum;
				sum += sp_matrix[i][j];
			}
			else
				sp_matrix[i][j] = 1-sum;
		}

	}

	rotation_matrix_initialized = true;
}

// Doesn't handle gated studies yet
void Reconstruction::SaveDICOM(const char* filename, const char* ProjDCM, int r)
{
	int i,j,k, t;
	float scale;

	RootDicomObj* proj = new RootDicomObj(ProjDCM);
	RootDicomObj* recon = new RootDicomObj();
	DicomObj* DO;
	DataElement* DE;

	ofstream fout;
	
	char temp_str[64];
	unsigned short temp_us;
//	unsigned short temp_at[2];
	unsigned short* temp_data;
	unsigned long len;

	time_t _Time;
	tm* timeinfo;
	time(&_Time);
	srand((unsigned int)_Time);
	timeinfo = localtime(&_Time);

	int n_slices, first_slice, last_slice, rows_cols, xy_offset;
	n_slices = nz-100;
	first_slice = 50;
	last_slice = n_slices + first_slice-1;
	rows_cols = 2*(r+5);
	xy_offset = (nxy-rows_cols)/2;

	// save essential fields, patient info, time, date
	// Transfer Syntax
	sprintf(temp_str,"1.2.840.10008.1.2.1");
	DE = new DataElement(0x0002,0x0010,"UI",strlen(temp_str),temp_str);
	recon->SetElement(DE);

	//Add stupid meta tag 
	len = recon->GetLength(0x0002,0x0010) + 8;
	recon->SetElement(new DataElement(0x0002,0x0000,"UL",sizeof(unsigned long),&len));
	
	// Image Type
	if(n_timeslots==1)
		sprintf(temp_str,"ORIGINAL\\PRIMARY\\RECON TOMO\\EMISSION");
	else
		sprintf(temp_str,"ORIGINAL\\PRIMARY\\RECON GATED TOMO\\EMISSION");
	DE = new DataElement(0x0008,0x0008,"CS",strlen(temp_str),temp_str);
	recon->SetElement(DE);

	// SOP Class UID
	sprintf(temp_str,"1.2.840.10008.5.1.4.1.1.20");
	DE = new DataElement(0x0008,0x0016,"UI",strlen(temp_str),temp_str);
	recon->SetElement(DE);
	
	// SOP Instance UID
//	sprintf(temp_str,"1.3.6.1.4.1.40450.1.1.%d",rand());
	sprintf(temp_str,"1.3.6.1.4.1.XXXXXXXXX.%d",rand());
	DE = new DataElement(0x0008,0x0018,"UI",strlen(temp_str),temp_str);
	recon->SetElement(DE);

	// Study UID*
	len = proj->GetValue(0x0020,0x000d,temp_str,sizeof(temp_str));
	DE = new DataElement(0x0020,0x000d,"UI",len,temp_str);
	recon->SetElement(DE);

	// Series UID*
	len = proj->GetValue(0x0020,0x000e,temp_str,sizeof(temp_str));
	DE = new DataElement(0x0020,0x000e,"UI",len,temp_str);
	recon->SetElement(DE);

	// Series Description*
	len = proj->GetValue(0x0008,0x1030,temp_str,sizeof(temp_str));
	DE = new DataElement(0x0008,0x1030,"LO",len,temp_str);
	recon->SetElement(DE);

	// Study Description*
	len = proj->GetValue(0x0008,0x103e,temp_str,sizeof(temp_str));
	DE = new DataElement(0x0008,0x103e,"LO",len,temp_str);
	recon->SetElement(DE);

	// Modality
	DE = new DataElement(0x0008,0x0060,"CS",2,"NM");
	recon->SetElement(DE);

	// Software Version
	sprintf(temp_str,"%d.%d",SPECT4D_MAJOR_VERSION,SPECT4D_MINOR_VERSION);
	recon->SetElement(new DataElement(0x0018,0x1020,"CS",strlen(temp_str),temp_str));

	// Patient Position*
	len = proj->GetValue(0x0018,0x5100,temp_str,sizeof(temp_str));
	DE = new DataElement(0x0018,0x5100,"CS",len,temp_str);
	recon->SetElement(DE);

	// Patient name*
	len = proj->GetValue(0x0010,0x0010,temp_str,sizeof(temp_str));
	DE = new DataElement(0x0010,0x0010,"PN",len,temp_str);
	recon->SetElement(DE);

	// Patient ID*
	len = proj->GetValue(0x0010,0x0020,temp_str,sizeof(temp_str));
	if(len==0)
	{
		//make up a patient ID
		sprintf(temp_str,"000000");
	}
	DE = new DataElement(0x0010,0x0020,"LO",len,temp_str);
	recon->SetElement(DE);

	// Patient birthdate*
	len = proj->GetValue(0x0010,0x0030,temp_str,sizeof(temp_str));
	DE = new DataElement(0x0010,0x0030,"DA",len,temp_str);
	recon->SetElement(DE);

	// Patient sex*
	len = proj->GetValue(0x0010,0x0040,temp_str,sizeof(temp_str));
	DE = new DataElement(0x0010,0x0040,"CS",len,temp_str);
	recon->SetElement(DE);

	// StudyDate*
	len = proj->GetValue(0x0008,0x0020,temp_str,sizeof(temp_str));
	DE = new DataElement(0x0008,0x0020,"DA",len,temp_str);
	recon->SetElement(DE);
	
	// StudyTime*
	len = proj->GetValue(0x0008,0x0030,temp_str,sizeof(temp_str));
	DE = new DataElement(0x0008,0x0030,"TM",len,temp_str);	
	recon->SetElement(DE);

	// Slice thickness
	sprintf(temp_str,"%1.4f",res);
	DE = new DataElement(0x0018,0x0050,"DS",strlen(temp_str),temp_str);
	recon->SetElement(DE);

	// Spacing between slices
	sprintf(temp_str,"%1.4f",-res);
	DE = new DataElement(0x0018,0x0088,"DS",strlen(temp_str),temp_str);
	recon->SetElement(DE);

	// Instance Number
	sprintf(temp_str,"1");
	DE = new DataElement(0x0020,0x0013,"IS",strlen(temp_str),temp_str);
	recon->SetElement(DE);

	// Samples per pixel
	temp_us = 1;
	DE = new DataElement(0x0028,0x0002,"US",sizeof(temp_us),&temp_us);
	recon->SetElement(DE);

	// Photometric Interpretation
	sprintf(temp_str,"MONOCHROME2");
	DE = new DataElement(0x0028,0x0004,"CS",strlen(temp_str),temp_str);
	recon->SetElement(DE);

	// Number of Frames
	sprintf(temp_str,"%d", n_slices*n_timeslots);
	DE = new DataElement(0x0028,0x0008,"IS",strlen(temp_str),temp_str);
	recon->SetElement(DE);

	// Rows & Columns
	temp_us = rows_cols;
	DE = new DataElement(0x0028,0x0010,"US",sizeof(temp_us),&temp_us);	// Rows
	recon->SetElement(DE);
	temp_us = rows_cols;
	DE = new DataElement(0x0028,0x0011,"US",sizeof(temp_us),&temp_us);	// Columns
	recon->SetElement(DE);

	// Pixel Spacing
	sprintf(temp_str,"%f\\%f",res,res);
	DE = new DataElement(0x0028,0x0030,"DS",strlen(temp_str),temp_str); // Pixel Spacing
	recon->SetElement(DE);

	// Bit storage
	temp_us = 16;
	DE = new DataElement(0x0028,0x0100,"US",sizeof(temp_us),&temp_us);	// Bits allocated
	recon->SetElement(DE);
	DE = new DataElement(0x0028,0x0101,"US",sizeof(temp_us),&temp_us);	// Bits stored
	recon->SetElement(DE);

	temp_us = 15;
	DE = new DataElement(0x0028,0x0102,"US",sizeof(temp_us),&temp_us);	// High bit
	recon->SetElement(DE);

	temp_us = 0;
	DE = new DataElement(0x0028,0x0103,"US",sizeof(temp_us),&temp_us);	// Pixel representation (?)
	recon->SetElement(DE);

	temp_us = 0;
	DE = new DataElement(0x0028,0x0106,"US",sizeof(temp_us),&temp_us);	// Smallest value
	recon->SetElement(DE);

	temp_us = 32767;
	DE = new DataElement(0x0028,0x0107,"US",sizeof(temp_us),&temp_us);	// Largest value
	recon->SetElement(DE);
	
	// Energy window description
	temp_us = 1;
	DE = new DataElement(0x0054,0x0011,"US",2,&temp_us);
	recon->SetElement(DE);
	DO = new DicomObj();
	sprintf(temp_str,"Tc-99m");
	DE = new DataElement(0x0054,0x0018,"SH",strlen(temp_str),temp_str);	// Isotope
	DO->SetElement(DE);
	DE = new DataElement(0x0054,0x0012,"SQ",0,NULL);
	DE->SetObject(DO);
	DE->CheckLength();
	recon->SetElement(DE);

	// Detector Information (should maybe be copied from the projection data?)
	temp_us = 1;
	recon->SetElement(new DataElement(0x0054,0x0021,"US",2,&temp_us)); // Number of detectors
	DO = new DicomObj();
	sprintf(temp_str,"APT2");
	DO->SetElement(new DataElement(0x0018,0x1180,"SH",strlen(temp_str),temp_str));	// Aperture
	sprintf(temp_str,"-1\\0\\0\\0\\-1\\0");
	DO->SetElement(new DataElement(0x0020,0x0037,"DS",strlen(temp_str),temp_str));	// Isotope
	DE = new DataElement(0x0054,0x0022,"SQ",0,NULL);
	DE->SetObject(DO);
	DE->CheckLength();
	recon->SetElement(DE);

	temp_data = new unsigned short[n_slices*n_timeslots];
	// Frame Increment Pointer
	if (n_timeslots>1)
	{
		temp_data[0] = 0x0054;
		temp_data[1] = 0x0070;
		temp_data[2] = 0x0054;
		temp_data[3] = 0x0080;
		DE = new DataElement(0x0028,0x0009,"AT",4*sizeof(unsigned short),temp_data);
	}
	else
	{
		temp_data[0] = 0x0054;
		temp_data[1] = 0x0080;
		DE = new DataElement(0x0028,0x0009,"AT",2*sizeof(unsigned short),temp_data);
	}
	recon->SetElement(DE);

	// Timeslot vector
	if(n_timeslots>1)
	{
		for(t=0;t<n_timeslots;t++)
			for(i=0; i<n_slices; i++)
				temp_data[t*n_slices+i] = t+1;
		DE = new DataElement(0x0054,0x0070,"US",n_slices*n_timeslots*sizeof(unsigned short),temp_data);	// Frame number
		recon->SetElement(DE);
		temp_us = n_timeslots;
		DE = new DataElement(0x0054,0x0071,"US",sizeof(unsigned short),&temp_us);
		recon->SetElement(DE);
	}
	// Slice Vector
	for(t=0;t<n_timeslots;t++)
		for(i=0; i<n_slices; i++)
			temp_data[t*n_slices+i] = i+1;
	DE = new DataElement(0x0054,0x0080,"US",n_slices*n_timeslots*sizeof(unsigned short),temp_data);	// Frame number
	recon->SetElement(DE);
	temp_us = n_slices;
	DE = new DataElement(0x0054,0x0081,"US",sizeof(unsigned short),&temp_us);
	recon->SetElement(DE);
	delete [] temp_data;


	// create mask
	int temp;
	Clear2();
	for(j=0;j<rows_cols;j++)
		for(k=0;k<rows_cols;k++)
		{
			temp = ((j-rows_cols/2)*(j-rows_cols/2)+(k-rows_cols/2)*(k-rows_cols/2));
			if (temp<=(r*r))
				data2[0][j][k]=1.0f;
		}


	// find maximum value
	scale = 0.0;
	for(t=0;t<n_timeslots;t++)
		for(i=first_slice;i<last_slice;i++)
			for(j=0;j<rows_cols;j++)
				for(k=0;k<rows_cols;k++)
				{
					if((data[t][i][j+xy_offset][k+xy_offset] * data2[0][j][k]) > scale)
						scale = data[t][i][j+xy_offset][k+xy_offset];
				}

	scale = 32000.0f/scale;

	// write data to dicom
	temp_data = new unsigned short[n_timeslots * n_slices * rows_cols * rows_cols];	
	for(t=0;t<n_timeslots;t++)
		for(i=0;i<n_slices;i++)
			for(j=0;j<rows_cols;j++)
				for(k=0;k<rows_cols;k++)
					temp_data[t*n_slices*rows_cols*rows_cols + i*rows_cols * rows_cols + j*rows_cols + k] = (unsigned short)(scale * data[t][first_slice+i][xy_offset+j][xy_offset+k] * data2[0][j][k]);
	DE = new DataElement(0x7fe0,0x0010,"OW",n_timeslots*n_slices*rows_cols*rows_cols*sizeof(unsigned short),temp_data);
	recon->SetElement(DE);
	delete [] temp_data;

	sprintf(temp_str,"0.0");
	DE = new DataElement(0x0028,0x1052,"DS",strlen(temp_str),temp_str);	// Rescale Intercept
	recon->SetElement(DE);
	sprintf(temp_str,"%1.6f",1/scale);
	DE = new DataElement(0x0028,0x1053,"DS",strlen(temp_str),temp_str); // Rescale Slope
	recon->SetElement(DE);

	fout.open(filename,ofstream::binary);
	recon->Write(fout);
	fout.close();

	delete recon;
	delete proj;

}


#include "pinhole.h"
void Reconstruction::CreateMMPAttenMap(const char* working_dir, char* pinhole_file, Projection* proj, float atten_scale)
{

	int angle, pinhole, i, j, k, m;	// current voxel

	float**** amap;
	float mu;

	char filename[64];
	ifstream fin;
	ofstream fout;
	ofstream f_info;

	float di, dj;
	int i2, j2;

#ifdef _WIN32
	__declspec(align(16)) float vox_pos[4];	// voxel position
	__declspec(align(16)) float temp_fl[4];
#elif defined __unix
	__attribute__((aligned(16))) float vox_pos[4];
	__attribute__((aligned(16))) float temp_fl[4];
#endif

	int n_proj = proj->GetNumProj();

	int n_pinholes;

	Pinhole** pinhole_geo = LoadPinholeGeometry(pinhole_file, 45.0f, n_pinholes);

	// create the attenuation map
	// allocate memory
	amap = new float***[n_pinholes];
	for(pinhole=0;pinhole<n_pinholes;pinhole++)
		amap[pinhole] = new float**[nxy];
	for(pinhole=0;pinhole<n_pinholes;pinhole++)
		for(i=0;i<nxy;i++)
			amap[pinhole][i] = new float*[nxy];
	for(pinhole=0;pinhole<n_pinholes;pinhole++)
		for(i=0;i<nxy;i++)
			for(j=0;j<nxy;j++)
				amap[pinhole][i][j] = new float[nxy];

	// precalculate x, y, z
	float *src_xy, *src_z;
	src_z = new float[nxy];
	for(i=0;i<nxy;i++)
		src_z[i] = (i-(nxy-1)/2.0f)*res;
	src_xy = new float[nxy];
	for(i=0;i<nxy;i++)
		src_xy[i] = (i-(nxy-1)/2.0f)*res;

	__m128 cos_2a;
	__m128 pinhole_norm;
	__m128 current_loc;
	__m128 delta;
	__m128 apex;	// apex of acceptance cone
	__m128 temp1;
	__m128 temp2;
	__m128 temp3;

	float in_cone;	// used as a boolean, but needs to work with SSE
	float dl;		// length of voxel traverse (in voxels)

	for(angle=0;angle<n_proj;angle++)
	{
#ifdef _WIN32
		cout << "Creating map for detector " << proj->detector[angle] << ", angular view " << proj->angular_view[angle] << " (" << proj->GetAngle(angle) << "\xB0) (";
#else
		cout << "Creating map for detector " << proj->detector[angle] << ", angular view " << proj->angular_view[angle] << " (" << proj->GetAngle(angle) << "\xC2\xB0) (";
#endif
		// rotate the original
		RotateOut(-proj->GetAngle(angle),0,proj->GetZShift(angle));

		for(pinhole=0;pinhole<n_pinholes;pinhole++)
		{
			cout << pinhole+1;
			if((pinhole+1) < n_pinholes)
				cout << ", ";
			// load pinhole data
			pinhole_geo[pinhole]->GetNormal(temp_fl);
			pinhole_norm = _mm_load_ps(temp_fl);;

			// locate apex of cone
			temp1 = _mm_set_ss(pinhole_geo[pinhole]->GetDiameter() / (2.0f * tan(pinhole_geo[pinhole]->GetConeAngle())));
			temp1 = _mm_shuffle_ps(temp1,temp1,_MM_SHUFFLE(1,0,0,0));
			temp1 = _mm_mul_ps(pinhole_norm,temp1);
			pinhole_geo[pinhole]->GetLocation(temp_fl);
			apex = _mm_load_ps(temp_fl);
			apex = _mm_add_ps(apex,temp1);

			// angle is opened up a bit so the edge of the attenuation map doesn't limit the projection
			cos_2a = _mm_set_ss(0.5); //45 deg map (cos^2 = 0.5)

			// ray trace each pixel
			for(i=0;i<nxy;i++)
			{
				for(j=0;j<nxy;j++)
				{
					for(k=0;k<nxy;k++)
					{

						amap[pinhole][i][j][k] = 0.0f;

						// check that voxel is in the cone
						current_loc = _mm_set_ps(0.0f,src_z[i],src_xy[j],src_xy[k]);
						delta = _mm_sub_ps(apex,current_loc);
						//temp2 = _mm_dp_ps(delta,delta,0x71);
						temp2 = _mm_setzero_ps();
						temp3 = _mm_mul_ps(delta,delta);
						temp2 = _mm_add_ss(temp2, temp3);
						temp3 = _mm_shuffle_ps(temp3,temp3,_MM_SHUFFLE(0,3,2,1));
						temp2 = _mm_add_ss(temp2, temp3);
						temp3 = _mm_shuffle_ps(temp3,temp3,_MM_SHUFFLE(0,3,2,1));
						temp2 = _mm_add_ss(temp2, temp3);

						//temp1 = _mm_dp_ps(delta,pinhole_norm,0x71);
						temp1 = _mm_setzero_ps();
						temp3 = _mm_mul_ps(delta,pinhole_norm);
						temp1 = _mm_add_ss(temp1, temp3);
						temp3 = _mm_shuffle_ps(temp3,temp3,_MM_SHUFFLE(0,3,2,1));
						temp1 = _mm_add_ss(temp1, temp3);
						temp3 = _mm_shuffle_ps(temp3,temp3,_MM_SHUFFLE(0,3,2,1));
						temp1 = _mm_add_ss(temp1, temp3);

						temp1 = _mm_mul_ss(temp1, temp1);
						temp1 = _mm_div_ss(temp1, temp2);
						temp1 = _mm_cmplt_ss(cos_2a, temp1);
										
						_mm_store_ss(&in_cone,temp1);
						if(!in_cone)
							continue;

						// set starting location and delta
						temp1 = _mm_shuffle_ps(delta,delta,_MM_SHUFFLE(0,0,0,0));
						delta = _mm_div_ps(delta,temp1);
						current_loc = _mm_set_ps(0.0f,float(i),float(j),float(k));	// current location in voxel coordinates

						mu = 0.0;

						for(m=k; m<nxy; m++)
						{
							// increment position
							current_loc = _mm_add_ps(current_loc,delta);
							_mm_store_ps(vox_pos,current_loc);
							i2 = int(floor(vox_pos[2]));
							j2 = int(floor(vox_pos[1]));
							di = vox_pos[2] - i2;
							dj = vox_pos[1] - j2;

							
							// add cumulative attenuation
							mu += (1-di)*(1-dj)*data2[i2][j2][m] +
								  (1-di)*   dj *data2[i2][j2+1][m] +
								  di *(1-dj)*data2[i2+1][j2][m] +
								  di *   dj *data2[i2+1][j2+1][m]; // 
						}
						//temp2 = _mm_dp_ps(delta,delta,0x71);
						temp2 = _mm_setzero_ps();
						temp3 = _mm_mul_ps(delta,delta);
						temp2 = _mm_add_ss(temp2, temp3);
						temp3 = _mm_shuffle_ps(temp3,temp3,_MM_SHUFFLE(0,3,2,1));
						temp2 = _mm_add_ss(temp2, temp3);
						temp3 = _mm_shuffle_ps(temp3,temp3,_MM_SHUFFLE(0,3,2,1));
						temp2 = _mm_add_ss(temp2, temp3);

						temp2 = _mm_sqrt_ss(temp2);
						_mm_store_ss(&dl,temp2);
						amap[pinhole][i][j][k] = exp(-mu * res * dl * atten_scale * 0.1f); // mu is in cm^-1 and res is in mm
					}
				}
			}
		}
		// save attenuation map for projection angle
		sprintf(filename, "%s" __SL "att_map_%d_%d.bin",working_dir,proj->detector[angle],proj->angular_view[angle]);		
		fout.open(filename,fstream::binary);
		for(pinhole=0;pinhole<n_pinholes;pinhole++)
			for(i=0;i<nxy;i++)
				for(j=0;j<nxy;j++)
					fout.write(reinterpret_cast<char*>(amap[pinhole][i][j]),nxy*sizeof(float));
		fout.close();

		cout << ")" << endl;
	}

	for(pinhole=0;pinhole<n_pinholes;pinhole++)
		for(i=0;i<nxy;i++)
			for(j=0;j<nxy;j++)
				delete [] amap[pinhole][i][j];
	for(pinhole=0;pinhole<n_pinholes;pinhole++)
		for(i=0;i<nxy;i++)
			delete [] amap[pinhole][i];
	for(pinhole=0;pinhole<n_pinholes;pinhole++)
	{	
		delete [] amap[pinhole];
		delete pinhole_geo[pinhole];
	}

	delete [] amap;
	delete [] pinhole_geo;

}


void resort(int* list, int n)
{
	int i;

	int* list1;
	int* list2;
	int list1sz, list2sz;

	if(n>2)
	{
		if(n%2)
		{
			list1sz = (n+1)/2;
			list2sz = n - list1sz;
		}
		else
			list1sz = list2sz = n/2;

		list1 = new int[list1sz];
		list2 = new int[list2sz];

		// take the list apart
		for(i=0;i<list1sz;i++)
			list1[i] = list[i];
		for(i=0;i<list2sz;i++)
			list2[i] = list[list1sz+i];

		// resort each half
		resort(list1,list1sz);
		resort(list2,list2sz);

		// and put them back together
		for(i=0;i<list1sz;i++)
			list[2*i] = list1[i];
		for(i=0;i<list2sz;i++)
			list[2*i+1] = list2[i];
	}
}
