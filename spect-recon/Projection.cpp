#include "Projection.h"
#include "RootDicomObj.h"

#include "../lib/fft.h"

#include<cstdlib>
#include<cstring>
#include<cmath>
#include<iostream>
#include<fstream>

/********************************************************************************************
/ Projection
********************************************************************************************/
Projection::Projection()	// default constructon, empty projection, is this ever called???
{
	int i;

	n_timeslots = 1;
	ny = 0;
	nz = 0;
	res = 0;
	num_proj = 0;
	angles = 0;
	z_shift = 0;
	detector = 0;
	angular_view = 0;
	timeslot = 0;
	frame_time = 0;
	_data = 0;
	data = 0;

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
	std::ifstream fin;

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
	catch(std::bad_alloc& ba)
	{
		std::cout << "std::bad_alloc caught in Projection::Projection(int, int, int, etc... )  --  " << ba.what() << std::endl;
		exit(-1);
	}

	if(filename)
	{
		fin.open(filename,std::ios::binary);
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

	if(DCM->FindElement(0x0010,0x0010)==0) // look for mandatory patient name field
	{
		std::cout << "No valid DICOM data loaded." << std::endl;
		exit(-1);
	}

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
		std::cout << "Fixing aperture name!" << std::endl;
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
	catch(std::bad_alloc& ba)
	{
		std::cout << "std::bad_alloc caught in Projection::Projection(const char*)  --  " << ba.what() << std::endl;
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
	std::cout << "******Projection Info************" << std::endl;
	std::cout << "Number of projections: " << num_proj << std::endl;
	std::cout << nz << " x " << ny << " (" << res << " mm)" << std::endl;
	std::cout << "Aperture: " << aperture << std::endl;
	std::cout << "Scan Length: " << GetScanLength() << " mm" << std::endl;
	std::cout << "Scan Duration: " << GetScanDuration() << " min" << std::endl;
	if(n_timeslots>1)
		std::cout << "Gates: " << GetNumTimeSlots() << std::endl;
	std::cout << "*********************************" << std::endl << std::endl;
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
	catch(std::bad_alloc& ba)
	{
		std::cout << "std::bad_alloc caught in Projection::Projection(int, int, int, etc... )  --  " << ba.what() << std::endl;
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

	std::ofstream fout;
	
	fout.open(filename,std::fstream::binary);
	for(i=0;i<num_proj;i++)
		for(j=0;j<nz;j++)
			fout.write(reinterpret_cast<char*>(data[i][j]),ny*sizeof(float));
	fout.close();

}

void Projection::LoadBIN(const char* filename)
{
	int i,j;

	std::ifstream fout;
	fout.open(filename,std::fstream::binary);
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
		std::cout << "Projections must be sorted by time slot before creating subsets" << std::endl;
		return 0;
	}

	if(num_proj%n_sets)
	{
		std::cout << "Total projections must divide evenly..." << std::endl;
		return 0;
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

	std::ofstream fout;
	std::ifstream header_info;
	
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

	DE = new DataElement(0x0054,0x0022,"SQ",0,0);
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

	DE = new DataElement(0x0054,0x0052,"SQ",0,0);
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

	fout.open(filename,std::ofstream::binary);
	proj->Write(fout);
	fout.close();

	delete proj;

}

void Projection::resort(int* list, int n)
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
                       
