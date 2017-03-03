#include "Reconstruction.h"
#include "AttenMap.h"

#include "RootDicomObj.h"
#include "Pinhole.h"
#include "fft.h"

#include "version.h"

#include <cstring>
#include <cmath>

/* Gaussian Rotator Parameters */
#define ROT_SUBPIX      65
#define ROT_MAT_SIZE    3 // must be an odd number (3 or 5)
#define ROT_FWHM        1.0f

#define __SL "/"


/********************************************************************************************
/ Reconstruction
********************************************************************************************/

float** Reconstruction::sp_matrix = 0;
bool Reconstruction::rotation_matrix_initialized = false;
ProjTable** Reconstruction::proj_table = 0;
int Reconstruction::n_proj_tables = 0;

// constructor
Reconstruction::Reconstruction(int size_xy, int size_z, int timeslots, float new_res, char** filename)
{
	int t,i,j;
	std::ifstream f_in;
	std::ifstream::pos_type size;
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
	catch(std::bad_alloc& ba)
	{
		std::cout << "Reconstruction(int,int,float,char*):std::bad_alloc caught: " << ba.what() << std::endl;
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
			f_in.open(filename[t],std::ios::in|std::ios::binary|std::ios::ate);
			if(f_in.is_open())
			{
				size = f_in.tellg();
				f_in.seekg(0,std::ios::beg);

				if( (size % (nxy*nxy*sizeof(float))) != 0)
				{
					std::cout << "Data file doesn't contain an integral number of slices." << std::endl;
				}
				else
				{
					slices = int(size)/(nxy*nxy*sizeof(float));

					if(slices>nz)
						std::cout << "Data file exceeds size of reconstruction." << std::endl;
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
	std::cout << "**********Reconstruction Info***********" << std::endl;
	std::cout << nz << "x" << nxy << "x" << nxy << " (" << res << " mm)" << std::endl;
	std::cout << "Time slots: " << n_timeslots << std::endl;
	std::cout << "****************************************" << std::endl << std::endl;

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

	int i_offset =  int(std::floor(z_offset/res));
	float i_fr = (z_offset/res) - i_offset;

	//std::cout << "(" << z_offset << ")";

	if (i_offset < 0)
	{
		//std::cout << std::endl << "warning: clipping in z-dim" << " (" << z_offset << ")" << std::endl;
		i_min = -i_offset;
	}
	else
		i_min = 0;

	if(i_offset + nxy >= nz)
	{
		//std::cout << std::endl << "warning: clipping in z-dim" << " (" << z_offset << ")" << std::endl;
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

	int i_offset = int(std::floor(z_offset/res));
	float i_fr = (z_offset/res) - i_offset;

	//std::cout << "(" << z_offset << ")";

	if (i_offset < 0)
	{
		//std::cout << std::endl << "warning: clipping in z-dim" << " (" << z_offset << ")" << std::endl;
		i_min = -i_offset;
	}
	else
		i_min = 0;

	if(i_offset + nxy >= nz)
	{
		//std::cout << std::endl << "warning: clipping in z-dim" << " (" << z_offset << ")" << std::endl;
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
								//std::cout << t;
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
	std::ofstream f;
	f.open("c:\\SPECT\\filter.txt");
	f << "XY Filter (" << len_xy << ")" << std::endl;
	for(i=0;i<len_xy;i++)
		f << filt_xy[i] << ", ";
	f << std::endl << std::endl;
	f << "Z Filter (" << len_z << ")" << std::endl;
	for(i=0;i<len_xy;i++)
		f << filt_xy[i] << ", ";
	f << std::endl;
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
	std::ofstream fout;

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
		

		fout.open(_filename,std::ofstream::binary);
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
	std::ifstream fin;

	//TODO: check for file existence and size compatibility

	fin.open(filename,std::ofstream::binary);
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
	std::ifstream fin;

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
		//std::cout << "Loading " << filename << std::endl;
		try
		{
			proj_table[n_proj_tables] = new ProjTable(filename);
		}
		catch (std::bad_alloc &ba)
		{
			std::cout << "std::bad_alloc caught: " << ba.what() << std::endl;
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

	std::cout << "Projecting ";

	// clear the projections
	proj->Clear();

	for(angle = 0;angle<proj->num_proj;angle++)	// calling this variable angle is confusing, it's the index of the projection, not the actual projection angle
	{
		//std::cout << "(" << proj->angles[angle] << ")";
		std::cout << ".";

		x_shift = (45.0f - float(proj->RadiusOfRotation[proj->detector[angle]-1]));
		y_shift = float(proj->TransaxialApertureOffset[proj->detector[angle]-1]);
		z_shift = float(proj->AxialApertureOffset[proj->detector[angle]-1]);

		//std::cout << proj->detector[angle] << " : (" << x_shift << ", " << y_shift << ", " << z_shift << ")" << std::endl;

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

	std::cout << std::endl;
	
}

/******************************************
/ Backprojects the projection data into the reconstruction
/ Only backprojects into timeslot 0
******************************************/
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

	std::cout << "Backprojecting ";
	Clear(0);
	for(angle=0; angle<proj->num_proj; angle++)
	{
		//std::cout << "(" << proj->angles[angle] << ")";
		std::cout << ".";

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
			catch(std::bad_alloc &ba)
			{
				std::cout << "std::bad_alloc caught: " << ba.what() << std::endl;
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

	std::cout << std::endl;

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

	std::ofstream fout;
	
	char temp_str[64];
	unsigned short temp_us;
//	unsigned short temp_at[2];
	unsigned short* temp_data;
	unsigned long len;

	time_t _Time;
//	tm* timeinfo;
	time(&_Time);
	srand((unsigned int)_Time);
//	timeinfo = localtime(&_Time);

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
	DE = new DataElement(0x0054,0x0012,"SQ",0,0);
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
	DE = new DataElement(0x0054,0x0022,"SQ",0,0);
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

	fout.open(filename,std::ofstream::binary);
	recon->Write(fout);
	fout.close();

	delete recon;
	delete proj;

}

void Reconstruction::CreateMMPAttenMap(const char* working_dir, char* pinhole_file, Projection* proj, float atten_scale)
{

	int angle, pinhole, i, j, k, m;	// current voxel

	float**** amap;
	float mu;

	char filename[64];
	std::ifstream fin;
	std::ofstream fout;
	std::ofstream f_info;

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

	Pinhole** pinhole_geo = Pinhole::LoadPinholeGeometry(pinhole_file, 45.0f, n_pinholes);

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
		std::cout << "Creating map for detector " << proj->detector[angle] << ", angular view " << proj->angular_view[angle] << " (" << proj->GetAngle(angle) << "\xB0) (";
#else
		std::cout << "Creating map for detector " << proj->detector[angle] << ", angular view " << proj->angular_view[angle] << " (" << proj->GetAngle(angle) << "\xC2\xB0) (";
#endif
		// rotate the original
		RotateOut(-proj->GetAngle(angle),0,proj->GetZShift(angle));

		for(pinhole=0;pinhole<n_pinholes;pinhole++)
		{
			std::cout << pinhole+1;
			if((pinhole+1) < n_pinholes)
				std::cout << ", ";
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
							i2 = int(std::floor(vox_pos[2]));
							j2 = int(std::floor(vox_pos[1]));
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
		fout.open(filename,std::ios::binary);
		for(pinhole=0;pinhole<n_pinholes;pinhole++)
			for(i=0;i<nxy;i++)
				for(j=0;j<nxy;j++)
					fout.write(reinterpret_cast<char*>(amap[pinhole][i][j]),nxy*sizeof(float));
		fout.close();

		std::cout << ")" << std::endl;
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

