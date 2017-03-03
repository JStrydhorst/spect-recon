// Fast Fourier Transforms 

#include "fft.h"

extern "C" void cdft(int, int, double*, int*, double*);

#include <cmath>
#include <algorithm>

void fft1(float *data_r, float *data_i, int nx, int isign)
{
	int i;
	double *data = new double[2 * nx];
	int *ip = new int[2+int(sqrt(nx))+1];
	ip[0] = 0;
	double *w = new double[nx/2];
	
	for(i=0;i<nx;i++)
	{
		data[2*i] = data_r[i];
		data[2*i+1] = data_i[i];
	}
	cdft(2*nx, isign, data, ip, w);
	for(i=0;i<nx;i++)
	{
		data_r[i] = float(data[2*i]);
		data_i[i] = float(data[2*i+1]);
	}

	delete [] data;

	delete[] ip;
	delete[] w;
}

void fft2(float **data_r, float **data_i, int nx, int ny, int isign)
{
	int i,j;
	double *x_data;
	double **y_data;

	y_data = new double*[nx];
	for(i=0;i<nx;i++)
		y_data[i] = new double[2*ny];

	x_data = new double[2*nx];

	int *ip = new int[2 + int(sqrt(std::max(nx,ny))) + 1];
	ip[0] = 0;
	double *w = new double[std::max(nx,ny) / 2];

	// transform along x-axis
	for(i=0;i<ny;i++)
	{
		// copy row to temp data
		for(j=0;j<nx;j++)
		{
			x_data[2*j] = data_r[i][j];
			x_data[2*j+1] = data_i[i][j];
		}

		cdft(2 * nx, isign, x_data, ip,w);
		
		// copy to temporary structure
		for(j=0;j<nx;j++)
		{
			y_data[j][2*i] = x_data[2*j];
			y_data[j][2*i+1] = x_data[2*j+1];
		}
	}

	ip[0] = 0;

	// transform along y_axis and put back in original format
	for(j=0;j<nx;j++)
	{
		cdft(2 * ny, isign, y_data[j], ip, w);

		for(i=0;i<ny;i++)
		{
			data_r[i][j] = float(y_data[j][2*i]);
			data_i[i][j] = float(y_data[j][2*i+1]);
		}
	}

	delete [] x_data;

	for(i=0;i<nx;i++)
		delete [] y_data[i];
	delete [] y_data;

	delete[] ip;
	delete[] w;

}

void fft3(float ***data_r, float ***data_i, int nx, int ny, int nz, int isign)
{
	int i, j, k;

	int nn;
	nn = std::max(nx, std::max(ny, nz));

	double *temp = new double[nn * 2];

	int *ip = new int[2 + int(sqrt(nn)) + 1];
	ip[0] = 0;
	double *w = new double[nn / 2];

	// transform in x direction
	for (i = 0;i < nz;i++)
		for (j = 0;j < ny;j++)
		{
			for (k = 0;k < nx;k++)
			{
				temp[2 * k] = data_r[i][j][k];
				temp[2 * k + 1] = data_i[i][j][k];
			}
			cdft(2 * nx, isign, temp, ip, w);
			for (k = 0;k < nx;k++)
			{
				data_r[i][j][k] = float(temp[2 * k]);
				data_i[i][j][k] = float(temp[2 * k + 1]);
			}
		}

	ip[0] = 0;
	// transform in y direction
	for (i = 0;i < nz;i++)
		for (k = 0;k < nx;k++)
		{
			for (j = 0;j < ny;j++)
			{
				temp[2 * j] = data_r[i][j][k];
				temp[2 * j + 1] = data_i[i][j][k];
			}
			cdft(2 * ny, isign, temp, ip, w);
			for (j = 0;j < ny;j++)
			{
				data_r[i][j][k] = float(temp[2 * j]);
				data_i[i][j][k] = float(temp[2 * j + 1]);
			}
		}

	ip[0] = 0;
	// transform in z direction
	for(j=0;j<ny;j++)
		for(k=0;k<nx;k++)
		{
			for(i=0;i<nz;i++)
			{
				temp[2*i] = data_r[i][j][k];
				temp[2*i+1] = data_i[i][j][k];
			}
			cdft(2 * nz, isign, temp, ip, w);
			for(i=0;i<nz;i++)
			{
				data_r[i][j][k] = float(temp[2*i]);
				data_i[i][j][k] = float(temp[2*i+1]);
			}
		}

	delete[] temp;
  	delete[] ip;
	delete[] w;

}
