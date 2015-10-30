// fourier transforms that handle the data in a slightly more convient form - arrays of reals and imaginaries
void fft1(float *data_r, float *data_i, int nx, int isign);
void fft2(float **data_r, float **data_i, int nx, int ny, int isign);	// 2D fourier transform
void fft3(float ***data_r, float ***data_i, int nx, int ny, int nz, int isign);	// 3D fourier transform

