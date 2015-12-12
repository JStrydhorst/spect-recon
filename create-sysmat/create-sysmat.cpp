#include <iostream>
#include <fstream>
#include <smmintrin.h>

using namespace std;

#include "../lib/pinhole.h"
#include "../lib/proj_table.h"

#define N_SUBPIX 10
#define N_SUBVOX 2

void intersect_plane_line(__m128 &p0, __m128 &norm, __m128 &l0, __m128 &l, __m128 &intersect) ;
bool in_cone(__m128 &point, __m128 &apex, __m128 &normal, float &cos2alpha);

void SaveBin(char* filename, float** image, int nx, int ny);

void main()
{
	char filename[260];

	int ph;
	int i,j,k;	// voxel indices (i=z, j=y, k=x)
//	int vox_index;
	int m,n;	// pixel indices (m=z, n=y)
//	int pix_index;
	int sub_i, sub_j, sub_k;
	int sub_m, sub_n;

	int m_min, m_max;
	int n_min, n_max;

	ProjTable* proj_table;

	// reconstruction geometry
	int nz = 128, nxy = 128;
	float res = 0.6f;
/*	int nz = 32, nxy = 32;
	float res = 2.4f; */
/*	int nz = 76, nxy = 76;
	float res = 1.0f;*/

	// projection geometry
	int proj_nyz = 256;	// 256 x 256 detector grid
	float proj_res = 1.0f;
	
	// pinhole geometry
	float focal_length = 131.0f;	// focal length of nanospect collimator
	float RoR = 45.0f;
	char aperture[] = "APT2";
	int n_pinholes;
	sprintf(filename,"%s.pin",aperture);
	Pinhole** pinholes = LoadPinholeGeometry(filename, RoR, n_pinholes);
	__m128 plane1 = _mm_set_ps(0.0,0.0,0.0,40.0f);	// back of collimator is nominally 45 mm from axis
	__m128 plane2 = _mm_set_ps(0.0,0.0,0.0,50.0f);

	//float mu = 3.31f;	// attenuation of 140keV photons /mm


	// working variables
	float dt;
	float omega;	// transmission, solid angle

	float pinhole_rad;
	__m128 pinhole_loc, normal;
	__m128 src, dest;

	__m128 apex1, apex2;
	float cos_2a;
	__m128 temp1, temp2;

	__m128 q, l, intersect;
	float d;

	float sum;

	__declspec(align(16)) float fl_array[4];

	// precalculate x, y, z
	float *src_x, *src_y, *src_z;
	src_x = new float[nxy];
	for(i=0;i<nxy;i++)
		src_x[i] = (i-(nxy-1)/2.0f)*res;
	src_y = new float[nxy];
	for(i=0;i<nxy;i++)
		src_y[i] = (i-(nxy-1)/2.0f)*res;
	src_z = new float[nz];
	for(i=0;i<nz;i++)
		src_z[i] = (i-(nz-1)/2.0f)*res;

	//float subvox[] = {-res/4, res/4};
	float subvox[N_SUBVOX];
	for(i=0;i<N_SUBVOX;i++)
		subvox[i] = (i-(N_SUBVOX-1)/2.0f)*(res/N_SUBVOX);

	// precalculate detector position lookup table
	float *dest_y, *dest_z;
	dest_z = new float[proj_nyz];
	for(i=0;i<proj_nyz;i++)
		dest_z[i] = (i-(proj_nyz-1)/2.0f) * proj_res;
	dest_y = new float[proj_nyz];
	for(i=0;i<proj_nyz;i++)
		dest_y[i] = (i-(proj_nyz-1)/2.0f) * proj_res;
	
	float subpix[N_SUBPIX];
	for(i=0;i<N_SUBPIX;i++)
		subpix[i] = (i-(N_SUBPIX-1)/2.0f)*(proj_res/N_SUBPIX);

	// create rotation mask
	float** rot_mask;
	rot_mask = new float*[nxy];
	for(i=0;i<nxy;i++)
		rot_mask[i] = new float[nxy];
	
	unsigned int rad_sq = 57; // 2 pixel margin allows 5 pixel gaussian matrix (usually use 3 anyways)
	rad_sq *= rad_sq;

	float xy_offset = (nxy-1)/2.0f;
	for(i=0;i<nxy;i++)
		for(j=0;j<nxy;j++)
			if( ((i-xy_offset)*(i-xy_offset)+(j-xy_offset)*(j-xy_offset)) < rad_sq)
				rot_mask[i][j] = 1;
			else
				rot_mask[i][j] = 0;

	for(ph=0;ph<n_pinholes;ph++)
	{
		cout << "Pinhole: " << ph+1 << endl;
		cout << "(" << pinholes[ph]->location.m128_f32[0] << ", " << pinholes[ph]->location.m128_f32[1] << ", "<< pinholes[ph]->location.m128_f32[2] << ")" << endl; 

		// create a table for each pinhole
		proj_table = new ProjTable();

		pinhole_loc = pinholes[ph]->location;
		normal = pinholes[ph]->normal;

		// calculate the square of the cone opening angle for later use
		cos_2a = cos(pinholes[ph]->cone_angle) * cos(pinholes[ph]->cone_angle);

		// calculate effective pinhole radius
		pinhole_rad = (pinholes[ph]->diameter)/2.0f; // for tungsten and 140keV, penetration is minimal

		// locate apexes of cones
		dt = pinholes[ph]->diameter / (2.0f * tanf(pinholes[ph]->cone_angle));
		temp1 = _mm_set_ss(dt);
		temp1 = _mm_shuffle_ps(temp1,temp1,_MM_SHUFFLE(1,0,0,0));
		temp1 = _mm_mul_ps(normal,temp1);
		apex1 = _mm_add_ps(pinhole_loc,temp1);
		apex2 = _mm_sub_ps(pinhole_loc,temp1);

		for(j=0;j<nxy;j++)
		{
			cout << ".";
			for(k=0;k<nxy;k++)
				if(rot_mask[j][k])
				{
					for(i=0;i<nxy;i++)
					{
						// calculate scan range for m and n
						src = _mm_set_ps(0.0f, src_z[i], src_y[j],src_x[k]);

						temp1 = _mm_set_ps(0.0f, 0.0f, 0.0f, (focal_length+RoR));
						temp1 = _mm_sub_ss(temp1, src);	// src to det dist
						temp2 = _mm_sub_ps(pinhole_loc,src); // src to pinhole vector
						temp1 = _mm_div_ss(temp1, temp2);
						temp1 = _mm_shuffle_ps(temp1,temp1,_MM_SHUFFLE(1,0,0,0));
						temp1 = _mm_mul_ps(temp1,temp2);
						temp1 = _mm_add_ps(temp1, src);
						_mm_store_ps(fl_array,temp1);

						m_min = int(fl_array[2]/proj_res) + (proj_nyz/2) - 20; // could make the size of the scan depend on the magnification...
						m_max = m_min + 40;
						n_min = int(fl_array[1]/proj_res) + (proj_nyz/2) - 20;  
						n_max = n_min + 40;

						if(m_min<21)
							m_min = 21;
						if(m_max>234)
							m_max = 234;
						if(n_min<21)
							n_min = 21;
						if(n_max>234)
							n_max = 234;

						for(m=m_min;m<m_max;m++)
							for(n=n_min;n<n_max;n++)
							{
								sum = 0.0f;

								for(sub_i=0;sub_i<N_SUBVOX;sub_i++)
									for(sub_j=0;sub_j<N_SUBVOX;sub_j++)
										for(sub_k=0;sub_k<N_SUBVOX;sub_k++)
											for(sub_m=0;sub_m<N_SUBPIX;sub_m++)
												for(sub_n=0;sub_n<N_SUBPIX;sub_n++)
												{
													src = _mm_set_ps(0.0f, src_z[i]+subvox[sub_i], src_y[j]+subvox[sub_j], src_x[k]+subvox[sub_k]);
													//src = _mm_set_ps(0.0f, src_z[i], src_y[j],src_x[k]);
													dest = _mm_set_ps(0, dest_z[m]+subpix[sub_m], dest_y[n]+subpix[sub_n], (focal_length+RoR));

													q = _mm_sub_ps(dest, src);

													// is plane 1 intersection in cone 1?
													temp1 = _mm_sub_ss(plane1,src);
													temp1 = _mm_div_ss(temp1,q);
													temp1 = _mm_shuffle_ps(temp1,temp1,_MM_SHUFFLE(1,0,0,0));
													temp1 = _mm_mul_ps(temp1,q);
													temp2 = _mm_add_ps(src,temp1);
													if(!in_cone(temp2,apex1,normal,cos_2a))
														continue;	// ray hits surface of collimator

													// is plane 2 intersection in cone 2?
													temp1 = _mm_sub_ss(plane2,src);
													temp1 = _mm_div_ss(temp1,q);
													temp1 = _mm_shuffle_ps(temp1,temp1,_MM_SHUFFLE(1,0,0,0));
													temp1 = _mm_mul_ps(temp1,q);
													temp2 = _mm_add_ps(src,temp1);
													if(!in_cone(temp2,apex2,normal,cos_2a))
														continue;	// ray hits surface of collimator
														
													// does ray pass through pinhole?
													temp1 = _mm_dp_ps(q,q,0x71);
													temp1 = _mm_sqrt_ss(temp1);
													temp1 = _mm_shuffle_ps(temp1,temp1,_MM_SHUFFLE(1,0,0,0));
													l = _mm_div_ps(q,temp1);
													intersect_plane_line(pinhole_loc,normal,src,l,intersect);
													temp1 = _mm_sub_ps(intersect,pinhole_loc);
													temp1 = _mm_dp_ps(temp1,temp1,0x71);
													temp1 = _mm_sqrt_ss(temp1);
													_mm_store_ss(&d,temp1);
													if(d > pinhole_rad)
														continue;

													// calculate solid angle
													/*
													temp1 = _mm_set_ps(0.0f, 0.0f, 0.0f, (focal_length+RoR));		// hardcoded detector location
													temp2 = _mm_dp_ps(q,q,0x71);
													temp2 = _mm_sqrt_ss(temp2);
													temp1 = _mm_div_ss(temp1, temp2);					// temp 1 = cos(theta)
													temp1 = _mm_mul_ss(temp1, _mm_mul_ss(temp1, temp1)); // cos^3
													_mm_store_ss(&omega, temp1);
													*/
													temp1 = _mm_set_ps(0.0f, 0.0f, 0.0f, (focal_length+RoR));		// hardcoded detector location
													temp2 = _mm_dp_ps(q,q,0x71);
													temp2 = _mm_sqrt_ss(temp2);
													temp1 = _mm_div_ss(temp1, temp2);					// purely 1/r^2 with no cos theta dependance
													temp2 = _mm_div_ss(q,temp2);						// q_x / |q|
													temp1 = _mm_mul_ss(temp1, temp1);					//
													temp1 = _mm_mul_ss(temp1, temp2);
													_mm_store_ss(&omega, temp1);

													sum += omega;
												}

								if(sum>1.0f)
									proj_table->AddEntry(i,j,k,m,n,sum/(N_SUBVOX*N_SUBVOX*N_SUBVOX*N_SUBPIX*N_SUBPIX));
												
							}
							

					}
				}
		}

		cout << endl;
		cout << "Table length: " << proj_table->n << " (" << proj_table->n*9/1024 << " kB)" << endl;

		sprintf_s(filename, "%s_%d_%d_%3.1f_%d_%d_%3.1f_%d.pta",aperture,nxy,nz,res,proj_nyz,proj_nyz,proj_res,ph);

		proj_table->SaveProjTable(filename);

		delete proj_table;
		
	} // ph

	for(i=0;i<n_pinholes;i++)
		delete pinholes[i];
	delete [] pinholes;
}

void intersect_plane_line(__m128 &p0, __m128 &norm, __m128 &l0, __m128 &l, __m128 &intersect) 
{
	__m128 num, denom;
	__m128 d;

	num = _mm_sub_ps(p0,l0);
	num = _mm_dp_ps(num,norm,0x71);
	denom = _mm_dp_ps(l,norm,0x71);
	d = _mm_div_ss(num,denom);
	d = _mm_shuffle_ps(d,d,_MM_SHUFFLE(1,0,0,0));
	intersect = _mm_add_ps(l0,_mm_mul_ps(d,l));
}

bool in_cone(__m128 &point, __m128 &apex, __m128 &normal, float &cos2alpha)
{
	float cos2theta;	// cos of angle between r and normal, squared
	__m128 r = _mm_sub_ps(point,apex);	// apex to point vector
	__m128 temp1, temp2;
	temp1 = _mm_dp_ps(r,normal,0x71);
	temp1 = _mm_mul_ss(temp1, temp1);
	temp2 = _mm_dp_ps(r,r,0x71);
	temp1 = _mm_div_ss(temp1,temp2);
	_mm_store_ss(&cos2theta,temp1);

	if(cos2theta > cos2alpha)
		return true;
	else
		return false;
}

void SaveBin(char* filename, float** image, int nx, int ny)
{
	ofstream f;

	f.open(filename,ios::out|ios::binary);
	for(int i=0;i<ny;i++)
		f.write(reinterpret_cast<char*>(image[i]),nx*sizeof(float));
	f.close();
}