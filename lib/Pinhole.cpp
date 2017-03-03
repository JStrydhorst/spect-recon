#include "Pinhole.h"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <fstream>
	
Pinhole::Pinhole(float x, float y, float z, float dia, float angle, float fy, float fz)
{
	//fVector _norm;
	__m128 temp;
	__m128 temp2;

	location = _mm_set_ps(0,z,y,x);

	cone_angle = angle;
	diameter = dia;

	normal = _mm_set_ps(0,z-fz,y-fy,x);
	//temp = _mm_dp_ps(normal,normal,0x7F);
	temp = _mm_setzero_ps();
	temp2 = _mm_mul_ps(normal,normal);
	temp = _mm_add_ss(temp, temp2);
	temp2 = _mm_shuffle_ps(temp2,temp2,_MM_SHUFFLE(0,3,2,1));
	temp = _mm_add_ss(temp, temp2);
	temp2 = _mm_shuffle_ps(temp2,temp2,_MM_SHUFFLE(0,3,2,1));
	temp = _mm_add_ss(temp, temp2);

	temp = _mm_sqrt_ps(temp);
	temp = _mm_shuffle_ps(temp,temp,_MM_SHUFFLE(0,0,0,0));
	normal = _mm_div_ps(normal,temp);
}

void* Pinhole::operator new(size_t size)
{
	void *temp,*ptr;
	temp = malloc(size+15+sizeof(void*));
	ptr = (void*)( ((std::uintptr_t)temp + sizeof(void*) + (size_t)15) & ~0xF);
	*((void**)ptr-1) = temp;
	return ptr;
}

void Pinhole::operator delete(void* ptr)
{
	void* temp;
	temp = *((void**)ptr-1);
	free(temp);
}

Pinhole** Pinhole::LoadPinholeGeometry(char* filename, float COR, int &n_pinholes) // loads the pinhole geometry from a file
{
	int i;
	char temp_str[64];

	Pinhole** pin_array;

	float y, z;
	float dia, cone_angle;
	float py, pz;

	// open pinhole description file
	std::ifstream fin;
	fin.open(filename);

	// parse file
	fin >> temp_str;
	if (temp_str[0] == '#')
		fin.ignore(256,'\n'); // ignore first line if it starts with #

	fin >> temp_str;
	if(temp_str[0] == '[')
		std::cout << "\nLoading pinhole definition " << temp_str << std::endl;
	else
	{
		std::cout << "\nInvalid pinhole definition file." << std::endl;
		fin.close();
		return NULL;
	}

	fin >> n_pinholes;
	pin_array = new Pinhole*[n_pinholes];

	for(i=0;i<n_pinholes;i++)
	{
		fin >> y >> z >> dia >> cone_angle >> py >> pz;
		cone_angle *= 4.0f * atan(1.0f) / 180.0f;
		pin_array[i] = new Pinhole(COR, y, z, dia, cone_angle, py, pz);
	}

	fin.close();

	return pin_array;
}
