#ifndef __PINHOLE_H
#define __PINHOLE_H
	
#include <emmintrin.h>

/********************************************************************************************
/ Pinhole class
********************************************************************************************/
class Pinhole
{
public:
	Pinhole(float x, float y, float z, float dia, float angle, float fy, float fz);

	static Pinhole** LoadPinholeGeometry(char* filename, float COR, int &n_pinholes);
	
	void* operator new(size_t);	// needed to align memory allocation for __m128 variables
	void operator delete(void*);

	float GetConeAngle() { return cone_angle; }
	float GetDiameter() { return diameter; }
	void GetLocation(float* loc) { _mm_store_ps(loc,location); }
	void GetNormal(float* norm) { _mm_store_ps(norm,normal); }

	__m128 location;
	__m128 normal;
	float cone_angle;	// opening angle (cone half-angle)
	float diameter;
};

#endif // __PINHOLE_H
