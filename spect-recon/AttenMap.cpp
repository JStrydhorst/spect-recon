#include "AttenMap.h"

#include <cstdlib>
#include <iostream>
#include <fstream>

AttenMap::AttenMap()
{
	num_proj = 0;
	nz = 0;
	nxy = 0;
	atten_map = 0;
}

AttenMap::AttenMap(const char *filename, int new_num_proj, int new_nz, int new_nxy)
{	
	int angle, i, j;
	std::ifstream fin;

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
	catch(std::bad_alloc& ba)
	{
		std::cout << "bad_alloc caught in AttenMap::AttenMAP  --  " << ba.what() << std::endl;
		exit(-1);
	}

	fin.open(filename,std::fstream::binary);
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
