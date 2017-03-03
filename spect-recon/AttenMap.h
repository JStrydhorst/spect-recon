#ifndef __ATTENMAP_H
#define __ATTENMAP_H

/********************************************************************************************
/ Attenuation map
********************************************************************************************/
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

#endif // __ATTENMAP_H
