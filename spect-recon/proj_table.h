typedef unsigned char uchar;

/********************************************************************************************
/ ProjTable class
********************************************************************************************/
class ProjTable
{
public:
	ProjTable(char* filename = NULL);
	~ProjTable();

	ProjTable* Abridge(float threshold);

	void AddEntry(uchar new_ri, uchar new_rj, uchar new_rk, uchar new_pj, uchar new_pk, float new_weight);
	void SaveProjTable(char* filename);
	int LoadProjTable(char* filename);
	void ProjTableStats();	// dumps projection table stats to screen

	friend class Reconstruction;

	unsigned int n;					// number of entries in table

private:
	unsigned int mem_alloc;			// number of entries allocated
	unsigned char *ri, *rj, *rk;	// reconstruction voxel
	unsigned char *pi, *pj;			// projection pixel
	float *weight;					// weighting of projected voxel
};

