#ifndef __RECONSTRUCTION_H
#define __RECONSTRUCTION_H

#include "Projection.h"
#include "ProjTable.h"

/********************************************************************************************
/ Reconstruction class
*********************************************************************************************/
class Reconstruction
{
public:
	// constructors/destructor
	Reconstruction(int size_xy, int size_z, int timeslots, float new_res, char **files=0);
	Reconstruction(const char* filename);
	~Reconstruction();

	void Clear(int timeslot);	// erases the contents of data
	void Project(Projection* proj, const char* mmp_atten = 0);
	void Backproject(const Projection* proj, const char* mmp_atten = 0);
	void Update(const Reconstruction* scale, const Reconstruction* norm, const Reconstruction* MAPPartial, int timeslot, float threshold);

	void UpdateMAPSpatialPartial(const Reconstruction* param, int timeslot, float beta); //, float delta);		// updates the partial derivative of the energy function
//	void UpdateCTMAPPartial(const Reconstruction* param, const Reconstruction* CT, int timeslot, float beta, float delta);
	void UpdateMAPTemporalPartial(const Reconstruction* param, int timeslot, float beta);

	// make the images pretty (updates all timeslots)
	void Filter(float sigma);		// 3D Gaussian blur, with sigma (in mm) specified
	void CleanUpStrayPixels(float limit); // sets pixels above limit to limit, clears pixels not sampled at all
	void Mask(int radius);					// clears pixels outside a certain radius

	// projection table handling
	int LoadProjTables(const Projection* proj, const char* dir);
	void UnloadProjTables();
	bool ProjTableLoaded();

	// extract information about the reconstruction
	void ReconstructionInfo();
	float Activity();						// returns total activity
	float MSE(Reconstruction* param);		// returns mean squared error relative to supplied recon
	float VOIActivity(int timeslot = 0);	// returns the activity within a fixed VOI hardcoded in the function

	// export to a file
	void SaveBIN(const char* filename, int mask=64);		// can mask the edges of the recon when saving. Use with caution.
	void LoadBIN(const char* filename, bool duplicate_gate0);
	void SaveDICOM(const char* filename, const char* ProjDCM, int r);

	// test functions and digital phantoms
	void SetSingleVoxel(int k, int j, int i) { data[0][i][j][k] = 1.0;}
	void SetPlane(int k);
	void CylindricalPhantom(float radius, int first, int last, float offset);
	void Crosshairs();
	void HotRodPhantom();
	void ColdRodPhantom();
	void TestRotation();

	void CreateMMPAttenMap(const char* working_dir, char* pinhole_file, Projection* proj, float atten_scale=1.0f);

private:
	void InitRotationMatrix();
	void RotateOut(float angle, int timeslot, float z_offset, float xshift=0.0f, float yshift=0.0f); // puts data from data into data2
	void RotateIn(float angle,  int timeslot, float z_offset,float xshift=0.0f, float yshift=0.0); // adds data from data2 into data
	void Clear2();	// erases the contents of data2
	void CopyGate(int dest_timeslot, int src_timeslot);	// copies the content of one time gate to another (is this ever used???)

	inline float phi(float u, float delta);	// penalty function for MAP reconstruction

	int n_timeslots;
	int nxy, nz;
	float res;
	float *_data;		// reconstruction data stored in a consecutive block
	float ****data;		// stores original
	float *_data2;
	float ***data2;		// data structure used for rotations and accumulations and other fun stuff
	float **rot_mask;	// 

	static ProjTable** proj_table;
	static int n_proj_tables;

	static float **sp_matrix;	// matrix used for gaussian rotation
	static bool rotation_matrix_initialized;

};

#endif // __RECONSTRUCTION.H
