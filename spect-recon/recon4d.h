#define SPECT4D_MAJOR_VERSION	1
#define SPECT4D_MINOR_VERSION	3

class ProjTable;

/********************************************************************************************
/ Projection class
*********************************************************************************************/
class Projection
{
public:
	Projection();
	Projection(int projections, int size_y, int size_z, float res, char* aperture_name, const char* filename = NULL); // blank or binary file
	Projection(const char* filename);		// dicom file
	Projection(const Projection &param);	// copy parameters from another object
	~Projection();

	void ProjectionInfo();	// outputs the projection info to stdout
	void AddTimeslots(int n);
	void SetTimeslot(int timeslot);		// changes the timeslot setting so the correct timeslot is projected

	// manipulation
	void Clear();
	void Ratio(const Projection* param);	// divides param by the current projection
	void Sum(const Projection* param, float k);	// sum projection + k * param
	void Filter(float fwhm);				// applies a Gaussian filter with the specified FWHM
	void ApplyOffsets(int dir);			// shifts Projections based on the detector position offsets (0->offsets ideal projections, 1->idealizes offset projections)
	Projection* CreateNorm(float halflife); // halflife specified in minutes, halflife of 0 -> no decay
	//void ApplyCos3Scaling();			// applies cos3 scaling to projections to simulate sensitivity map of scanner
	void ApplyDecay(float halflife);	// scales the projectionn to account for decay


	// Simulation
	void Scale(float f);		// scales all of the projections by a factor of f
	void Poisson(long seed);				// converts each pixel into a poisson variable with the given lambda

	// File I/O
	void SaveBIN(const char* filename);		// Saves the projections to a .bin file
	void LoadBIN(const char* filename);		// Replaces the currect projection data with data from the specified file
	void SaveDICOM(const char* filename, bool scale);	// scale determines if the values are scaled or integers.

	int GetNumProj() { return num_proj; };
	int GetNumTimeSlots() { return n_timeslots; };
	float GetAngle(int index) { return angles[index]; };
	float GetZShift(int index) { return z_shift[index]; };
	float GetScanLength();
	float GetScanDuration();

	// Subsets...
	Projection** SplitTimeSlots();	// returns a projection set with a single time slot
	Projection** CreateOrderedSubsets(int n_sets);	// returns an array of projections for the specified timeslot...

	friend class Reconstruction;

private:
	int n_timeslots;
	int ny, nz;			// size of each projection, in pixels
	float res; // resolution of pixels
	int num_proj;		// total number of projections
	float *angles;		// array containing the angles of each projection (degrees)
	float *z_shift;		// z_shift of projection
	unsigned short *detector;		// array of detector heads
	unsigned short* angular_view;	// index for each projection angle used to load the correct attenuation map
	unsigned short *timeslot;
	float *frame_time;	// in minutes

	float *_data;
	float ***data;		// projection data array of floats [projection][row (axial)][col (transaxial)]

	double AxialDetectorOffset[4];
	double TransaxialDetectorOffset[4];

	double AxialApertureOffset[4];
	double TransaxialApertureOffset[4];

	double RadiusOfRotation[4];

	char aperture[16];	// aperture name
};

/********************************************************************************************
/ Reconstruction class
*********************************************************************************************/
class Reconstruction
{
public:
	// constructors/destructor
	Reconstruction(int size_xy, int size_z, int timeslots, float new_res, char **files=NULL);
	Reconstruction(const char* filename);
	~Reconstruction();

	void Clear(int timeslot);	// erases the contents of data
	void Project(Projection* proj, const char* mmp_atten = NULL);
	void Backproject(const Projection* proj, const char* mmp_atten = NULL);
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

