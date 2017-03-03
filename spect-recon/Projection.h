#ifndef __PROJECTION_H
#define __PROJECTION_H

/********************************************************************************************
/ Projection class
*********************************************************************************************/
class Projection
{
public:
	Projection();
	Projection(int projections, int size_y, int size_z, float res, char* aperture_name, const char* filename = 0); // blank or binary file
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

	void resort(int*, int);
};

#endif // __PROJECTION_H
