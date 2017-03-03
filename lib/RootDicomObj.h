#ifndef __ROOTDICOMOBJ_H
#define __ROOTDICOMOBJ_H

#include "DicomObj.h"

class RootDicomObj : public DicomObj
{
public:
	RootDicomObj();
	RootDicomObj(const char* filename, bool hdr_only = false);	// constructor called for root Dicom objects
	~RootDicomObj() {};					// destructor

	// overrides the Write function from the DicomObj class
	void Write(std::ofstream& f);				// write Dicom object to file
	int Write(const char* filename);

	static int depth;					// only used for indenting when writing
	static bool VRMap_Initialized;		// tracks if the VR_Map has been intilized yet
};

#endif
