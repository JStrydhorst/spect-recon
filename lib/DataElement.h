#ifndef __DATAELEMENT_H
#define __DATAELEMENT_H

#include "HDLListObj.h"

#include<string>
#include<fstream>
#include<map>

enum valrep_t {VR_AE,
                                VR_AS,
                                VR_AT,
                                VR_CS,
                                VR_DA,
                                VR_DS,
                                VR_DT,
                                VR_FL,
                                VR_FD,
                                VR_IS,
                                VR_LO,
                                VR_LT,
                                VR_OB,
                                VR_OF,
                                VR_OW,
                                VR_PN,
                                VR_SH,
                                VR_SL,
                                VR_SQ,
                                VR_SS,
                                VR_ST,
                                VR_TM,
                                VR_UI,
                                VR_UL,
                                VR_UN,
                                VR_US,
                                VR_UT};

extern std::map<std::string, valrep_t> VRMap;


class DataElement : public HDLListObj
{
public:
	DataElement(std::ifstream& f);	// creates a single data element from an ifstream
	DataElement(unsigned short newGroup, unsigned short newElement, std::string newVR,
				unsigned long newLength, const void* newValue);
	~DataElement();

	// virtual functions from abstract base class
	void Print(std::ostream& os=std::cout);			// displays the data contained in the DataElement
	void Write(std::ofstream& f);				// writes the DataElement to a file
	unsigned long GetSize();
	int CheckLength();
	int Validate();							// checks the element against a Dicom dictionary

	int Modify(void* newVal, unsigned long newLen);

	int SetObject(HDLListObj*);						// inserts Dicom Objects into SQ objects

	unsigned short GetGroup() { return Group; }
	unsigned short GetElement() { return Element; }
	unsigned long GetLength();
	unsigned long GetValue(void* buffer, unsigned long buf_size); // copies the Value to buffer and returns Length

	HDLListObj* GetSQObject(int n);		// returns the DicomObject in the SQ DataElement indexed by n;
	// int GetNumSQObj();	// returns the total number of sequenced objects is function needed?

private:
	unsigned short Group;
	unsigned short Element;
	std::string VR;
	unsigned long Length;
	char* Value;
};

#endif // __DATAELEMENT_H
