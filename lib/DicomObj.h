#ifndef __DICOMOBJ_H
#define __DICOMOBJ_H

#include "HDLListObj.h"
#include "DataElement.h"

#include<iostream>
#include<fstream>

class DicomObj : public HDLListObj
{
public:
	DicomObj(){ Length = 0;};
	DicomObj(std::ifstream& f);				// constructor called for nested Dicom objects
	~DicomObj() {};

	// virtual functions from abstract base class
	void Print(std::ostream& os=std::cout);		// Prints the contents of the Object
	void Write(std::ofstream& f);			// write Dicom object to file with (FFFE, E0000) tag and length.
	void ReadFromFile(std::istream& f);		// reads and create elements from a text file;
	unsigned long GetSize();			// Length of the Dicom Object (including tags)
	int CheckLength();					// verifies that the length field is correct
	int Validate();

	int SetElement(DataElement*);	// inserts a data element into the Dicom Object (in the right place)
	int DeleteElement(unsigned short Group, unsigned short Element);
	int ModifyElement(unsigned short Group, unsigned short Element, void* newVal, unsigned long newLen);
	int FindElement(unsigned short Group, unsigned short Element);
	uint32_t GetLength(unsigned short Group, unsigned short Element);	// Length of the specified Element
	uint32_t GetValue(unsigned short Group, unsigned short Element, void* buff, unsigned long buf_size);

	DicomObj* GetSQObject(unsigned short Group, unsigned short Element, int n=0);

private:
	unsigned long Length;
};

#endif // __DICOMOBJ_H
