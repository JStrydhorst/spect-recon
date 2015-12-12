// dicom.h

// currently handles only expicit VR encoding
// implicit VR would require a library

#define _CRT_SECURE_NO_WARNINGS

#ifndef _DICOM_H
#define _DICOM_H

#include <cstdlib>
#include <cstring>

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "hdllist.h"

using namespace std;

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

static map<string, valrep_t> VRMap;

void InitVRMap()
{
	VRMap["AE"] = VR_AE;
	VRMap["AS"] = VR_AS;
	VRMap["AT"] = VR_AT;
	VRMap["CS"] = VR_CS;
	VRMap["DA"] = VR_DA;
	VRMap["DS"] = VR_DS;
	VRMap["DT"] = VR_DT;
	VRMap["FL"] = VR_FL;
	VRMap["FD"] = VR_FD;
	VRMap["IS"] = VR_IS;
	VRMap["LO"] = VR_LO;
	VRMap["LT"] = VR_LT;
	VRMap["OB"] = VR_OB;
	VRMap["OF"] = VR_OF;
	VRMap["OW"] = VR_OW;
	VRMap["PN"] = VR_PN;
	VRMap["SH"] = VR_SH;
	VRMap["SL"] = VR_SL;
	VRMap["SQ"] = VR_SQ;
	VRMap["SS"] = VR_SS;
	VRMap["ST"] = VR_ST;
	VRMap["TM"] = VR_TM;
	VRMap["UI"] = VR_UI;
	VRMap["UL"] = VR_UL;
	VRMap["UN"] = VR_UN;
	VRMap["US"] = VR_US;
	VRMap["UT"] = VR_UT;
}

class DataElement : public HDLListObj
{
public:
	DataElement(ifstream& f);	// creates a single data element from an ifstream
	DataElement(unsigned short newGroup, unsigned short newElement, string newVR,
				unsigned long newLength, const void* newValue);
	~DataElement();

	// virtual functions from abstract base class
	void Print(ostream& os=cout);			// displays the data contained in the DataElement
	void Write(ofstream& f);				// writes the DataElement to a file
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
	string VR;
	unsigned long Length;
	char* Value;
};

class DicomObj : public HDLListObj
{
public:
	DicomObj(){ Length = 0;};
	DicomObj(ifstream& f);				// constructor called for nested Dicom objects
	~DicomObj() {};

	// virtual functions from abstract base class
	void Print(ostream& os=cout);		// Prints the contents of the Object
	void Write(ofstream& f);			// write Dicom object to file with (FFFE, E0000) tag and length.
	void ReadFromFile(istream& f);		// reads and create elements from a text file;
	unsigned long GetSize();			// Length of the Dicom Object (including tags)
	int CheckLength();					// verifies that the length field is correct
	int Validate();

	int SetElement(DataElement*);	// inserts a data element into the Dicom Object (in the right place)
	int DeleteElement(unsigned short Group, unsigned short Element);
	int ModifyElement(unsigned short Group, unsigned short Element, void* newVal, unsigned long newLen);
	int FindElement(unsigned short Group, unsigned short Element);
	unsigned long GetLength(unsigned short Group, unsigned short Element);	// Length of the specified Element
	unsigned long GetValue(unsigned short Group, unsigned short Element, void* buff, unsigned long buf_size);

	DicomObj* GetSQObject(unsigned short Group, unsigned short Element, int n=0);

private:
	unsigned long Length;
};

class RootDicomObj : public DicomObj
{
public:
	RootDicomObj();
	RootDicomObj(const char* filename, bool hdr_only = false);	// constructor called for root Dicom objects
	~RootDicomObj() {};					// destructor

	// overrides the Write function from the DicomObj class
	void Write(ofstream& f);				// write Dicom object to file
	int Write(const char* filename);

	static int depth;					// only used for indenting when writing
	static bool VRMap_Initialized;		// tracks if the VR_Map has been intilized yet
};

int RootDicomObj::depth = 0;
bool RootDicomObj::VRMap_Initialized = false;

DataElement::DataElement(ifstream& f)
{
	char temp[4] = {0,0,0,0};
	unsigned char SQEndTag[] = {0xFE, 0xFF, 0xDD, 0xE0};

	DicomObj* newDCMObj;
	unsigned long BytesRead;

	f.read((char*)&Group,2);
	f.read((char*)&Element,2);

	f.read(temp,2);
	VR = temp;

	Length = 0;
	switch(VRMap[VR])
	{
		case VR_OB:		// if VR is OB, OW, OF, SQ, UT, or UN, skip two bytes, then read 4 byte length
		case VR_OW:
		case VR_OF:
		case VR_UT:
		case VR_UN:
			f.seekg(2,ios::cur);
			f.read((char*)&Length,4);
			break;
		case VR_SQ:
			Value = NULL;
			f.seekg(2,ios::cur);
			f.read((char*)&Length,4);
			// this is where we'll be creating new Dicom Objects for SQ data
			// for each Dicom Object in the sequence
			if(Length==0xFFFFFFFF)
			{
				// check the tag before passing to the DicomObj creator
				f.read(temp,4);
				while(memcmp(temp,SQEndTag,4))
				{
					f.seekg(-4,ios::cur);
					newDCMObj = new DicomObj(f);
					AddObj(newDCMObj);
					f.read(temp,4);
				}
				// skip the next four bytes and continue
				f.seekg(4,ios::cur);
			}
			else
			{
				BytesRead = 0;
				while(BytesRead < Length)
				{
					newDCMObj = new DicomObj(f);
					AddObj(newDCMObj);
					BytesRead += newDCMObj->GetSize();
				}
			}
			return;
		default:
			f.read((char*)&Length,2);
			break;
	}

	if(Length)
	{
		Value = new char[Length];
		if(Value)
			f.read(Value,Length);
		else
			cout << "Error: unable to allocate memory for data element" << endl;
	}
	else
		Value = NULL;

	return;
}

DataElement::DataElement(unsigned short newGroup, unsigned short newElement, string newVR,
				unsigned long newLength, const void* newValue)
{
	Group = newGroup;
	Element = newElement;
	VR = newVR;
	Length = newLength;

	if(newLength % 2)
		Length++;

	if(Length)
	{
		Value = new char[Length];
		if(Value)
			memcpy(Value, newValue, newLength);
		else
			return; // doesn't try to pad if memory allocation failed
	}
	else
		Value = NULL;

	if(newLength % 2) // pad odd field lengths with spaces or zeros
	{
		switch(VRMap[VR])
		{
		case VR_AE:
		case VR_AS:
		case VR_CS:
		case VR_DA:
		case VR_DS:
		case VR_DT:
		case VR_IS:
		case VR_LO:
		case VR_LT:
		case VR_PN:
		case VR_SH:
		case VR_ST:
		case VR_UT:
			Value[newLength] = 0x20;
			break;
		default:
			Value[newLength] = 0;
		}
	}
}


DataElement::~DataElement()
{
	if(Value)
		delete Value;
}


void DataElement::Print(ostream& os)
{
	for(int i=0;i<RootDicomObj::depth;i++)
		os << "    ";
	os << "(";
	os.width(4);
	os.fill('0');
	os << hex << Group << ",";
	os.width(4);
	os << Element << ")  ";

	os << VR;

	os.width(9);
	os.fill(' ');
	if(Length == 0xFFFFFFFF)
		os << "undef";
	else
		os << dec << Length << "  ";

	switch(VRMap[VR])
	{
	case VR_AE:
	case VR_AS:
	case VR_CS:
	case VR_DA:
	case VR_DS:
	case VR_DT:
	case VR_IS:
	case VR_LO:
	case VR_LT:
	case VR_PN:
	case VR_SH:
	case VR_ST:
	case VR_TM:
	case VR_UI:
	case VR_UT:
		if (Length < 80)
			os.write(Value, Length);
		else {
			os.write(Value, 80);
			os << "...";
		}
		break;
	case VR_AT:	// write the attribute tags as pairs of 16-bit values
		for(unsigned int i=0;i<(Length>>2);i++) {
			os << "(";
			os.width(4);
			os.fill('0');
			os << hex << *(unsigned short*)(Value + (i*4)) << ",";
			os.width(4);
			os << *(unsigned short*)(Value + (i*4 + 2)) << ") ";
		}
		break;
	case VR_FD:
		os << *(double*)Value;
		break;
	case VR_FL:
		os << *(float*)Value;
		break;
	case VR_SL:
		os << *(long*)Value;
		break;
	case VR_SS:
		os << *(short*)Value;
		break;
	case VR_UL:
		os << *(unsigned long*)Value;
		break;
	case VR_US:
		for(unsigned int i=0;i<(Length-2);i+=2)
			os << *(unsigned short*)(Value+i) << "\\";
		os << *(unsigned short*)(Value+Length-2);
		break;
	case VR_SQ:
		os << endl;

		// do some indentinig
		RootDicomObj::depth++;
		for(int i=0;i<RootDicomObj::depth;i++)
			os << "    ";
		os << "+++++";

		FirstObj();
		while(currentObj)
		{
			os << endl;
			((DicomObj*)currentObj)->Print(os);
			for(int i=0;i<RootDicomObj::depth;i++)
				os << "    ";
			os << "+++++";
			if(!IncrObjPtr())
				break;
		}
		// back out of the indenting
		RootDicomObj::depth--;
		break;
	default:	// OB, OF, OW, and UN don't get printed
		break;
	}
	os << endl;
}

// Writes the currect data element to a file
void DataElement::Write(ofstream &f)
{
	const unsigned char SQEndTag[] = {0xFE,0xFF,0xDD,0xE0};

	f.write((char*)&Group,2);
	f.write((char*)&Element,2);
	f << VR;

	switch(VRMap[VR])
	{
	case VR_OB:		// if VR is OB, OW, OF, SQ, UT, or UN, write two bytes, then write 4 byte length
	case VR_OW:
	case VR_OF:
	case VR_UT:
	case VR_UN:
	case VR_SQ:
		f.put(0);
		f.put(0);
		f.write((char*)&Length,4);
		break;
	default:
		f.write((char*)&Length,2);
		break;
	}
	
	switch(VRMap[VR])
	{
	case VR_SQ:
		FirstObj();
		while(currentObj)
		{
			currentObj->Write(f);
			if(!IncrObjPtr())
				break;
		}

		if(Length == 0xFFFFFFFF)
		{
			f.write((char*)SQEndTag, 4);
			f.put(0);
			f.put(0);
			f.put(0);
			f.put(0);
		}
		break;
	default:
		if(Length)
			f.write(Value,Length);
		break;
	}
}


// Returns the actual size of the element when written to a file,
// including all tags. SQs with unspecified lengths have to be added up
// but SQs with specified lengths return the stored Length. If an
// underlying data structure has changed, CheckLength should be called
// first on the RootDicomObject to ensure that the value stored in the
// Length fields of the SQ elements are correct
unsigned long DataElement::GetSize()
{
	unsigned long Size = 0;
	
	switch(VRMap[VR])
	{
		case VR_OB:		// if VR is OB, OW, OF, SQ, UT, or UN, skip two bytes, then read 4 byte length
		case VR_OW:
		case VR_OF:
		case VR_UT:
		case VR_UN:
			return Length + 12;
		case VR_SQ:
			if(Length == 0xFFFFFFFF)
			{
				// read the size of the underlying elements
				FirstObj();
				while(currentObj)
				{
					Size += currentObj->GetSize();
					if(!IncrObjPtr())
						break;
				}
				return Size + 20;
			}
			else
				return Length + 12;
		default:
			return Length + 8;
	}
}

// copies the contents of the Value field into buffer
// returns Length (number of bytes copied)
unsigned long DataElement::GetValue(void* buffer, unsigned long buf_size)
{
	if (VRMap[VR] == VR_SQ)	// SQ data isn't stored in value
		return -1;

	if (buf_size < Length)
		return 0;

	memcpy(buffer, Value, Length);

	return Length;

}

// modifies the contents of the Value field (and sometimes the Length)
// returns 0 if sucessful, -1 otherwise
int DataElement::Modify(void* newVal, unsigned long newLen)
{
	if(Length == newLen)
	{
		memcpy(Value, newVal, newLen);
		return 0;
	}
	else
	{
		delete Value;

		Length = newLen + (newLen % 2); // pad to an even number of bytes if necessary
		Value = new char[Length];
		if(Value)
		{
			memcpy(Value, newVal, newLen);	// copy data

			if(newLen%2)					// pad if needed
			{
				switch(VRMap[VR])
				{
				case VR_AE:
				case VR_AS:
				case VR_CS:
				case VR_DA:
				case VR_DS:
				case VR_DT:
				case VR_IS:
				case VR_LO:
				case VR_LT:
				case VR_PN:
				case VR_SH:
				case VR_ST:
				case VR_UT:
					Value[newLen] = 0x20;
					break;
				default:
					Value[newLen] = 0;
				}
			}
		}
		else
			return -1; // unable to allocate memory for Value
	}

	return 0;
}

// Checks the lengths of the SQ DataElements by
// adding up the lengths of the underlying
// returns 1 if a length had to be fixed, 0 otherwise
// Does nothing if the DataElement isn't an SQ
int DataElement::CheckLength()
{
	unsigned long Size = 0;
	int changes = 0;

	// verifies the lengths of SQ objects
	if(VRMap[VR]==VR_SQ)
	{
		// call CheckLength on all underlying objects
		FirstObj();
		while(currentObj)
		{
			changes += currentObj->CheckLength(); // called on all DicomObj under the SQ DataElement
			if(!IncrObjPtr())
				break;
		}

		if(Length==0xFFFFFFFF)
			return changes;	// don't modify the length field if the length is undefined
		
		// verify the length by adding up the lengths of nested DicomObjects
		FirstObj();
		while(currentObj)
		{
			Size += currentObj->GetSize(); // called on DicomObjects
			if(!IncrObjPtr())
				break;
		}

		if(Size != Length)
		{
			// cout << "Updating length of (" << hex << Group << "," << hex << Element << ")" << endl;
			Length = Size;
			changes++;
		}

	}
	return changes;
}


int DataElement::SetObject(HDLListObj * newObj)
{
	if(VRMap[VR] != VR_SQ)
		return -1;	// objects can only be inserted into SQ elements

	AddObj(newObj);
	return 0;

}


// checks that the data element is valid
// returns 0 if OK, -1 if there is an error
int DataElement::Validate()
{
	int i;
	int len;
	int VM = 0;	// value multiplicity
	char* data;
	char* pch;

	bool invalid_data = false;

	switch(VRMap[VR])
	{
	case VR_AE:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if(len>16)
				invalid_data = true;
			while(pch[len-1]==' ') // trim trailing spaces
				len--;
			i=0;
			while(pch[i]==' ')	// trim leading spaces
				i++;
			// check for invalid characters
			for(;i<len;i++)
				if(pch[i]<' ' || pch[i]>'~')
					invalid_data = true;
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_AS:
		VM = 1;
		if(Length!=4)
			return -1;
		for(i=0;i<3;i++)
			if(Value[i]<'0' || Value[i]>'9')
				invalid_data = true;
		if(!(Value[3]=='D' || Value[3]=='W' || Value[3]=='M' || Value[3]=='Y'))
			invalid_data = true;
		break;

	case VR_AT:
		if(Length%4)
			invalid_data = true;
		VM = Length/4;
		break;

	case VR_CS:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if(len>16)
				invalid_data = true;
			while(pch[len-1]==' ') // trim trailing spaces
				len--;
			i=0;
			while(pch[i]==' ')	// trim leading spaces
				i++;
			// check for invalid characters
			for(;i<len;i++)
			{
				if(pch[i]==' ' || pch[i]=='_')
					continue;
				if(pch[i]>='0' && pch[i]<='9')
					continue;
				if(pch[i]>='A' && pch[i]<='Z')
					continue;
				invalid_data = true;
			}
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_DA:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if((len!=8) && (len!=10))
				invalid_data = true;
			// check for invalid characters
			while(pch[len-1]==' ') // trim trailing spaces
				len--;
			for(i=0;i<len;i++)
			{
				if(pch[i]>='0' && pch[i]<='9')
					continue;
				if(pch[i]=='.')
					continue;
				invalid_data = true;
			}
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_DS:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if(len>16)
				invalid_data = true;
			// check for invalid characters
			while(pch[len-1]==' ') // trim trailing spaces
				len--;
			i=0;
			while(pch[i]==' ')	// trim leading spaces
				i++;
			for(;i<len;i++)
			{
				if(pch[i]>='0' && pch[i]<='9')
					continue;
				if(pch[i]=='.' || pch[i]=='+' || pch[i]=='-')
					continue;
				if(pch[i]=='e' || pch[i]=='E')
					continue;
				invalid_data = true;
			}
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_DT:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if(len>26)
				invalid_data = true;
			while(pch[len-1]==' ') // trim trailing spaces
				len--;
			for(i=0;i<len;i++)
			{
				if(pch[i]>='0' && pch[i]<='9')
					continue;
				if(pch[i]=='.' || pch[i]=='+' || pch[i]=='-')
					continue;
				invalid_data = true;
			}
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_FL:
		if(Length%4)
			invalid_data = true;
		VM = Length/4;
		break;

	case VR_FD:
		if(Length%8)
			invalid_data = true;
		VM = Length/8;
		break;

	case VR_IS:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if(len>12)
				invalid_data = true;
			// check for invalid characters
			while(pch[len-1]==' ') // trim trailing spaces
				len--;
			i=0;
			while(pch[i]==' ')	// trim leading spaces
				i++;
			for(;i<len;i++)
			{
				if(pch[i]>='0' && pch[i]<='9')
					continue;
				if(pch[i]=='+' || pch[i]=='-')
					continue;
				invalid_data = true;
			}
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_LO:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if(len>64)
				invalid_data = true;
			// check for invalid characters
			for(i=0;i<len;i++)
			{
				if(pch[i]>=' ' && pch[i]<='~')
					continue;
				if(pch[i]=='\x1b')
					continue;
				invalid_data = true;
			}
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_LT:	// this validation doesn't support extended character sets
		VM=1;
		if(Length>10240)
			invalid_data = true;
		for(i=0;i<int(Length);i++)
		{
			if(Value[i]>=' ' && Value[i]<='~')
				continue;
			if(Value[i]=='\x0a' || Value[i]=='\x0c' || Value[i]=='\x0d' || Value[i]=='\x1b')
				continue;
			invalid_data = true;
		}
		break;

	case VR_OB:
	case VR_OF:
	case VR_OW:
		break;

	case VR_PN:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if(len>64)
				invalid_data = true;
			while(pch[len-1]==' ') // trim trailing spaces
				len--;
			// check for invalid characters
			for(i=0;i<len;i++)
			{
				if(pch[i]>=' ' && pch[i]<='~')
					continue;
				if(pch[i]=='\x1b')
					continue;
				invalid_data = true;
			}
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_SH:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if(len>16)
				invalid_data = true;
			// check for invalid characters
			for(i=0;i<len;i++)
			{
				if(pch[i]>=' ' && pch[i]<='~')
					continue;
				if(pch[i]=='\x1b')
					continue;
				invalid_data = true;
			}
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_SL:
		if(Length%4)
			invalid_data = true;
		VM = Length/4;
		break;

	case VR_SQ:
		// call validate on each object in the sequence
		FirstObj();
		while(currentObj)
		{
			if(currentObj->Validate()==-1)
				invalid_data = true;
			if(!IncrObjPtr())
				break;
		}
		break;

	case VR_SS:
		if(Length%2)
			invalid_data = true;
		VM = Length/2;
		break;

	case VR_ST:	// this validation doesn't support extended character sets
		VM=1;
		if(Length>1024)
			invalid_data = true;
		for(i=0;i<int(Length);i++)
		{
			if(Value[i]>=' ' && Value[i]<='~')
				continue;
			if(Value[i]=='\x0a' || Value[i]=='\x0c' || Value[i]=='\x0d' || Value[i]=='\x1b')
				continue;
			invalid_data = true;
		}
		break;

	case VR_TM:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if(len>16)
				invalid_data = true;
			while(pch[len-1]==' ') // trim trailing spaces
				len--;
			// check for invalid characters
			for(i=0;i<len;i++)
			{
				if(pch[i]>='0' && pch[i]<='9')
					continue;
				if(pch[i]=='.')
					continue;
				invalid_data = true;
			}
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_UI:
		data = new char[Length+2];
		strncpy(data,Value,Length);
		data[Length] = 0;
		pch = strtok(data,"\\");
		while(pch!=NULL)
		{
			len = strlen(pch);
			// check length
			if(len>64)
				invalid_data = true;
			// check for invalid characters
			for(i=0;i<len;i++)
			{
				if(pch[i]>='0' && pch[i]<='9')
					continue;
				if(pch[i]=='.')
					continue;
				invalid_data = true;
			}
			pch = strtok(NULL,"\\");
			VM++;
		}
		delete [] data;
		break;

	case VR_UL:
		if(Length%4)
			invalid_data = true;
		VM = Length/4;
		break;

	case VR_UN:
		break;

	case VR_US:
		if(Length%2)
			invalid_data = true;
		VM = Length/2;
		break;

	case VR_UT:	// this validation doesn't support extended character sets
		VM=1;
		for(unsigned int i=0;i<Length;i++)
		{
			if(Value[i]>=' ' && Value[i]<='~')
				continue;
			if(Value[i]=='\x0a' || Value[i]=='\x0c' || Value[i]=='\x0d' || Value[i]=='\x1b')
				continue;
			invalid_data = true;
		}
		break;
		
	default:
		break;

	}

	if(invalid_data)
	{
		cout.fill('0');
		cout << "Inavlid Data in (" << hex << setw(4) << Group << "," << hex << setw(4) << Element << ")" << endl;
		return -1;
	}

	return 0;
}

unsigned long DataElement::GetLength()
{
	CheckLength();
	return Length;
}

HDLListObj* DataElement::GetSQObject(int n)
{
	if(VRMap[VR] == VR_SQ)
	{
		FirstObj();

		while(n)
		{
			if(!IncrObjPtr())
				return NULL;	// nth element doesn't exist
			n--;
		}

		return currentObj;
	}
	else
		return NULL;
}

RootDicomObj::RootDicomObj()
{
	if(!VRMap_Initialized)
	{
		InitVRMap();
		VRMap_Initialized = true;
	}
}

RootDicomObj::RootDicomObj(const char* filename, bool hdr_only)
{
	ifstream f;
	DataElement* newDataElement;
	char temp[4];
	unsigned short temp_us;

	f.open(filename,fstream::binary);

	if(!VRMap_Initialized)
	{
		InitVRMap();
		VRMap_Initialized = true;
	}

	f.seekg(128,ios::beg); // first 128 bytes in a DICOM file are unused
	f.read(temp,4);
	// verify that the first 4 bytes are DICM
	if(memcmp(temp,"DICM",4))
	{
		cout << "Not a valid DICOM file!" << endl;
		f.close();
		return;
	}

	while(f.peek() != EOF) // process the entire file
	{
		if(hdr_only)
		{
			// check next tag, skip 07FE,0010
			f.read(reinterpret_cast<char*>(&temp_us),2);
			if(temp_us == 0x7FE0)
			{
				f.read(reinterpret_cast<char*>(&temp_us),2);
				if(temp_us == 0x0010)
				{
					f.close();
					return;
				}
				else
					f.seekg(-2,ios::cur);
			}
			f.seekg(-2,ios::cur);

		}
		newDataElement = new DataElement(f);
		AddObj(newDataElement);
	}

	CheckLength();

	f.close();
}


/******************************
/ Writes a DICOM file to disk
/ Always returns 0.
/******************************/
void RootDicomObj::Write(ofstream& f)
{
	CheckLength();

	for(int i=0; i<128; i++)
		f.put(0);
	f.write("DICM",4);

	FirstObj();
	while(currentObj)
	{
		currentObj->Write(f);
		if(!IncrObjPtr())
			break;
	}
}

// wrapper for the write funtion that takes a filename
int RootDicomObj::Write(const char* filename)
{
	ofstream fout;

	fout.open(filename,ios::binary);
	if(!fout.is_open())
		return -1;
	Write(fout);
	fout.close();

	return 0;
}

/******************************
/ Creates a new DicomObj from
/ the data supplied in the input
/ stream f
/******************************/
DicomObj::DicomObj(ifstream &f)
{
	unsigned long BytesRead = 0;
	DataElement* newDataElement;
	char temp[4];
	const unsigned char DicomObjTag[4] = {0xFE,0xFF,0x00,0xE0};
	const unsigned char ItemEndTag[4] = {0xFE,0xFF,0x0D,0xE0};


	// read tags
	f.read(temp,4);
	if(memcmp(temp,DicomObjTag,4)) 
		cout << "Nested Dicom Object Tag error\n";
	// read ObjLength
	f.read((char*)&Length,4);

	if(Length == 0xFFFFFFFF)
	{
		// create data elements until an FFFE,E00D tag
		f.read(temp,4);
		while(memcmp(temp, ItemEndTag,4))
		{
			f.seekg(-4,ios::cur);
			newDataElement = new DataElement(f);
			AddObj(newDataElement);
			f.read(temp,4);
		}
		f.seekg(4,ios::cur);
	}
	else
	{
		while(BytesRead < Length)
		{
			newDataElement = new DataElement(f);
			AddObj(newDataElement);
			BytesRead += newDataElement->GetSize();
		}

		if(BytesRead != Length)
			cout << "DicomObj length error\n";
	}
}

/***************************************
/ Calls the Print() function on each
/ object below itself in the heirarchy
/***************************************/
void DicomObj::Print(ostream& os)
{
	FirstObj();
	while(currentObj)
	{
		currentObj->Print(os);
		if(!IncrObjPtr())
			break;
	}
}

/************************************
/ Writes a DicomObj to the file specified
/ Writes the Item Tag, then the length,
/ then writes the elements below itself.
/ If the length is unspecified, the Item
/ Delimitation Tag is written after the
/ Element list.
/************************************/
void DicomObj::Write(ofstream& f)
{
	const unsigned char DicomObjTag[4] = {0xFE,0xFF,0x00,0xE0};
	const unsigned char ItemEndTag[4] = {0xFE,0xFF,0x0D,0xE0};

	f.write((char*)DicomObjTag,4);
	f.write((char*)&Length,4);

	FirstObj();
	while(currentObj)
	{
		currentObj->Write(f);
		if(!IncrObjPtr())
			break;
	}

	if(Length == 0xFFFFFFFF)
	{
		f.write((char*)ItemEndTag,4);
		f.put(0);
		f.put(0);
		f.put(0);
		f.put(0);
	}
}

void DicomObj::ReadFromFile(istream& in)
{
	unsigned short Group, Element;
	string VR;
	unsigned long ElementLen;
	char Data[1024];
	//char* temp_str;	// used for long fields of text
	float temp_fl;
	double temp_dbl;
	long temp_lo;
	unsigned long temp_ul;
	short temp_sh;
	unsigned short temp_us;
	unsigned short SQ_Group, SQ_Element;

	stringstream ss;
	unsigned int n_AT;
	unsigned short* AT;

	DataElement* newElement;
	DicomObj* SQObj;

	bool done;

	while(1)
	{
		// Read Group, Element, VS, and skip any whitespace
		in >> hex >> Group >> Element;
		// check for EOF or object delimiter
		if(in.eof() || Group==0xFFFE)
			break;

		in >> VR >> ws;

		switch(VRMap[VR])
		{
		case VR_AE:
		case VR_AS:
		case VR_CS:
		case VR_DA:
		case VR_DS:
		case VR_DT:
		case VR_IS:
		case VR_LO:
		case VR_PN:
		case VR_SH:
		case VR_TM:
		case VR_UI:
			in.getline(Data,1024);
			ElementLen = strlen(Data);
			newElement = new DataElement(Group,Element,VR,ElementLen,Data);
			SetElement(newElement);
			break;

		case VR_AT:
			// read data
			in.getline(Data,1024);
			// count attribute tags and set up storage
			n_AT = ((unsigned int)(in.gcount())+1)/5;
			ss << Data;
			AT = new unsigned short[n_AT];
			for(unsigned int i=0;i<n_AT;i++)
				ss >> hex >> AT[i];			
			ElementLen = 2*n_AT;
			newElement = new DataElement(Group,Element,VR,ElementLen,AT);
			SetElement(newElement);
			break;

		case VR_FL:
			in >> dec >> temp_fl;
			ElementLen  = 4;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_fl);
			SetElement(newElement);
			break;

		case VR_FD:
			in >> dec >> temp_dbl;
			ElementLen  = 8;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_dbl);
			SetElement(newElement);
			break;

		case VR_SL:
			in >> dec >> temp_lo;
			ElementLen  = 4;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_lo);
			SetElement(newElement);
			break;

		case VR_SQ:
			// Create SQ Element
			newElement = new DataElement(Group,Element,VR,0,NULL);
			done = false;
			while(!done)
			{
				in >> hex >> SQ_Group >> SQ_Element;
				if(SQ_Group != 0xFFFE)
				{
					cout << "Error importing SQ Element" << endl;
					done = true;
					break;
				}

				switch(SQ_Element)
				{
				case 0xE000:
					SQObj = new DicomObj();
					SQObj->ReadFromFile(in);
					newElement->SetObject(SQObj);
					break;

				case 0xE0DD:
					done = true;
					break;

				default:
					cout << "Error" << endl;
					done = true;
					break;
				}
			}

			SetElement(newElement);
			break;

		case VR_SS:
			in >> dec >> temp_sh;
			ElementLen  = 2;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_sh);
			SetElement(newElement);
			break;

		case VR_UL:
			in >> dec >> temp_ul;
			ElementLen  = 4;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_ul);
			SetElement(newElement);
			break;

		case VR_US:
			in >> dec >> temp_us;
			ElementLen  = 2;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_us);
			SetElement(newElement);
			break;

		case VR_LT:
		case VR_ST:
		case VR_UT:
			//could contain CR/LF/FF/ESC, can't be imported by this function
			cout << "Import of LT/ST/UT data from text file not supported." << endl;
			break;

		case VR_OB:
		case VR_OF:
		case VR_OW:
		case VR_UN:
			cout << "Import of OB/OF/OW/UN elements from text file not supported." << endl;
			break;

		default:
			cout << "Unknown VR." << endl;
			break;
		}

	}

	CheckLength();

}


/***************************************
/ Calls the Validate() function on each
/ object below itself in the heirarchy
/***************************************/
int DicomObj::Validate()
{
	bool invalid_obj=false;

	FirstObj();
	while(currentObj)
	{
		if(currentObj->Validate()==-1)
			invalid_obj = true;
		if(!IncrObjPtr())
			break;
	}

	if(invalid_obj)
		return -1;

	return 0;
}
/***************************************************
/ Finds and sets the currentObj pointer to point to
/ the element at the next level of the heirarcy that
/ mathches the group and element specified
/ returns 1 if the object is found, returns 0 if
/ if the object is not found.
/***************************************************/
int DicomObj::FindElement(unsigned short Group, unsigned short Element)
{
	FirstObj();

	if(!currentObj)
	{
		cout << "No elements in object." << endl;
		return 0;
	}

	// Increment the current object pointer until we're at the right group
	while(dynamic_cast<DataElement*>(currentObj)->GetGroup() < Group)
		if(!IncrObjPtr())
			return 0;	// reached the end of the Dicom Object

	// If we're at the right group, increment the current object pointer until we're at the right element
	while( (dynamic_cast<DataElement*>(currentObj)->GetGroup() == Group) && (dynamic_cast<DataElement*>(currentObj)->GetElement() < Element) )
		if(!IncrObjPtr())
			return 0;	// reached the end of the Dicom Object

	// If we've found the object, return 1
	if( (dynamic_cast<DataElement*>(currentObj)->GetGroup() == Group) && (dynamic_cast<DataElement*>(currentObj)->GetElement() == Element) )
		return 1;

	// We've past the point where it should be and we haven't found if yet. Return 0.
	return 0;
}

/**************************************************
/ Inserts a new DataElement into a Dicom Item
/ If the element already exists, the existing
/ element is deleted first. Always returns 0.
/**************************************************/
int DicomObj::SetElement(DataElement* newElement)
{
	// if no objects exist yet
	if(!currentObj)
	{
		AddObj(newElement);
		return 0;
	}

	// the object is already there and we're replacing it
	if(FindElement(newElement->GetGroup(), newElement->GetElement()))
	{
		DeleteObj();
		InsertObj(newElement);
		return 0;
	}
	
	// the object is not there and we've reached the end of the list
	// (the current Group is less than the new group, but the Find function returned 0)
	if( dynamic_cast<DataElement*>(currentObj)->GetGroup() < newElement->GetGroup() )
	{
		AddObj(newElement);
		return 0;
	}

	// (or the Group is correct but the current element is less than the new one)
	if( (dynamic_cast<DataElement*>(currentObj)->GetGroup() == newElement->GetGroup()) 
		   && (dynamic_cast<DataElement*>(currentObj)->GetElement() < newElement->GetElement()) )
	{
		AddObj(newElement);
		return 0;
	}

	InsertObj(newElement);
	return 0;
}

/******************************************
/ Finds and deletes the specified element
/******************************************/
int DicomObj::DeleteElement(unsigned short Group, unsigned short Element)
{
	if(FindElement(Group, Element))
	{
		DeleteObj();
		return 0;
	}
	else
		return -1;
}

/******************************************
/ Finds the element specified, then calls
/ returns the length of the object
/******************************************/
unsigned long DicomObj::GetLength(unsigned short Group, unsigned short Element)
{
	// find the element
	if(FindElement(Group, Element))
		return dynamic_cast<DataElement*>(currentObj)->GetLength();	// copy the data into the buffer
	else
		return -1;	// (gggg,eeee) pair not found
}

/******************************************
/ Finds the element specified, then calls
/ the function to copy the data from the
/ element into the buffer specified. Returns
/ the number of bytes copied, or -1 if the
/ DataElement isn't found.
/******************************************/
unsigned long DicomObj::GetValue(unsigned short Group, unsigned short Element, void* buffer, unsigned long buf_size)
{
	// find the element
	if(FindElement(Group, Element))
		return dynamic_cast<DataElement*>(currentObj)->GetValue(buffer, buf_size);	// copy the data into the buffer
	else
		return -1;	// (gggg,eeee) pair not found
}

/**************************************
/ Returns the amount of space occupied
/ on disk by the DicomObj, including
/ the tags.
/ CheckLength should be called to update
/ the Length field before calling GetSize
/ if any objects or elements have changed
/**************************************/
unsigned long DicomObj::GetSize()
{
	unsigned long Size = 0;
	if(Length == 0xFFFFFFFF)
	{
		FirstObj();
		while(currentObj)
		{
			Size+= currentObj->GetSize();
			if(!IncrObjPtr())
				break;
		}
		return Size + 16;
	}
	else
		return Length + 8;
}

// change the value stored in a data element
// returns 0 if things are OK, -1 if there's an error (Element not found or not modified sucessfully)
int DicomObj::ModifyElement(unsigned short Group, unsigned short Element, void *newVal, unsigned long newLen)
{
	if(FindElement(Group, Element))
		return dynamic_cast<DataElement*>(currentObj)->Modify(newVal, newLen);  // Modify returns 0 unless there's a problem

	return -1;
}

// Checks the length of the DICOM object by calling GetSize on all Data Elements
// returns 0 if all Length fields are correct, returns 1 if the Length were
// updated. This should be called before saving a DicomObject to a file to ensure
// that all the SQ elements have the correct length field.
int DicomObj::CheckLength()
{
	unsigned long Size = 0;
	int changes = 0;

	// call CheckLength on all lower levels of heirarchy first
	FirstObj();
	while(currentObj)
	{
		changes += currentObj->CheckLength(); // called on DataElements
		if(!IncrObjPtr())
			break;
	}

	FirstObj();
	while(currentObj)
	{
		Size += currentObj->GetSize(); // called on DataElements
		if(!IncrObjPtr())
			break;
	}

	if(Size != Length)
	{
		// cout << "Updating length of DicomObj." << endl;
		Length = Size;
		changes++;
	}

	return changes;
}

DicomObj* DicomObj::GetSQObject(unsigned short Group, unsigned short Element, int n)
{
	if(FindElement(Group,Element))
	{
		return dynamic_cast<DicomObj*>(dynamic_cast<DataElement*>(currentObj)->GetSQObject(n));
	}
	else
		return NULL;
}

#endif

#undef _CRT_SECURE_NO_WARNINGS