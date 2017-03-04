// dicom.h
#include "RootDicomObj.h"

#include<cstring>

std::map<std::string, valrep_t> VRMap;

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

int RootDicomObj::depth = 0;

bool RootDicomObj::VRMap_Initialized = false;

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
	std::ifstream f;
	DataElement* newDataElement;
	char temp[4];
	unsigned short temp_us;

	f.open(filename,std::fstream::binary);
	if(!f.is_open())
	{
		std::cout << "File " << filename << " not found." << std::endl;
		return;
	}

	if(!VRMap_Initialized)
	{
		InitVRMap();
		VRMap_Initialized = true;
	}

	f.seekg(128,std::ios::beg); // first 128 bytes in a DICOM file are unused
	f.read(temp,4);
	// verify that the first 4 bytes are DICM
	if(std::memcmp(temp,"DICM",4))
	{
		std::cout << "Not a valid DICOM file!" << std::endl;
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
					f.seekg(-2,std::ios::cur);
			}
			f.seekg(-2,std::ios::cur);

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
******************************/
void RootDicomObj::Write(std::ofstream& f)
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
	std::ofstream fout;

	fout.open(filename,std::ios::binary);
	if(!fout.is_open())
		return -1;
	Write(fout);
	fout.close();

	return 0;
}

