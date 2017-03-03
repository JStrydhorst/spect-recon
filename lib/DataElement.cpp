// DataElement.cpp

// currently handles only expicit VR encoding
// implicit VR would require a library

#include "DataElement.h"
#include "DicomObj.h"
#include "RootDicomObj.h"

#include<cstring>
#include<iomanip>

DataElement::DataElement(std::ifstream& f)
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
			f.seekg(2,std::ios::cur);
			f.read((char*)&Length,4);
			break;
		case VR_SQ:
			Value = NULL;
			f.seekg(2,std::ios::cur);
			f.read((char*)&Length,4);
			// this is where we'll be creating new Dicom Objects for SQ data
			// for each Dicom Object in the sequence
			if(Length==0xFFFFFFFF)
			{
				// check the tag before passing to the DicomObj creator
				f.read(temp,4);
				while(std::memcmp(temp,SQEndTag,4))
				{
					f.seekg(-4,std::ios::cur);
					newDCMObj = new DicomObj(f);
					AddObj(newDCMObj);
					f.read(temp,4);
				}
				// skip the next four bytes and continue
				f.seekg(4,std::ios::cur);
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
			std::cout << "Error: unable to allocate memory for data element" << std::endl;
	}
	else
		Value = NULL;

	return;
}

DataElement::DataElement(unsigned short newGroup, unsigned short newElement, std::string newVR,
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
		delete [] Value;
}


void DataElement::Print(std::ostream& os)
{
	for(int i=0;i<RootDicomObj::depth;i++)
		os << "    ";
	os << "(";
	os.width(4);
	os.fill('0');
	os << std::hex << Group << ",";
	os.width(4);
	os << Element << ")  ";

	os << VR;

	os.width(9);
	os.fill(' ');
	if(Length == 0xFFFFFFFF)
		os << "undef";
	else
		os << std::dec << Length << "  ";

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
			os << std::hex << *(unsigned short*)(Value + (i*4)) << ",";
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
		os << std::endl;

		// do some indentinig
		RootDicomObj::depth++;
		for(int i=0;i<RootDicomObj::depth;i++)
			os << "    ";
		os << "+++++";

		FirstObj();
		while(currentObj)
		{
			os << std::endl;
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
	os << std::endl;
}

// Writes the currect data element to a file
void DataElement::Write(std::ofstream &f)
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
		std::cout.fill('0');
		std::cout << "Inavlid Data in (" << std::hex << std::setw(4) << Group << ",";
		std::cout << std::hex << std::setw(4) << Element << ")" << std::endl;
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


