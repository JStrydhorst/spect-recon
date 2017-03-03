// DicomObj.cpp

#include "DicomObj.h"

#include<cstring>
#include<sstream>

/******************************
/ Creates a new DicomObj from
/ the data supplied in the input
/ stream f
******************************/
DicomObj::DicomObj(std::ifstream &f)
{
	unsigned long BytesRead = 0;
	DataElement* newDataElement;
	char temp[4];
	const unsigned char DicomObjTag[4] = {0xFE,0xFF,0x00,0xE0};
	const unsigned char ItemEndTag[4] = {0xFE,0xFF,0x0D,0xE0};


	// read tags
	f.read(temp,4);
	if(std::memcmp(temp,DicomObjTag,4)) 
		std::cout << "Nested Dicom Object Tag error\n";
	// read ObjLength
	f.read((char*)&Length,4);

	if(Length == 0xFFFFFFFF)
	{
		// create data elements until an FFFE,E00D tag
		f.read(temp,4);
		while(std::memcmp(temp, ItemEndTag,4))
		{
			f.seekg(-4,std::ios::cur);
			newDataElement = new DataElement(f);
			AddObj(newDataElement);
			f.read(temp,4);
		}
		f.seekg(4,std::ios::cur);
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
			std::cout << "DicomObj length error\n";
	}
}

/***************************************
/ Calls the Print() function on each
/ object below itself in the heirarchy
***************************************/
void DicomObj::Print(std::ostream& os)
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
************************************/
void DicomObj::Write(std::ofstream& f)
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

void DicomObj::ReadFromFile(std::istream& in)
{
	unsigned short Group, Element;
	std::string VR;
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

	std::stringstream ss;
	unsigned int n_AT;
	unsigned short* AT;

	DataElement* newElement;
	DicomObj* SQObj;

	bool done;

	while(1)
	{
		// Read Group, Element, VS, and skip any whitespace
		in >> std::hex >> Group >> Element;
		// check for EOF or object delimiter
		if(in.eof() || Group==0xFFFE)
			break;

		in >> VR >> std::ws;

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
				ss >> std::hex >> AT[i];			
			ElementLen = 2*n_AT;
			newElement = new DataElement(Group,Element,VR,ElementLen,AT);
			SetElement(newElement);
			break;

		case VR_FL:
			in >> std::dec >> temp_fl;
			ElementLen  = 4;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_fl);
			SetElement(newElement);
			break;

		case VR_FD:
			in >> std::dec >> temp_dbl;
			ElementLen  = 8;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_dbl);
			SetElement(newElement);
			break;

		case VR_SL:
			in >> std::dec >> temp_lo;
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
				in >> std::hex >> SQ_Group >> SQ_Element;
				if(SQ_Group != 0xFFFE)
				{
					std::cout << "Error importing SQ Element" << std::endl;
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
					std::cout << "Error" << std::endl;
					done = true;
					break;
				}
			}

			SetElement(newElement);
			break;

		case VR_SS:
			in >> std::dec >> temp_sh;
			ElementLen  = 2;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_sh);
			SetElement(newElement);
			break;

		case VR_UL:
			in >> std::dec >> temp_ul;
			ElementLen  = 4;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_ul);
			SetElement(newElement);
			break;

		case VR_US:
			in >> std::dec >> temp_us;
			ElementLen  = 2;
			newElement = new DataElement(Group,Element,VR,ElementLen,&temp_us);
			SetElement(newElement);
			break;

		case VR_LT:
		case VR_ST:
		case VR_UT:
			//could contain CR/LF/FF/ESC, can't be imported by this function
			std::cout << "Import of LT/ST/UT data from text file not supported." << std::endl;
			break;

		case VR_OB:
		case VR_OF:
		case VR_OW:
		case VR_UN:
			std::cout << "Import of OB/OF/OW/UN elements from text file not supported." << std::endl;
			break;

		default:
			std::cout << "Unknown VR." << std::endl;
			break;
		}

	}

	CheckLength();

}


/***************************************
/ Calls the Validate() function on each
/ object below itself in the heirarchy
***************************************/
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
***************************************************/
int DicomObj::FindElement(unsigned short Group, unsigned short Element)
{
	FirstObj();

	if(!currentObj)
	{
		std::cout << "No elements in object." << std::endl;
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
**************************************************/
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
******************************************/
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
******************************************/
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
******************************************/
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
**************************************/
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
		return 0;
}
