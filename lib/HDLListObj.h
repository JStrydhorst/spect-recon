#ifndef __HDLLISTOBJ_H
#define __HDLLISTOBJ_H

#include<iostream>

class HDLListObj
{
public:
	HDLListObj();
	virtual ~HDLListObj();

//	int CountObj() { return numObj;} // returns the number of objects

	// virtual functions to be implemented in the derived classes
	virtual void Print(std::ostream& os=std::cout) = 0;	// Prints to object to a specified stream
	virtual void Write(std::ofstream&) = 0;			// Writes the object to a file
	virtual unsigned long GetSize() = 0;
	virtual int CheckLength() = 0;
	virtual int Validate() = 0;

protected:
	void FirstObj();			// sets the currentObj pointer to the first object
	void LastObj();				// sets the currentObj pointer to the last object
	HDLListObj* IncrObjPtr();	// increments the currentObj to the next object
	HDLListObj* DecrObjPtr();	// returns 0 if at the end or at the beginning
	void AddObj(HDLListObj*);	// adds an object to the end of the list
	int InsertObj(HDLListObj*);	// inserts object before the currentObj
	void DeleteObj();			// deletes the current object
	HDLListObj* ExtractObj();	// removes the current object from the list and returns the pointer to it

	HDLListObj* currentObj;		// pointer to object at next level of heirarchy

private:
	HDLListObj* nextObj;	// next object at the same level
	HDLListObj* prevObj;	// previous object at the same level
	
	int numObj;			// number of elements at next level of heirarchy
};

#endif // __HDLLISTOBJ_H
