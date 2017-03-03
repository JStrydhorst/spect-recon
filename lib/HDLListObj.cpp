// hdllist.h
// double linked list class

// Abstract base class for a heirarchical dynamic linked list

#include "HDLListObj.h"

// default construction
HDLListObj::HDLListObj()
{
	currentObj = NULL;
	nextObj = NULL;
	prevObj = NULL;
	numObj = 0;
}

HDLListObj::~HDLListObj()
{
	// deletes all objects underneath itself
	if(currentObj)
	{
		LastObj();
		while(DecrObjPtr())
			delete currentObj->nextObj;
		delete currentObj;
	}
}

void HDLListObj::FirstObj()
{
	if(currentObj)
		while(DecrObjPtr()){};
}

void HDLListObj::LastObj()
{
	if(currentObj)
		while(IncrObjPtr()){};
}

void HDLListObj::AddObj(HDLListObj* newListObj)
{
	LastObj();		// go to the end of the list
	if(currentObj)	// if this isn't the first element we're adding
	{
		currentObj->nextObj = newListObj;
		newListObj->prevObj = currentObj;
	}
	else
		newListObj->prevObj = NULL;

	newListObj->nextObj = NULL;	// the is now the last object in the list

	currentObj = newListObj;	// set the pointer to to point to the new object

	numObj++;
}

int HDLListObj::InsertObj(HDLListObj* newListObj)
{
	if(!currentObj)
	{
		std::cout << "Error: InsertObj, currentObj not valid";
		return -1;
	}

	newListObj->prevObj = currentObj->prevObj;
	newListObj->nextObj = currentObj;
	if(currentObj->prevObj)	// adding between two objects
		currentObj->prevObj->nextObj = newListObj;
	currentObj->prevObj = newListObj;

	currentObj = newListObj;

	numObj++;

	return 0;

}

void HDLListObj::DeleteObj()
{
	HDLListObj* temp;

	if(!currentObj)
		return;

	if(currentObj->prevObj) // if there's an object preceding it
		currentObj->prevObj->nextObj = currentObj->nextObj;

	if(currentObj->nextObj)	// if there's an object following it
	{
		currentObj->nextObj->prevObj = currentObj->prevObj;
		temp = currentObj->nextObj;
	}
	else
		temp = currentObj->prevObj;

	delete currentObj;
	
	numObj--;

	currentObj = temp;
}

HDLListObj* HDLListObj::ExtractObj()
{
	HDLListObj* temp;

	if(!currentObj)
		return NULL;

	temp = currentObj;

	// remove the current object from the set
	if(currentObj->prevObj) // if there's an object preceding this one
		currentObj->prevObj->nextObj = currentObj->nextObj;

	if(currentObj->nextObj)	// if there's an object after this one
	{
		currentObj->nextObj->prevObj = currentObj->prevObj;
		currentObj = currentObj->nextObj;
	}
	else
		currentObj = currentObj->prevObj;

	numObj--;

	// clear the next and prev pointers on the newly independent object
	temp->prevObj = NULL;
	temp->nextObj = NULL;
	return temp;
}

HDLListObj* HDLListObj::IncrObjPtr()
{
	if(currentObj->nextObj)
		return currentObj = currentObj->nextObj;
	else return NULL;
}

HDLListObj* HDLListObj::DecrObjPtr()
{
	if(currentObj->prevObj)
		return currentObj = currentObj->prevObj;
	else return NULL;
}

