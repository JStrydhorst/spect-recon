#include <cstdlib>
#include <iostream>
#include <fstream>

#include "proj_table.h"

using namespace std;

ProjTable::ProjTable(char* filename)
{
	// initially allocate space for 2^20 (~1M) entries
	if(filename)
		LoadProjTable(filename);
	else
	{
		n = 0;
		mem_alloc = 1<<20;
		try
		{
			ri = new uchar[mem_alloc];
			rj = new uchar[mem_alloc];
			rk = new uchar[mem_alloc];
			pi = new uchar[mem_alloc];
			pj = new uchar[mem_alloc];
			weight = new float[mem_alloc];
		}
		catch(bad_alloc& ba)
		{
			cout << "ProjTable(char*): bad_alloc caught: " << ba.what() << endl;
		}
	}
}

ProjTable::~ProjTable()
{
	delete [] ri;
	delete [] rj;
	delete [] rk;
	delete [] pi;
	delete [] pj;
	delete [] weight;
}

void ProjTable::AddEntry(uchar new_ri, uchar new_rj, uchar new_rk, uchar new_pi, uchar new_pj, float new_weight)
{
	unsigned int i;

	uchar* temp_uchar_p;
	float* temp_float_p;

	if(n>=mem_alloc)
	{
		// allocate some more memory first
		mem_alloc += 1<<20;

		try
		{

			temp_uchar_p = new uchar[mem_alloc];
			for(i=0;i<n;i++)
				temp_uchar_p[i] = ri[i];
			delete [] ri;
			ri = temp_uchar_p;

			temp_uchar_p = new uchar[mem_alloc];
			for(i=0;i<n;i++)
				temp_uchar_p[i] = rj[i];
			delete [] rj;
			rj = temp_uchar_p;

			temp_uchar_p = new uchar[mem_alloc];
			for(i=0;i<n;i++)
				temp_uchar_p[i] = rk[i];
			delete [] rk;
			rk = temp_uchar_p;

			temp_uchar_p = new uchar[mem_alloc];
			for(i=0;i<n;i++)
				temp_uchar_p[i] = pi[i];
			delete [] pi;
			pi = temp_uchar_p;

			temp_uchar_p = new uchar[mem_alloc];
			for(i=0;i<n;i++)
				temp_uchar_p[i] = pj[i];
			delete [] pj;
			pj = temp_uchar_p;

			temp_float_p = new float[mem_alloc];
			for(i=0;i<n;i++)
				temp_float_p[i] = weight[i];
			delete [] weight;
			weight = temp_float_p;
		}
		catch(bad_alloc& ba)
		{
			cout << "ProjTable::AddEntry(): bad_alloc caught: "<< ba.what() << endl;
			return;
		}

	}
	
	// add new entry
	ri[n] = new_ri;
	rj[n] = new_rj;
	rk[n] = new_rk;
	pi[n] = new_pi;
	pj[n] = new_pj;
	weight[n] = new_weight;
	n++;
}


void ProjTable::SaveProjTable(char* filename)
{
	ofstream f;

	f.open(filename,ios::binary);

	f.write(reinterpret_cast<char*>(&n),sizeof(unsigned int));
	f.write(reinterpret_cast<char*>(ri),n*sizeof(unsigned char));
	f.write(reinterpret_cast<char*>(rj),n*sizeof(unsigned char));
	f.write(reinterpret_cast<char*>(rk),n*sizeof(unsigned char));
	f.write(reinterpret_cast<char*>(pi),n*sizeof(unsigned char));
	f.write(reinterpret_cast<char*>(pj),n*sizeof(unsigned char));
	f.write(reinterpret_cast<char*>(weight),n*sizeof(float));

	f.close();
	
}

int ProjTable::LoadProjTable(char* filename)
{

	unsigned int i;

	ifstream f;

	f.open(filename,ios::binary);

	f.read(reinterpret_cast<char*>(&n),sizeof(unsigned int));

	mem_alloc = n;
	try
	{
		ri = new uchar[n];
		rj = new uchar[n];
		rk = new uchar[n];
		pi = new uchar[n];
		pj = new uchar[n];
		weight = new float[n];
	}
	catch(bad_alloc&)
	{
		cout << "Error allocating memory." << endl;
		return -1;
	}

	f.read(reinterpret_cast<char*>(ri),n*sizeof(unsigned char));
	f.read(reinterpret_cast<char*>(rj),n*sizeof(unsigned char));
	f.read(reinterpret_cast<char*>(rk),n*sizeof(unsigned char));
	f.read(reinterpret_cast<char*>(pi),n*sizeof(unsigned char));
	f.read(reinterpret_cast<char*>(pj),n*sizeof(unsigned char));
	f.read(reinterpret_cast<char*>(weight),n*sizeof(float));

	f.close();
	
	for(i=0;i<n;i++)
		weight[i] /= 800;	// kludge to fix projection table scale

	return 0;

}

void ProjTable::ProjTableStats()
{
	unsigned int i;
	float min, max;
	min = weight[0];
	max = weight[0];

	cout << "Size of table: " << n << endl;
	for(i=0;i<n;i++)
	{
		if (weight[i] < min) min = weight[i];
		if (weight[i] > max) max = weight[i];
		//if (weight[i] > 1.5)
			//cout << int(ri[i]) << "," << int(rj[i]) << "," << int(rk[i]) << " --> " << int(pi[i]) << "," << int(pj[i]) << "     " << weight[i] << endl;
	}

	cout << "Min value: " << min << endl;
	cout << "Max value: " << max << endl;

	cout << endl;

}

ProjTable* ProjTable::Abridge(float threshold)
{
	unsigned int i;
	ProjTable *newTable = new ProjTable();

	for (i=0;i<n;i++)
		if(weight[i]>threshold)
			newTable->AddEntry(ri[i],rj[i],rk[i],pi[i],pj[i],weight[i]);

	return newTable;
}