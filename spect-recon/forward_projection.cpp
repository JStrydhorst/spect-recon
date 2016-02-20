#define _CRT_SECURE_NO_WARNINGS

#include <cstdlib>
#include <iostream>

#include "recon4d.h"

using namespace std;

#ifdef _WIN32
#define __SL "\\"
#define __HOME "c:"
#elif defined __unix
#define __SL "/"
#define __HOME getenv("HOME")
#endif

int main(int argc, char* argv[])
{
//	int t;

	char working_dir[260];
	char sysmat_dir[260];
	char filename[260];
//	char** source_files;

	int recon_z = 128;
	int recon_xy = 128;
	float recon_res = 0.6f;

	int n_timeslots = 1;

	//Projection* proj_30;
	//Projection* proj_20;
	//Projection* proj_10;

	// parse command line
	sprintf(working_dir,"%s", argv[1]);

	// create projection data
	sprintf(filename, "%s"__SL"SPECT.dcm",working_dir);	// default data file is SPECT.dcm; should verify existence
	Projection* proj = new Projection(filename);

	recon_z = int( proj->GetScanLength() / recon_res) + recon_xy + 1;	// extended to allow rotation to ends, assumes system matrix is a cube
	// create and load reconstruction
	/*
	source_files = new char*[n_timeslots];
	for(t=0;t<n_timeslots;t++)
	{
		source_files[t] = new char[260];
		sprintf(source_files[t], "%s"__SL"%s_%d.bin",working_dir, argv[2],t+1);
	}*/

	Reconstruction* recon = new Reconstruction(recon_xy,recon_z,recon_res,"c:\\SPECT\\2011_11_Abs_in_vivo\\rat1\\rat_volume.bin");		// reconstruction
	recon->SaveBIN("c:\\SPECT\\2011_11_Abs_in_vivo\\rat1\\rat_volume_padded.bin");
	/*
	for(t=0;t<n_timeslots;t++)
		delete [] source_files[t];
	delete [] source_files;
	*/

	// load system matrix
	sprintf(sysmat_dir,"%s"__SL"SPECT"__SL"system_matrix",__HOME);	// default pinhole directory
	cout << "System matrix directory: " << sysmat_dir << endl << endl;
	cout << "Loading projection tables..." << endl << endl;
	if(recon->LoadProjTables(proj,sysmat_dir) == -1)
	{
		cout << "Error loading projection tables." << endl;
		return -1;
	}

	// proj->AddTimeslots(n_timeslots);

	// project noise free data
	//for(int t=0;t<n_timeslots;t++)
	recon->Project(proj);
	sprintf_s(filename,"%s"__SL"projections.dcm",working_dir)
	recon->SaveDICOM(filename);

	recon->UnloadProjTables();
	delete recon;
	/*
	// save
	sprintf(filename,"%s"__SL"proj_noise_free.bin",working_dir);
	proj->SaveBIN(filename);
	sprintf(filename,"%s"__SL"proj_noise_free.dcm",working_dir);
	proj->SaveDICOM(filename,true);


	proj_30 = new Projection(*proj);
	proj_30->Scale(10.0f);	// 
//	proj_30->ApplyDecay(360.0f);
	proj_30->Poisson(2373468483l);	// add poisson noise
	sprintf(filename,"c:\\SPECT\\Xavier\\proj_noise_10.dcm");
	proj_30->SaveDICOM(filename,false);
	delete proj_30;

	delete proj;
	exit(0);

	proj_20 = new Projection(*proj);
	proj->Scale(0.66f);	// 
//	proj->ApplyDecay(360.0f);
	proj_20->Poisson(81962908);	// add poisson noise
	sprintf(filename,"c:\\SPECT\\MOBY\\1.0mm\\20_min\\SPECT.dcm");
	proj_20->SaveDICOM(filename,false);
	delete proj_20;

	proj_10 = new Projection(*proj);
	proj->Scale(0.33f);	// 
//	proj->ApplyDecay(360.0f);
	proj_10->Poisson(7523467);	// add poisson noise
	sprintf(filename,"c:\\SPECT\\MOBY\\1.0mm\\10_min\\SPECT.dcm");
	proj_10->SaveDICOM(filename,false);
	delete proj_10;
	*/

	delete proj;

}
