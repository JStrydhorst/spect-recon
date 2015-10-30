#define _CRT_SECURE_NO_WARNINGS

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <ctime>
#include <iostream>
#include <fstream>

#include "recon4d.h"

using namespace std;

#ifdef _WIN32
#define __SL "\\"
#define __HOME "c:"
#elif defined __unix
#define __SL "/"
#define __HOME getenv("HOME")
#endif

#define SPECT_OSEM_CPP_VER 1

int main(int argc, char* argv[])
{
	int arg;
	char filename[260], proj_file[260];
	bool atten_corr = false;
	bool scatter = false;
	bool create_atten_map = false;
	bool save_raw_data = false;
	bool save_current_recon = false;
	bool log_voi_activity = false;
	float pixel_limit = 0.0f;			// when set to zero, no limit is imposed on the pixel values
	float filter_fwhm = 0.0f, filter_sigma=0.0f;
	float update_threshold = 0.1f;		// pixels where the normalization falls below the threshold are not updated, default 0.1
	float half_life = 360.35f;			// default is the halflife (in minutes) of Tc-99m. Should obtain from DICOM file.
	float beta_s = 0.0f, delta=1.0f;			// weighting and parameter of spatial MAP prior
	float beta_t = 0.0f;				// weighting of temporal MAP prior
	float beta_CT = 0.0f; // delta_CT = 20.0f;				// CT prior, delta_CT defines the difference required for a CT boundary
	float atten_scale = 1.0f;			// for scaling the attenuation map

	time_t rawtime;
	tm* time_info;
	char dt_string[10];		// string for naming files and writing info
	ofstream recon_info;

	char working_dir[260];
	char sysmat_dir[260];
	char atten_map[260];

	char** file_list = new char*[1];

	// reconstruction parameters
	int recon_z;
	int recon_xy = 128;
	float recon_res = 0.6f;

	float scatter_k = 0.5f;

	// iteration and counters
	int t;
	int n_timeslots;
	int group, iter, subset;
	int n_groups = 4;
	int n_subsets[] = {6,4,3,1};
	int n_iter = 8;

	ifstream f;

	// Projection and recon data
	Projection* ProjectionData;		// original projection set file
	Projection** TrueProj;			// projection data divided up into time slots

	Reconstruction* OSEM_Recon;		// current reconstructions
	Reconstruction* OSEM_TempRecon;	// backprojection of ratios (only one, with one timeframe)
	Reconstruction** OSEM_Norms;	// normalization (one for each subset, only one timeframe)
//	Reconstruction* MLEM_Norm;		// normalization with all projections (only one of these)
	Reconstruction* MAP_Partials;	// partial derivatives used in MAP reconstruction
	Reconstruction* CT_Prior;

	Projection*** OSEM_Subsets;		// true projection data split up into subsets
	Projection** OSEM_TempProj;		// projection of current estimate
	Projection* NormProj;			// projection for creatings normalization
	Projection* Scatter;			// raw scatter data
	Projection** GatedScatter;		// gated scatter data
	Projection*** Scatter_Subsets;	// scatter chopped up into subsets

	cout << endl;


	// Set defaults and parse command line
	sprintf(working_dir,"%s" __SL "SPECT" __SL "temp",__HOME);			// default reconstruction directory
	sprintf(sysmat_dir,"%s" __SL "SPECT" __SL "system_matrix",__HOME);	// default pinhole directory

	if (argc>1)
	{
		arg = 1;
		if(strncmp(argv[1],"-",1)!=0)	// change working directory if the first argument isn't a switch
		{
			sprintf(working_dir,"%s", argv[1]);
			arg++;
		}

		while(arg<argc)
		{
			if(strcmp(argv[arg],"-a")==0)
				atten_corr = true;
			else if(strcmp(argv[arg],"-s")==0)
			{
				scatter = true;
				if((arg+1)<argc) // don't check for a k value if -s is the last argument
				{
					if(isdigit(argv[arg+1][0]) || (argv[arg+1][0]==0x2E)) // check if a non-default value for k has been entered
					{
						arg++;
						scatter_k = float(atof(argv[arg]));
						// check that scatter_k is reasonable?
						if (scatter_k == 0.0f)
						{
							cout << "Invalid scatter coefficient\n";
							scatter = false;
						}
					}
				}

			}
			else if(strcmp(argv[arg],"-m") == 0)	// use specified system matrix instead of default
			{
				arg++;
				if(arg<argc)
					sprintf(sysmat_dir,"%s",argv[arg]);

				// should check here if directory supplied is valid
			}
			else if(strcmp(argv[arg],"-l") == 0)	// limit pixel values
			{
				arg++;
				if(arg<argc)
					pixel_limit = float(atof(argv[arg]));
			}
			else if(strcmp(argv[arg],"-f") == 0)	// apply filtering
			{
				arg++;
				filter_fwhm = float(atof(argv[arg]));	// filtering is specified in terms of resolution (in mm)
				filter_sigma = filter_fwhm/2.35482f;	// convert to sigma
			}
			else if(strcmp(argv[arg],"-A")==0)		// create attenuation map instead of reconstruction
			{
				create_atten_map = true;
				arg++;
				sprintf(atten_map,"%s" __SL "%s",working_dir,argv[arg]);
			}
			else if(strcmp(argv[arg],"-t")==0)		// specify the update threshold
			{
				arg++;
				update_threshold = float(atof(argv[arg]));
			}
			else if(strcmp(argv[arg],"-b") == 0)	// include spatial smoothing MAP prior
			{
				arg++;
				beta_s = float(atof(argv[arg]));
				arg++;
				delta = float(atof(argv[arg]));
			}
			else if(strcmp(argv[arg],"-x") == 0)	// include spatial smoothing MAP with CT prior
			{
				arg++;
				beta_CT = float(atof(argv[arg]));
				arg++;
				delta = float(atof(argv[arg]));
				arg++;
				sprintf(atten_map,"%s" __SL "%s",working_dir,argv[arg]);
			}
			else if(strcmp(argv[arg],"-g") == 0)	// include temporal MAP prior
			{
				arg++;
				beta_t = float(atof(argv[arg]));
			}
			else if(strcmp(argv[arg],"-h") == 0)	// use different half-life
			{
				arg++;
				half_life = float(atof(argv[arg]));
			}
			else if(strcmp(argv[arg],"-B") == 0)	// scale the attenuation coefficient to use broad-beam attenuation
			{
				arg++;
				atten_scale = float(atof(argv[arg]));
			}
			else if(strcmp(argv[arg],"-low_res") == 0)	// use the low resolution reconstruction algorithm
			{
				recon_res = 1.0f;
				recon_xy = 76;
			}
			else if(strcmp(argv[arg],"-r") == 0)	// save the raw data in bin files
				save_raw_data = true;
			else if(strcmp(argv[arg],"-c") == 0)	// save the raw data in bin files
				save_current_recon = true;
			else if(strcmp(argv[arg],"-v") == 0)	// log VOI activity, VOI shape and location is hardcoded in recon4d.cpp
				log_voi_activity = true;
			else
			{
				cout << "Unsupported option `" << argv[arg] << "'" << endl;
				return 0;
			}

			arg++;
		}
	}

	// Create info file
	time(&rawtime);
	time_info = localtime(&rawtime);
	strftime(dt_string,10,"%m%d%H%M",time_info);
	sprintf(filename,"%s" __SL "recon_info_%s.txt",working_dir,dt_string);
	recon_info.open(filename);

	cout << "SPECT4D v" << SPECT4D_MAJOR_VERSION << "." << SPECT4D_MINOR_VERSION << "." << SPECT_OSEM_CPP_VER << endl;
	cout << "Working directory: " << working_dir << endl << endl;

	recon_info << "SPECT4D v" << SPECT4D_MAJOR_VERSION << "." << SPECT4D_MINOR_VERSION << "." << SPECT_OSEM_CPP_VER << endl;
	recon_info << "Working directory: " << working_dir << endl << endl;

	char *atten_dir;
	if(atten_corr | create_atten_map)
	{
		atten_dir = new char[260];

		if(atten_scale == 1.0f)
			sprintf(atten_dir, "%s" __SL "atten", working_dir);
		else
			sprintf(atten_dir, "%s" __SL "atten_%.2f", working_dir, atten_scale);
	}
	else
		atten_dir = NULL;

	sprintf(proj_file, "%s" __SL "SPECT.dcm",working_dir);	// default data file is SPECT.dcm; should verify existence
	ProjectionData = new Projection(proj_file);
	ProjectionData->ProjectionInfo();
	recon_info << "Projection Data:" << proj_file << endl;

	n_timeslots = ProjectionData->GetNumTimeSlots();
	TrueProj = ProjectionData->SplitTimeSlots();

	delete ProjectionData;

	OSEM_Subsets = new Projection**[n_timeslots];
	Scatter_Subsets = new Projection**[n_timeslots];

	recon_z = int( TrueProj[0]->GetScanLength() / recon_res) + recon_xy + 1;	// extended to allow rotation to ends, assumes system matrix is a cube


	if(create_atten_map)	// only uses data from first timeslot for projection angles, etc...
	{
		file_list[0] = atten_map;
		OSEM_Recon = new Reconstruction(recon_xy, recon_z, 1, recon_res, file_list);		// reconstruction
		OSEM_Recon->ReconstructionInfo();
	
		cout << "Creating attenuation maps..." << endl;

		sprintf(filename,"%s" __SL "%s",sysmat_dir,"APT2.pin");
	
		OSEM_Recon->CreateMMPAttenMap(atten_dir, filename, TrueProj[0], atten_scale);

		recon_info << "Created attenuation maps in: " << atten_dir << endl;
		recon_info.close();

		delete TrueProj;
		delete OSEM_Recon;

		return 0;
	}

	OSEM_Recon = new Reconstruction(recon_xy,recon_z,n_timeslots,recon_res);		// create reconstruction
	OSEM_Recon->ReconstructionInfo();


	recon_info << "Reconstruction dimensions: " << recon_xy << " x " << recon_xy << " x " << recon_z << ", " << recon_res << " mm voxels" << endl;
	if(n_timeslots>1)
		recon_info << "Time Slots: " << n_timeslots << endl;

	// Physics corrections
	if(atten_corr)
	{
		cout << "Attenuation correction enabled." << endl;
		recon_info << "Attenuation correction enabled." << endl;
		if(atten_scale != 1.0f)
			recon_info << "Using broad-beam attenuation, " << atten_scale << endl;
	}
	if(scatter)
	{
		cout << "Scatter compensation enabled. (k=" << scatter_k << ")" << endl;
		recon_info << "Scatter compensation enabled. (k=" << scatter_k << ")" << endl;
	}
	if(scatter && !atten_corr)
		cout << "Scatter compensation without attenuation correction is sort of dumb." << endl;

	// Recon algorithm options (filtering, voxel limits, support thresholds)
	if(filter_sigma > 0.0f)
	{
		cout << "Post filtering enabled: FWHM=" << filter_fwhm << " mm (sigma=" << filter_sigma << ")" << endl;
		recon_info << "Post filtering enabled: FWHM=" << filter_fwhm << " mm (sigma=" << filter_sigma << ")" << endl;
	}
	if(pixel_limit > 0.0f)
	{
		cout << "Voxel values limited to: " << pixel_limit << endl;
		recon_info << "Voxel values limited to: " << pixel_limit << endl;
	}
	if(update_threshold > 0.0f)
	{
		cout << "Voxels with support < " << update_threshold << " ignored." << endl;
		recon_info << "Voxels with support < " << update_threshold << " ignored." << endl;
	}
	
	// MAP and CT prior info
	if(beta_s > 0.0f)
	{
		cout << "Spatial smoothing prior, beta=" << beta_s << ", delta=" << delta << endl;
		recon_info << "Spatial smoothing prior, beta=" << beta_s << ", delta=" << delta << endl;
	}
	if(beta_CT > 0.0f)
	{
		cout << "Using CT prior: " << atten_map << ", beta=" << beta_s << ", delta=" << delta << "." << endl;
		recon_info << "Using CT prior: " << atten_map << ", beta=" << beta_s << ", delta=" << delta << "." << endl;
	}
	if(beta_t > 0.0f)
	{
		cout << "Temporal smoothing prior, beta=" << beta_t << ", delta=" << delta << endl;
		recon_info << "Temporal smoothing prior, beta=" << beta_t << endl;
	}


	OSEM_TempRecon = new Reconstruction(recon_xy,recon_z,1,recon_res);	// temp recon only need one
	
	// Load projection tables
	cout << "System matrix directory: " << sysmat_dir << endl << endl;
	cout << "Loading projection tables..." << endl << endl;
	if(OSEM_Recon->LoadProjTables(TrueProj[0],sysmat_dir) == -1)
	{
		cout << "Error loading projection tables." << endl;
		return -1;
	}

	if(beta_s > 0.0f || beta_t > 0.0f)
	{
		MAP_Partials = new Reconstruction(recon_xy,recon_z,1,recon_res);
	}
	else
		MAP_Partials = NULL;

	if(beta_CT > 0.0f)
	{
		file_list[0] = atten_map;
		CT_Prior = new Reconstruction(recon_xy,recon_z,1,recon_res,file_list);
	}
	else
		CT_Prior = NULL;


	// load scatter data
	if(scatter)
	{
		sprintf(filename, "%s" __SL "%s",working_dir,"scatter.dcm");
		Scatter = new Projection(filename);
		// smooth the scatter projections
		Scatter->Filter(3.0f);
		// divide into subsets
		GatedScatter = Scatter->SplitTimeSlots();  //TODO: Should maybe combine scatter data to get better statistics? Before filtering?
		delete Scatter;
	}

	// start with backprojected norm map to get rid of regions with no activity
	if(atten_corr)
		if(atten_scale==1.0f)
			sprintf(filename, "%s" __SL "norms" __SL "norm_with_att_1_0.bin", working_dir);
		else
			sprintf(filename, "%s" __SL "norms" __SL "norm_with_att_1_0_%.2f.bin", working_dir, atten_scale);
	else
			sprintf(filename, "%s" __SL "norms" __SL "norm_no_att_1_0.bin", working_dir);

	f.open(filename,ios::in);
	if(f.is_open())
		f.close();
	else
	{
		NormProj = TrueProj[0]->CreateNorm(half_life);
		OSEM_TempRecon->Backproject(NormProj,atten_dir);
		delete NormProj;
		OSEM_TempRecon->SaveBIN(filename);
	}
	OSEM_Recon->LoadBIN(filename, true);

	for(group = 0; group < n_groups; group++)
	{
		cout << "Group: " << group+1 << endl;
		// split the projections up into subsets
		for(t=0;t<n_timeslots;t++)
			OSEM_Subsets[t] = TrueProj[t]->CreateOrderedSubsets(n_subsets[group]);
		// create a projection for the forward projections (projections hold the shift and rotation information, stays the same over time slots, not over subsets)
		OSEM_TempProj = TrueProj[0]->CreateOrderedSubsets(n_subsets[group]);

		if(scatter)
			for(t=0;t<n_timeslots;t++)
				Scatter_Subsets[t] = GatedScatter[t]->CreateOrderedSubsets(n_subsets[group]);
				
		// create or load the normalization
		cout << "Generating normalization matrices..." << endl;
		OSEM_Norms = new Reconstruction*[n_subsets[group]];
		for(subset=0;subset<n_subsets[group];subset++)
		{
			// check if normalizations exist
			if(atten_corr)
				if(atten_scale==1.0f)
					sprintf(filename, "%s" __SL "norms" __SL "norm_with_att_%d_%d.bin", working_dir, n_subsets[group], subset);
				else
					sprintf(filename, "%s" __SL "norms" __SL "norm_with_att_%d_%d_%.2f.bin", working_dir, n_subsets[group], subset, atten_scale);
			else
				sprintf(filename, "%s" __SL "norms" __SL "norm_no_att_%d_%d.bin", working_dir, n_subsets[group], subset);

			f.open(filename,ios::in);

			if(f.is_open())
			{
				f.close();
				cout << ".";
				file_list[0] = atten_map;
				OSEM_Norms[subset] = new Reconstruction(recon_xy,recon_z,1,recon_res,file_list);
			}
			else
			{
				NormProj = OSEM_Subsets[0][subset]->CreateNorm(half_life); //, !dig_phantom);
				OSEM_Norms[subset] = new Reconstruction(recon_xy,recon_z,1,recon_res);
				OSEM_Norms[subset]->Backproject(NormProj, atten_dir);
				OSEM_Norms[subset]->SaveBIN(filename);	// save normalization to speed up recon later
				delete NormProj;
			}
		}
		cout << endl;

		// main reconstruction loop
		for(iter=0;iter<n_iter;iter++)
		{

			cout << "Iteration: " << iter+1 << endl;
			for(subset=0;subset<n_subsets[group];subset++)
			{
				cout << "Subset: " << subset+1 << endl;
				for(t=0;t<n_timeslots;t++)
				{
					if(n_timeslots>1)
						cout << "t=" << t << endl;
					OSEM_TempProj[subset]->SetTimeslot(t+1);	// timeslot numbering in projections starts with 1
					OSEM_Recon->Project(OSEM_TempProj[subset], atten_dir);
					//OSEM_TempProj[subset]->ApplyOffsets(0); these corrections might be applied by acquisition software. Hard to tell.

					// apply intrisinc camera resolution here?

					if(scatter)
						OSEM_TempProj[subset]->Sum(Scatter_Subsets[t][subset],scatter_k);
				
					OSEM_TempProj[subset]->Ratio(OSEM_Subsets[t][subset]);
					//OSEM_TempProj[subset]->ApplyOffsets(1);

					// apply camera resolution again?

					OSEM_TempRecon->Backproject(OSEM_TempProj[subset],atten_dir);

					if(MAP_Partials)
					{
						MAP_Partials->Clear(0);
						if((beta_s > 0.0f) && (iter > 0))
							MAP_Partials->UpdateMAPSpatialPartial(OSEM_Recon,t,beta_s);
						if((beta_t > 0.0f) && (iter > 0))
							MAP_Partials->UpdateMAPTemporalPartial(OSEM_Recon,t,beta_t);
					}

					OSEM_Recon->Update(OSEM_TempRecon, OSEM_Norms[subset], MAP_Partials, t,update_threshold);

					//save binary of current update
					if(save_current_recon)
					{
						sprintf(filename, "%s" __SL "current_recon.bin",working_dir);
						OSEM_Recon->SaveBIN(filename);
					}
#ifdef VOI_LOG
					// log activity in VOI
					if(log_voi_activity)
						recon_info << group << ", "<< iter << ", " << subset << ", " << t << ": " << OSEM_Recon->VOIActivity() << endl;
#endif
				}
			}
			cout << endl;

			//OSEM_Recon->CleanUpStrayPixels(MLEM_Norm, pixel_limit);
			if (pixel_limit > 0.0f)
				OSEM_Recon->CleanUpStrayPixels(pixel_limit);

			if (filter_sigma > 0.0f)
			{
				OSEM_Recon->Filter(filter_sigma);
				// mask out voxels outside of the system matrix, since some of these get activity as a result of the filtering operation
				if(recon_res==0.6)
					OSEM_Recon->Mask(63);
				else if(recon_res==0.1)
					OSEM_Recon->Mask(37);
			}
		}

	
		// delete normalization data
		for(subset=0;subset<n_subsets[group];subset++)
			delete OSEM_Norms[subset];
		delete [] OSEM_Norms;
		
		// delete projection data
		for(t=0;t<n_timeslots;t++)
		{
			for(subset=0;subset<n_subsets[group];subset++)
				delete OSEM_Subsets[t][subset];
			delete [] OSEM_Subsets[t];
		}

		for(subset=0;subset<n_subsets[group];subset++)
				delete OSEM_TempProj[subset];
		delete [] OSEM_TempProj;

		// delete scatter data
		if(scatter)
		{
			for(t=0;t<n_timeslots;t++)
			{
				for(subset=0;subset<n_subsets[group];subset++)
					delete Scatter_Subsets[t][subset];
				delete [] Scatter_Subsets[t];
			}
		}

	}

	// save the final results
	if(save_raw_data)
	{
		sprintf(filename, "%s" __SL "recons" __SL "OSEM_%s.bin", working_dir, dt_string);
		OSEM_Recon->SaveBIN(filename,52);
	}

	sprintf(filename, "%s" __SL "recons" __SL "OSEM_%s.dcm", working_dir, dt_string);
	OSEM_Recon->SaveDICOM(filename,proj_file,52);


	recon_info.close();

		if(scatter)
		{
			for(t=0;t<n_timeslots;t++)
				delete GatedScatter[t];
			delete [] Scatter_Subsets;
		}
		delete [] OSEM_Subsets;

	if(MAP_Partials)
		delete MAP_Partials;

//	delete MLEM_Norm;
	delete TrueProj;
	delete OSEM_TempRecon;
	delete OSEM_Recon;

	OSEM_Recon->UnloadProjTables();

	if(atten_dir)
		delete [] atten_dir;
	
}

