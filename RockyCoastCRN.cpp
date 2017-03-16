/*------------------------------------------------------------------------

	RockyCoastCRN.cpp
	
	Object to evolve the CRN concentration across a coastal profile as a 
	function of tides, topographic shielding, relative sea level rise, 
	cliff retreat and platform downwear. For comparison with measured
	10Be CRN measurements across a coastal platform.
	
	Please see/cite the following publications:
	
	Hurst, M.D., Rood, D.H., Ellis, M.A., Anderson, R.S. and 
	Dornbusch, U., (2016) Recent acceleration in coastal cliff retreat rates
	on the south coast of Great Britain. PNAS 113(47), 13336-13341,
	http://dx.doi.org/10.1073/pnas.1613044113

	Hurst, M. D., Rood, D. H., and Ellis, M. A. (2017) Controls on the 
	distribution of cosmogenic 10Be across shore platforms. Earth Surf. 
	Dynam. 5, 67-84, http://dx.doi.org/10.5194/esurf-5-67-2017
	
	Martin D. Hurst, British Geological Survey, December 2014

	Copyright (C) 2015, Martin Hurst
	
	
	Updates to allow compatibility with Hiro Matsumoto's shore platform model
	which has been coded up and added to RoBoCoP as part of MASTS project
	
	Matsumoto, H., Dickson, M. E., & Kench, P. S. (2016)
	An exploratory numerical model of rocky shore profile evolution. 
	Geomorphology, 268, 98-109, http://dx.doi.org/10.1016/j.geomorph.2016.05.017
	
	March 2017
	
	Copyright (C) 2017, Martin Hurst
	
	Developer can be contacted
	martin.hurst@glasgow.ac.uk
  
	Martin D. Hurst
	School of Geographical and Earth Sciences
	University of Glasgow
	Glasgow
	Scotland
	G12 8QQ
  
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

------------------------------------------------------------------------*/

#ifndef RockyCoastCRN_CPP
#define RockyCoastCRN_CPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <omp.h>
#include "RockyCoastCRN.hpp"
#include "Hiro.hpp"

using namespace std;

void RockyCoastCRN::Initialise()
{
	cout << "Warning: You have not specified which nuclides you wish to track" << endl;
	cout << "Defaulting to 10Be only" << endl;
	
	vector<int> WhichNuclides(1,10);
	Initialise(WhichNuclides);
}

void RockyCoastCRN::Initialise(vector<int> WhichNuclides)
{
  /* initialise an empty platform object
	retreatrate is the rate of cliff retreat (m/yr)
	beachwidth is the width of the beach (m), protecting the platform from CRN production
	JunctionElevation is the elevation of the platform/cliff junction
	platformgradient is the average slope of the platform surface (m/m) 
	cliffheight is the height of the adjacent cliff (m)
	tidalamplitude is the average tidal amplitude for diurnal tides. */
	
	cout << "Warning: You have initialised an empty RockCoastCRN object." << endl;
	cout << "This is only appropriate when reading in platform morphology from a file" << endl;
	
	NDV = -9999;
	RetreatRate1 = NDV;
	RetreatRate2 = NDV;
	ChangeTime = NDV;
	PlatformGradient = NDV;
	CliffHeight = NDV;
	BeachWidth = NDV;
	BermHeight = NDV;
	JunctionElevation = NDV;
	TidalAmplitude = NDV;
	
	XMin = -1;
	XMax = -1;
	double ZMax = 5;
  	double ZMin = -20;
  	
  	//Cliff is at the start of the vector
	CliffHeight = 10.;
	CliffPositionX = -1;
	CliffPositionInd = 0;
		
	//Setup CRN domain
	dX = 1.0;
	dZ = 0.1;
	dt = 1;
	
	NXNodes = (XMax-XMin)/dX + 1;
	NZNodes = (ZMax-ZMin)/dZ + 1;
	
	//Initialise the nuclides
	InitialiseNuclides(WhichNuclides);
	
	//Setup Surface Arrays
	vector<double> EmptyZ(NZNodes,0);
	vector<double> EmptyX(NXNodes,0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	X = EmptyX;
	Z = EmptyZ;
	PlatformElevation = EmptyXNDV;
	PlatformElevationOld = EmptyXNDV;
	SurfaceElevation = EmptyXNDV;
	vector< vector <double> > EmptySurfaceNs(NoNuclides,EmptyX);
	SurfaceN = EmptySurfaceNs;
	
	//Setup CRN Arrays 	
	vector< vector<double> > EmptyN(NXNodes,EmptyZ);
	vector< vector< vector<double> > > EmptyNs(NoNuclides,EmptyN);
	N = EmptyNs;

	for (int j=0; j<NZNodes; ++j) Z[j] = (((ZMax-ZMin)/2.)-j*((ZMax-ZMin)/(NZNodes-1)));	
	 
	//Initialise Geomagnetic Scaling as constant
	GeoMagScalingFactor = 1;

	//Set Sea level to zero
	SeaLevel = 0;
}

void RockyCoastCRN::Initialise(double retreatrate, double beachwidth, int beachtype, double bermheight, double platformgradient, double cliffheight, double junctionelevation, double tidalamplitude, double slr, int steppedplatform, double stepsize)
{
	/* initialise a platform object for single retreat rate scenarios
	retreatrate is the rate of cliff retreat (m/yr)
	retreat type is the style of retreat, dfaults to 0 for single retreat rate scenario
	beachwidth is the width of the beach (m), protecting the platform from CRN production
	JunctionElevation is the elevation of the platform/cliff junction
	platformgradient is the average slope of the platform surface (m/m) 
	cliffheight is the height of the adjacent cliff (m)
	tidalamplitude is the average tidal amplitude for diurnal tides. */
	int retreattype = 0;
  double changetime = 0;
	Initialise(retreatrate, retreatrate, retreattype, changetime, beachwidth, beachtype, bermheight, platformgradient, cliffheight, junctionelevation, slr, tidalamplitude, steppedplatform, stepsize);
}
	
void RockyCoastCRN::Initialise(double retreatrate1, double retreatrate2, int retreattype, double changetime, double beachwidth, int beachtype, double bermheight, double platformgradient, double cliffheight, double junctionelevation, double SLR, double tidalamplitude, int steppedplatform, double stepsize)
{
  /* initialise a platform object for two retreat rate scenarios
  retreatrate1 runs from 7.5ka until changetime, after which retreatrate2 continues to present
	retreatrate1 is the first rate of cliff retreat (m/yr)
	retreatrate2 is the 2nd rate of cliff retreat (m/yr)
	changetime is the time to switch from retreatrate1 to retreatrate2 (years BP)
	beachwidth is the width of the beach (m), protecting the platform from CRN production
	JunctionElevation is the elevation of the platform/cliff junction
	platformgradient is the average slope of the platform surface (m/m) 
	cliffheight is the height of the adjacent cliff (m)
	tidalamplitude is the average tidal amplitude for diurnal tides. */
	
	//Geometric parameters
	NXNodes = 201;											//Number of nodes in cross shore
	NZNodes = 201;											//Number of nodes in profile
	PlatformWidth = 1000.;						//width of model domain (m)
	PlatformDepth = 20.;							//Depth to which CRN will be tracked (needs to be large enough that sea level rise is ok)
	NDV = -9999;											//Place holder for no data
	
	//Assign Parameters
	RetreatRate1 = retreatrate1;
	RetreatRate2 = retreatrate2;
	RetreatType = retreattype;
	ChangeTime = changetime;
	PlatformGradient = platformgradient;
	CliffHeight = cliffheight;
	BeachWidth = beachwidth;
	MeanBeachWidth = BeachWidth;
	InitialBeachWidth = BeachWidth;
	BeachType = beachtype;
	BermHeight = bermheight;
	JunctionElevation = junctionelevation;
	SLRRate = SLR;
	TidalAmplitude = tidalamplitude;
	
	SteppedPlatform = steppedplatform;
	StepSize = stepsize;

	TidalPeriod = 12.42;
		
	//constant for Bruun profile beaches
	A = 0.125;
		
	//initialise tides, geomag and RSL data
	InitialiseTides(TidalAmplitude,TidalPeriod);
	//InitialiseGeomagData();
	//InitialiseRSLData();
	//Initialise Geomagnetic Scaling as constant
	GeoMagScalingFactor = 1;
	//Initialise Platform
	InitialisePlanarPlatformMorphology();
}

void RockyCoastCRN::Initialise(RoBoCoP RoBoCoPCoast, vector<int> WhichNuclides)
{
	/* initialise a RockyCoastCRN object for use with a RoBoCoP object */
	printf("\nRockyCoastCRN.Initialise: Initialised a RockyCoastCRN object for use with a RoBoCoP object\n");

	//set tidal amplitude based on RoBoCoP
	TidalAmplitude = RoBoCoPCoast.TidalAmplitude;

	//get max extent of shoreface from RoBoCoP
	XMin = RoBoCoPCoast.X[0];
	XMax = RoBoCoPCoast.X[RoBoCoPCoast.NoNodes-1];
	double ZMax = RoBoCoPCoast.Z[0];
	double ZMin = RoBoCoPCoast.Z[RoBoCoPCoast.NoNodes-1];
	//Cliff is at the start of the vector
	CliffHeight = ZMax;
	CliffPositionX = XMin;
	CliffPositionInd = 0;
	
	//Setup CRN domain based on RoBoCoP
	NXNodes = RoBoCoPCoast.NoNodes;
	NZNodes = RoBoCoPCoast.NoNodes;
	NDV = -9999;
	dX = (XMax-XMin)/(NXNodes-1);
	dt = RoBoCoPCoast.dt;
	
	//Initialise the nuclides
	InitialiseNuclides(WhichNuclides);
	
	//Setup Surface Arrays
	vector<double> EmptyZ(NZNodes,0);
	vector<double> EmptyX(NXNodes,0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	X = EmptyX;
	Z = EmptyZ;
	PlatformElevation = EmptyXNDV;
	PlatformElevationOld = EmptyXNDV;
	SurfaceElevation = EmptyXNDV;
	vector< vector <double> > EmptySurfaceNs(NoNuclides,EmptyX);
	SurfaceN = EmptySurfaceNs;
	
	//Setup CRN Arrays 	
	vector< vector<double> > EmptyN(NXNodes,EmptyZ);
	vector< vector< vector<double> > > EmptyNs(NoNuclides,EmptyN);
	N = EmptyNs;
	
	for (int j=0; j<NZNodes; ++j) Z[j] = (((ZMax-ZMin)/2.)-j*((ZMax-ZMin)/(NZNodes-1)));	

	//Get surface morphology from RoBoCoP
	int Ind = 0;
	for (int i=0; i<NXNodes; ++i)
	{
		X[i] = dX*i;
		while ((RoBoCoPCoast.X[Ind] < X[i]) && (Ind < RoBoCoPCoast.NoNodes)) ++Ind;
		if (RoBoCoPCoast.X[Ind] == X[i]) SurfaceElevation[i] = RoBoCoPCoast.Z[Ind];
		else SurfaceElevation[i] = RoBoCoPCoast.Z[Ind-1] + (RoBoCoPCoast.Z[Ind]-RoBoCoPCoast.Z[Ind-1])*((X[i]-RoBoCoPCoast.X[Ind-1])/(RoBoCoPCoast.X[Ind]-RoBoCoPCoast.X[Ind-1]));
	}

	// copy to PlatformElevation
	PlatformElevation = SurfaceElevation;

	//Initialise Geomagnetic Scaling as constant
	GeoMagScalingFactor = 1;
}

void RockyCoastCRN::Initialise(Hiro HiroCoast, vector<int> WhichNuclides)
{
	/* initialise a RockyCoastCRN object for use with a Hiro object */
	printf("\nRockyCoastCRN.Initialise: Initialised a RockyCoastCRN object for use with a Hiro Model object\n");

	//set tidal range based on HiroCoast
	TidalRange = HiroCoast.TidalRange;
	
	//Setup CRN domain based on HiroCoast
	NXNodes = HiroCoast.NXNodes;
	NZNodes = HiroCoast.NZNodes;
	NDV = -9999;
	dX = HiroCoast.dX;
	dZ = HiroCoast.dZ;
	dt = HiroCoast.dt;
	
	//Initialise the nuclides
	InitialiseNuclides(WhichNuclides);
	
	//Setup Surface Arrays
	vector<double> EmptyZ(NZNodes,0);
	vector<double> EmptyX(NXNodes,0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	X = EmptyX;
	Z = HiroCoast.Z;
	PlatformElevation = EmptyXNDV;
	PlatformElevationOld = EmptyXNDV;
	SurfaceElevation = EmptyXNDV;
	vector< vector <double> > EmptySurfaceNs(NoNuclides,EmptyX);
	SurfaceN = EmptySurfaceNs;
	
	//Setup CRN Arrays 	
	vector< vector<double> > EmptyN(NXNodes,EmptyZ);
	vector< vector< vector<double> > > EmptyNs(NoNuclides,EmptyN);
	N = EmptyNs;
	
	//Setup tides
	TidalRange = HiroCoast.TidalRange;
	
	//Initialise Geomagnetic Scaling as constant
	GeoMagScalingFactor = 1;
	
	//Get morphology using update function
	UpdateMorphology(HiroCoast);
}

void RockyCoastCRN::InitialiseNuclides(vector<int> WhichNuclides)
{
	// Spallation (a/g/yr) calibrated 10Be production rate (Lifton et al. 2014).
	Po_10Be_Spal = 4.007;

	// Spallation (a/g/yr) calibrated 14C production rate (add ref).
	Po_14C_Spal = 12.29;

	// Spallation (a/g/yr) calibrated 26Al production rate (add ref).
	Po_26Al_Spal = 4.007;

	// Spallation (a/g/yr) calibrated 36Cl production rate (add_ref).
	Po_36Cl_Spal = 4.007;

	// Total muogenic production rate (a/g/yr) following Braucher et al. (2013).
	Po_10Be_Muon_Fast = 0.028;
	Po_10Be_Muon_Slow = 0;

	// Total muogenic production rate (a/g/yr)
	Po_14C_Muon_Fast = 3.34;
	Po_14C_Muon_Slow = 0.44;

	// Total muogenic production rate (a/g/yr)
	Po_26Al_Muon_Fast = 3.34;
	Po_26Al_Muon_Slow = 0.44;

	// Total muogenic production rate (a/g/yr)
	Po_36Cl_Muon_Fast = 3.34;
	Po_36Cl_Muon_Slow = 0.44;

	// Attenuation Lengths
	//Spallogenic attenuation length (kg/m^2).
	Lamb_Spal = 1600.0;

	// Muogenic attenuation length (kg/m^2) following Braucher et al. (2013).
	Lamb_Muon = 42000.0;	

	// density of rock and water respectively
	rho_r = 2600.; 
	rho_w = 1035.;

	//Decay length scales
	z_rs = Lamb_Spal/rho_r;		//Decay length scale chalk spallation
	z_ws = Lamb_Spal/rho_w;		//Decay length scale sea water spallation
	z_rm = Lamb_Muon/rho_r;		//Decay length scale chalk muons
	z_wm = Lamb_Muon/rho_w;		//Decay length scale sea water muons

	// Half life of 10Be (Korschineck et al. 2010).
	Lambda_10Be = 4.99*pow(10.,-7);

	// Half life of 14C (ref).
	Lambda_14C = 4.99*pow(10.,-7);

	// Half life of 26Al (ref).
	Lambda_26Al = 4.99*pow(10.,-7);

	// Half life of 36Cl (ref).
	Lambda_36Cl = 4.99*pow(10.,-7);

	//Declare number of nuclides
	NoNuclides = Nuclides.size();
	
	//Setup empty vectors to store the required global production variables
	Nuclides = WhichNuclides;
	NoNuclides = Nuclides.size();
	
	vector<double> EmptyVec(NoNuclides);
	Po_Spal = EmptyVec;
	Po_Muon_Fast = EmptyVec;
	Po_Muon_Slow = EmptyVec;
	Lambda  = EmptyVec;
	
	for (int n=0; n<NoNuclides; ++n)
	{
		if (Nuclides[n] == 10)
		{
			Po_Spal[n] = Po_10Be_Spal;
			Po_Muon_Fast[n] = Po_10Be_Muon_Fast;
			Po_Muon_Slow[n] = Po_10Be_Muon_Slow;
			Lambda[n] = Lambda_10Be;
		}
		else if (Nuclides[n] == 14)
		{
			Po_Spal[n] = Po_14C_Spal;
			Po_Muon_Fast[n] = Po_14C_Muon_Fast;
			Po_Muon_Slow[n] = Po_14C_Muon_Slow;
			Lambda[n] = Lambda_14C;
		}
		else if (Nuclides[n] == 26)
		{
			Po_Spal[n] = Po_26Al_Spal;
			Po_Muon_Fast[n] = Po_26Al_Muon_Fast;
			Po_Muon_Slow[n] = Po_26Al_Muon_Slow;
			Lambda[n] = Lambda_26Al;
		}
		else if (Nuclides[n] == 36)
		{
			Po_Spal[n] = Po_36Cl_Spal;
			Po_Muon_Fast[n] = Po_36Cl_Muon_Fast;
			Po_Muon_Slow[n] = Po_36Cl_Muon_Slow;
			Lambda[n] = Lambda_36Cl;
		}
		else cout 	<< "WARNING: Unrecognised isotope, atomic"
						<<	"mass needs to be from a recognised cosmogenic"
						<< "isotope such as 10Be, 14C, 26Al, 36Cl";		
	}
}

void RockyCoastCRN::InitialisePlanarPlatformMorphology()
{
 	//Setup Surface Arrays
	vector<double> EmptyZ(NZNodes,0);
	vector<double> EmptyX(NXNodes,0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	X = EmptyX;
	Z = EmptyZ;
	PlatformElevation = EmptyXNDV;
	PlatformElevationOld = EmptyXNDV;
	SurfaceElevation = EmptyXNDV;
	vector< vector <double> > EmptySurfaceNs(NoNuclides,EmptyX);
	SurfaceN = EmptySurfaceNs;
	
	//Setup CRN Arrays 	
	vector< vector<double> > EmptyN(NXNodes,EmptyZ);
	vector< vector< vector<double> > > EmptyNs(NoNuclides,EmptyN);
	N = EmptyNs;
	
	// Create initial morphology
	for (int j=0; j<NZNodes; ++j) Z[j] = ((PlatformDepth/2.)-j*(PlatformDepth/(NZNodes-1)));
	for (int i=0; i<NXNodes; ++i) X[i] = (i*(PlatformWidth/(NXNodes-1)));
}

void RockyCoastCRN::InitialiseTides(double A, double T)
{
	/// TIDES 
	TidalAmplitude = A;
	TidalPeriod = T;
	for (double TT = 0; TT <= TidalPeriod; TT += 0.2) TideLevels.push_back(-TidalAmplitude*sin((2.*M_PI*TT)/(TidalPeriod)));
	NTidalValues = (double)TideLevels.size();
	WaterDepths.resize(NTidalValues);
	WaterLevels.resize(NTidalValues);
}

void RockyCoastCRN::InitialiseGeomagData()
{	
	/// GEOMAG VARIATION
	double indata;
	char Dummy[32];
	
	//read in geomag data from Lifton et al. (2014) model
	string GeoMagFilename = "Lifton_GeoMagModel_Sussex.data";
	ifstream GeoMagIn(GeoMagFilename.c_str());
	if (!GeoMagIn)
	{
		printf("RockyCoastCRN::%s: line %d GeoMag data file \"%s\" doesn't exist\n\n", __func__, __LINE__, GeoMagFilename.c_str());
		printf("Setting Geomagnetic scaling factor to 1");
		GeoMagScalingFactor = 1;	
	}
	else
	{
		GeoMagIn >> Dummy;
		GeoMagIn >> Dummy;
		while (!GeoMagIn.eof())
		{
		  GeoMagIn >> indata;
		  GeoMagTime.push_back(indata);
		  GeoMagIn >> indata;
		  GeoMagScalingFactors.push_back(indata);
		}
		GeoMagIn.close();
	}
}

void RockyCoastCRN::InitialiseRSLData()
{
	//AND RELATIVE SEALEVEL
	double indata;
	char Dummy[32];
	
	//read in RSL data from Bradley et al. (2011) model
	string GIAFilename = "Bradley_GIAModel_Sussex.data";
	ifstream GIAIn(GIAFilename.c_str());
	if (!GIAIn)
	{
	  printf("RockyCoastCRN::%s: line %d Relative Sea Level data file \"%s\" doesn't exist\n\n", __func__, __LINE__, GIAFilename.c_str());
	  printf("Setting Realtive Sea Level to zero");
	}
	else
	{
		GIAIn >> Dummy;
		GIAIn >> Dummy;
		while (!GIAIn.eof())
		{
		  GIAIn >> indata;
		  RSLTime.push_back(indata);
		  GIAIn >> indata;
		  RSLRate.push_back(indata);
		}
		GIAIn.close();
	}
}
	
void RockyCoastCRN::UpdateParameters( double RetreatRate1_Test, double RetreatRate2_Test, 
                                    double ChangeTime_Test, double BeachWidth_Test, 
                                    double JunctionElevation_Test)
{
  /*  Function to update model parameters for a new model run. 
  This is designed for use with the MCMC_RockyCoast object in order to save a 
  bit of time on reinitialising each time in the Markov chain 
  RetreatRate1_Test is the new RetreatRate1
  RetreatRate2_Test is the new RetreatRate2
  ChangeTime_Test is the new ChangeTime
  BeachWidth_Test is the new BeachWidth
  JunctionElevation_Test is the new JunctionElevation    */
  
  RetreatRate1 = RetreatRate1_Test;
  RetreatRate2 = RetreatRate2_Test;
  ChangeTime = ChangeTime_Test;
  BeachWidth = BeachWidth_Test;
  MeanBeachWidth = BeachWidth;
  InitialBeachWidth = BeachWidth;
  JunctionElevation = JunctionElevation_Test;

}

void RockyCoastCRN::RunModel(string outfilename, int WriteResultsFlag)
{
	/*  Main model loop. Iterates through time updating the CRN concentrations on the 
		platform and at depth. Cliff steps back following cliff retreat rates, and 
		platform downwears to maintain steady-state cliff profile. This currently only
		works for a constant and gradual rate of cliff retreat and platform downwear. 
		Future work can explore the influence of block removal on the platform, and 
		episodic cliff retreat scenarios.
		
		Retreat Type can be 
		  - 0 for single/constant retreat rates
		  - 1 for step change retreat rates at time ChangeTime
		  - 2 for linear/gradual change in retreat rates through time
		  
		WriteResultsFlag specifies whether to write results to file. Default is 1 (on).
		This should be turned off (set to 0) for running MCMC_RockyCoast analysis.
	*/

	//Filenames
	OutFileName = outfilename;

	//Setup Surface Arrays
	vector<double> EmptyZ(NZNodes,0);
	vector<double> EmptyX(NXNodes,0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	X = EmptyX;
	Z = EmptyZ;
	PlatformElevation = EmptyXNDV;
	PlatformElevationOld = EmptyXNDV;
	SurfaceElevation = EmptyXNDV;
	vector< vector <double> > EmptySurfaceNs(NoNuclides,EmptyX);
	SurfaceN = EmptySurfaceNs;
	
	//Setup CRN Arrays 	
	vector< vector<double> > EmptyN(NXNodes,EmptyZ);
	vector< vector< vector<double> > > EmptyNs(NoNuclides,EmptyN);
	N = EmptyNs;
	
	//set Sea level parameters
	//SLR = 0.0002;     //Rate of relative sea level rise (m/y)  
	SeaLevel = 0;			//Sea Level Tracker	
	MeanBeachWidth = BeachWidth;
	InitialBeachWidth = BeachWidth;
	
	//Time control
	Time = 0;			//time in years
	dt = 5;		    //time step
	XMax = 1000.; 	    //Distance from cliff
	
	//Write output control
	double WriteTime;
	double WriteInterval = 1000;
	
	//Work out start time from assumed cliff retreat rates.
	if (RetreatType == 0) MaxTime = XMax/RetreatRate1;
	else if (RetreatType == 1)
	{
	  double ChangeX = ChangeTime*RetreatRate2;
	  MaxTime = ChangeTime + (XMax-ChangeX)/RetreatRate1;
	}
	else if (RetreatType == 2) MaxTime = 7000;
	else 
	{
	  cout << "Retreat Type Unknown" << endl;
	  exit(EXIT_SUCCESS);
	}

  //limit start time to 10ka
	//MaxTime = 10000;
	
	cout << "Start time is " << MaxTime << " ka" << endl;
	Time = MaxTime;
	WriteTime = MaxTime;
	
	//setup CliffPositionX for tracking cliff position in X
	if (RetreatType == 0) CliffPositionX = MaxTime*RetreatRate1;
	else if (RetreatType == 1) CliffPositionX = ChangeTime*RetreatRate2 + (Time-ChangeTime)*RetreatRate1; 		      
	else if (RetreatType == 2) CliffPositionX = MaxTime*(RetreatRate1+RetreatRate2)/2.;
	else 
	{
	  cout << "Retreat Type Unknown" << endl;
	  exit(EXIT_SUCCESS);
	}
	
	//Update XMax to match CliffPositionX (i.e. the start point)
	XMax = CliffPositionX;
	
	//Get surface elevations
	CliffPositionInd = 0;
	ZTrackInd = 0;
	
	double TempBeachThickness;
	
	for (int i=0; i<NXNodes; ++i)
	{
	  if (X[i] > CliffPositionX) 
	  {
	    //update surface elevation
	    PlatformElevation[i] = JunctionElevation-(X[i]-CliffPositionX)*PlatformGradient;
	    if (SteppedPlatform == 1) PlatformElevation[i] = round(PlatformElevation[i]/StepSize)*StepSize;
	    
	    //update beach thickness
	    if ((X[i]-CliffPositionX) < BeachWidth) TempBeachThickness = JunctionElevation+BermHeight-PlatformElevation[i];
	    else TempBeachThickness = (JunctionElevation+BermHeight-A*pow((X[i]-CliffPositionX-BeachWidth),2./3.)) - PlatformElevation[i];
	    if (TempBeachThickness > 0) BeachThickness[i] = TempBeachThickness;
	    else BeachThickness[i] = 0;	    
	  }
	  else CliffPositionInd = i;
	  
	  SurfaceElevation[i] = PlatformElevation[i]+BeachThickness[i];
	}
  PlatformElevationOld = PlatformElevation;
  
  
  //Set Z tracker to keep track of elevation index at the platform/cliff junction
  for (int j=0; j<NZNodes; ++j) if (Z[j] > SeaLevel+JunctionElevation+BermHeight) ZTrackInd = j;
	
	//////////////////////////////////////////////////////////
	//MAIN MODEL LOOP
	//////////////////////////////////////////////////////////

	while (CliffPositionX > 0)
	{
		//Get retreat rate depending on style of retreat
		GetRetreatRate();
		
		//Update beachwidth
		if (BeachType == 1) GetSinWaveBeachWidth(Time);
		else if (BeachType == 2) GetThinningBeachWidth(Time);
		
		//Get geomag scaling factor
		GeoMagScalingFactor = GetGeoMagScalingFactor(Time);
		
		//Get sea level rise from local record
		//SLR = GetSeaLevelRise(Time);
		
    //update CRN concentrations
    UpdateCRNs();
    
    //update morphology
    UpdateEquillibriumMorphology();
    
    //write output?
    if ((WriteResultsFlag != 0) && (Time <= WriteTime))
    {
      WriteProfile(OutFileName, Time);
      WriteCRNProfile(OutFileName, Time);
      WriteTime -= WriteInterval;
    }
    
		//update cliff position and time
		CliffPositionX -= RetreatRate*dt;
		Time -= dt;
		SeaLevel += SLRRate*dt;
		if (CliffPositionX < X[CliffPositionInd]) CliffPositionInd -= 1;
		//if (SeaLevel < Z[ZTrackInd]) ZTrackInd += 1;
    if (CliffPositionInd < 0) CliffPositionInd = 0;
    
	}
	
	WriteProfile(OutFileName, Time);
	WriteCRNProfile(OutFileName, Time);
	
  CliffPositionInd = 0;
  
}

void RockyCoastCRN::GetRetreatRate()
{
  /*
  Function to update get the retreat rate based on the retreat scenario
  where 0 is constant retreat rate, 1 is step change in retreat rate and
  2 is linear change in retreat rate
  
  Martin Hurst
  February 2016
  */
  
  //temp parameters
  double dTime, Factor,dRR;
  if (RetreatType == 0) RetreatRate = RetreatRate1;
	else if (RetreatType == 1)
	{
	  if (Time > ChangeTime) RetreatRate = RetreatRate1;
	  else RetreatRate = RetreatRate2;
	}
	else if (RetreatType == 2)
	{
	  dTime = MaxTime-Time;
	  Factor = dTime/MaxTime;
	  dRR = RetreatRate1-RetreatRate2;
	  RetreatRate = RetreatRate1-dRR*Factor;
	}
	else 
	{
	  printf("RockyCoastCRN::%s: line %d Unknown retreat type!\n\n", __func__, __LINE__);
	  exit(EXIT_SUCCESS);
	}
}

void RockyCoastCRN::UpdateCRNs()
{
	/*
	Function to update the concentrations of cosmogenic radionuclides at and below 
	the platform surface. Currently for 10Be only.

	Martin Hurst
	February 2016
	*/

	//Temp parameters
	vector<double> EmptyNVec(NXNodes,0);
	vector <vector <double> > P_Spal(NoNuclides,EmptyNVec);
	P_Muon_Fast = P_Spal;
	P_Muon_Slow = P_Spal;
	
	// variation in production on the platform at depth as a function of a tidal period
	// First get a Tidal Sequence of Water Levels, clipping for negative depths
	for (int i=0, NN=WaterLevels.size(); i<NN;++i) WaterLevels[i] = SeaLevel+TideLevels[i];

	//LOOP Through the concentrations arrays
	for (int j=0; j<NXNodes; ++j)
	{
		//Get topographic shielding factor
		TopoShieldingFactor = GetTopographicShieldingFactor(X[j]-CliffPositionX, CliffHeight);
		
		//FOR EACH NUCLIDE OF INTEREST SET SURFACE PRODUCTION RATES
		for (int n=0; n<NoNuclides; ++n)
		{
			P_Spal[n][j] = Po_Spal[n];
			P_Muon_Fast[n][j] = Po_Muon_Fast[n];
			P_Muon_Slow[n][j] = Po_Muon_Slow[n];
		}
			
		//Sort out the water shielding
		if (SurfaceElevation[j] < SeaLevel+0.5*TidalRange)
		{
			//if under water calculate production modified for water depth
			for (int n=0; n<NoNuclides; ++n)
			{
				P_Spal[n][j] = 0;
				P_Muon_Fast[n][j] = 0;
				P_Muon_Slow[n][j] = 0;
			}
			
			//Loop over the tidal cycle and total the production
			for (int a=0, NN=WaterLevels.size(); a<NN; ++a)
			{
				if (WaterLevels[a] >= SurfaceElevation[j]) WaterDepths[a] = WaterLevels[a]-SurfaceElevation[j];
				else WaterDepths[a] = 0;

				//Calculate Production for this profile
				//FOR EACH NUCLIDE OF INTEREST
				for (int n=0; n<NoNuclides; ++n)
				{
					P_Spal[n][j] += GeoMagScalingFactor*TopoShieldingFactor*Po_Spal[n]*exp(-WaterDepths[a]/z_ws);
					P_Muon_Fast[n][j] += TopoShieldingFactor*Po_Muon_Fast[n]*exp(-WaterDepths[a]/z_wm);
					P_Muon_Slow[n][j] += TopoShieldingFactor*Po_Muon_Slow[n]*exp(-WaterDepths[a]/z_wm);
				}
			}
			//FOR EACH NUCLIDE OF INTEREST
			for (int n=0; n<NoNuclides; ++n)
			{
				//find mean production rate at surface
				//normalise by the tidal duration
				P_Spal[n][j] /= NTidalValues;
				P_Muon_Fast[n][j] /= NTidalValues;
				P_Muon_Slow[n][j] /= NTidalValues;
			}
		}
		
		bool Top = false;
		
		for (int i=0; i<NZNodes; ++i)
		{
			if (Z[i] > PlatformElevation[j])
			{
				for (int n=0; n<NoNuclides; ++n)
				{
					N[n][j][i] = 0;
				}
			}
			else if (Z[i] > PlatformElevation[j]-20.)
			{
				for (int n=0; n<NoNuclides; ++n)
				{
					//update concentrations at depth
					//This is kept as SurfaceElevation not Platform Elevation for now
					//NB This assumes that material density of the beach is the same as the bedrock!
					N[n][j][i] += dt*P_Spal[n][j]*exp((Z[i]-SurfaceElevation[j])/z_rs);	//spallation
					N[n][j][i] += dt*P_Muon_Fast[n][j]*exp((Z[i]-SurfaceElevation[j])/z_rm);	//muons
					N[n][j][i] += dt*P_Muon_Slow[n][j]*exp((Z[i]-SurfaceElevation[j])/z_rm);	//muons
					//remove atoms due to radioactive decay
					N[n][j][i] -= dt*Lambda[n];
				}
				if (Top == false)
				{
					for (int n=0; n<NoNuclides; ++n)
					{
						//update concentration at the platform surface, accounting for beach cover
						SurfaceN[n][j] = N[n][j][i];
						Top = 1;
					}
				}
			}
		}
	}
}

void RockyCoastCRN::UpdateEquillibriumMorphology()
{
  /*
  Function to update the morphology of the platform. Moves the cliff following the 
  prescribed retreat rate and updates the platform surface.
  
  Martin Hurst
  February 2016
  */
  
  //temp parameters
  double TempPlatformElevChange, TempBeachThickness;
	
  //Update surface elevations
  PlatformElevationOld = PlatformElevation;
  
	for (int i=CliffPositionInd; i<NXNodes; ++i)
	{
		if (X[i] >= CliffPositionX)
		{
			if (PlatformElevation[i] == NDV)
			{
			  PlatformElevation[i] = SeaLevel+JunctionElevation-(X[i]-CliffPositionX)*PlatformGradient;
			  if (SteppedPlatform == 1) PlatformElevation[i] = round(PlatformElevation[i]/StepSize)*StepSize;
			}
			else
			{
				if (SteppedPlatform == 0) 
				{
				  TempPlatformElevChange = dt*(RetreatRate*PlatformGradient-SLRRate);
				  if (TempPlatformElevChange < 0) TempPlatformElevChange = 0;
				  PlatformElevation[i] -= TempPlatformElevChange;
				}
				else if (SteppedPlatform == 1) 
				{
				  PlatformElevation[i] = SeaLevel+JunctionElevation-(X[i]-CliffPositionX)*PlatformGradient;
				  PlatformElevation[i] = round(PlatformElevation[i]/StepSize)*StepSize;
				}
			}
			
	    //update beach thickness
	    if ((X[i]-CliffPositionX) < BeachWidth) TempBeachThickness = JunctionElevation+BermHeight-PlatformElevation[i];
	    else TempBeachThickness = BermHeight-A*pow((X[i]-CliffPositionX-BeachWidth),2./3.);

	    if (TempBeachThickness > 0) BeachThickness[i] = TempBeachThickness;
	    else BeachThickness[i] = 0;	    
	    
	    SurfaceElevation[i] = PlatformElevation[i]+BeachThickness[i];
		}
	}
}

//void RockyCoastCRN::UpdateMorphology(RoBoCoP RoBoCoPCoast)
//{
//  // Find cliff position
//  XMin = RoBoCoPCoast.X[0];
//  XMax = RoBoCoPCoast.X[RoBoCoPCoast.NoNodes-1];
//  
//  CliffPositionX = XMin;
//  CliffPositionInd = 0;
//  
//  PlatformElevationOld = PlatformElevation;
//  
//  // Add nodes to front of RockyCoastCRN as required
//  vector<double> EmptyZ(NZNodes,0.0);
//  while (XMin < X[0]) 
//  {
//    X.insert(X.begin(),X[0]-dX);
//    SurfaceElevation.insert(SurfaceElevation.begin(), NDV);
//    PlatformElevation.insert(PlatformElevation.begin(), NDV);
//    PlatformElevationOld.insert(PlatformElevationOld.begin(), NDV);
//    SurfaceN.insert(SurfaceN.begin(),0);
//    N.insert(N.begin(),EmptyZ);
//    ++NXNodes;
//  }
//  
//  // Get surface morphology from RoBoCoP
//  int Ind = 0;
//  SurfaceElevation[0] = CliffHeight;
//	for (int i=1; i<NXNodes; ++i)
//  {
//    while ((RoBoCoPCoast.X[Ind] < X[i]) && (Ind < RoBoCoPCoast.NoNodes)) ++Ind;
//    if (RoBoCoPCoast.X[Ind] == X[i]) SurfaceElevation[i] = RoBoCoPCoast.Z[Ind];
//    else SurfaceElevation[i] = RoBoCoPCoast.Z[Ind-1] + (RoBoCoPCoast.Z[Ind]-RoBoCoPCoast.Z[Ind-1])*((X[i]-RoBoCoPCoast.X[Ind-1])/(RoBoCoPCoast.X[Ind]-RoBoCoPCoast.X[Ind-1]));
//  }	
//  // Copy Z to SurfaceElevation
//  // Do we need both?
//  PlatformElevation = SurfaceElevation;
//  
//  // Copy sea level
//  SeaLevel = RoBoCoPCoast.SeaLevel;
//}
void RockyCoastCRN::UpdateMorphology(Hiro HiroCoast)
{
	//Get number of nodes in coastal morphology
	X = HiroCoast.X;
	Zx = HiroCoast.Zx;
	SurfaceInd = HiroCoast.ZInd;
	
	int TempXSize = X.size();
	while (TempXSize > NXNodes)
	{
		for (int n=0; n<NoNuclides; ++n)
		{
			N[n].push_back(N[n][NXNodes-1]);
			SurfaceN[n].push_back(SurfaceN[n][NXNodes-1]);
		}
		PlatformElevation.push_back(PlatformElevation[NXNodes-1]);
		--TempXSize;
	}
	NXNodes = X.size();
	 
	PlatformElevationOld = PlatformElevation;
	PlatformElevation = Zx;
	SurfaceElevation = Zx;
	
	for (int j=0; j<NXNodes; ++j)
	{
		for (int n=0; n<NoNuclides; ++n)
		{
			SurfaceN[n][j] = N[n][j][SurfaceInd[j]];
		}
	}
	
	//Cliff is on the right, find it
	//This will need updating once we have sea level rise
	CliffHeight = Z[NZNodes-1];
	CliffPositionX = 0;
	for (int j=NXNodes-1;j>0; --j)
	{
		if (Zx[j] < Zx[NXNodes-1])
		{
			CliffPositionX = X[j];
			break;
		}
	}
	
	//update Sea levl
	SeaLevel = HiroCoast.SeaLevel;
	TidalRange = HiroCoast.TidalRange;
}

void RockyCoastCRN::UpdateMorphology(vector<double> XCoast, vector<double> ZCoast)
{
	//Get number of nodes in coastal morphology
	int NoCoastNodes = XCoast.size();
	
	// Find cliff position
	XMin = XCoast[0];
	XMax = XCoast[NoCoastNodes-1];
  
	CliffPositionX = XCoast[0];
	CliffPositionInd = 0;
  
	PlatformElevationOld = PlatformElevation;
  
	// may come back and add dynamic growth of CRN model later
//	// Add nodes to front of RockyCoastCRN as required
//	vector<double> EmptyZ(NZNodes,0.0);
//	while (XMin < X[0]) 
//	{
//		X.insert(X.begin(),X[0]-dX);
//		SurfaceElevation.insert(SurfaceElevation.begin(), NDV);
//		PlatformElevation.insert(PlatformElevation.begin(), NDV);
//		PlatformElevationOld.insert(PlatformElevationOld.begin(), NDV);
//		SurfaceN.insert(SurfaceN.begin(),0);
//		N.insert(N.begin(),EmptyZ);
//		++NXNodes;
//	}
  
	// Get surface morphology from X and Z vectors
	int Ind = 0;
	SurfaceElevation[0] = CliffHeight;
	for (int i=1; i<NXNodes; ++i)
	{
		 while ((XCoast[Ind] < X[i]) && (Ind < NoCoastNodes)) ++Ind;
		 if (XCoast[Ind] == XCoast[i]) SurfaceElevation[i] = ZCoast[Ind];
		 else SurfaceElevation[i] 	= ZCoast[Ind-1] 
		 									+ (ZCoast[Ind]-ZCoast[Ind-1])*((X[i]-XCoast[Ind-1])/(XCoast[Ind]-XCoast[Ind-1]));
	}
	
	// Copy Z to SurfaceElevation
	// Do we need both?
	PlatformElevation = SurfaceElevation;
}

double RockyCoastCRN::GetTopographicShieldingFactor(double X, double CliffHeight)
{
	/* 
	Function to calculate topographic shielding factor due to presence of a cliff.
	Cliff is assumed to be straight in profile, of fixed height, and vertical
	 
	Martin Hurst
	October 2014
	*/    

	if (X == CliffPositionX-dX) return 0.5;
	else if (CliffPositionX-X > 5*CliffHeight) return 1.0;
	else if (X<CliffPositionX)
	{
		//##### Cliff Shielding parameters #####
		double d_theta_phi = (M_PI/180.)*5.0;				//azimuth and angle stepping
		double FMax = 2.0*M_PI*Po_10Be_Spal/(3.3);		//Maximum Intensity

		double Viewshed;
		double F = 0;

		for (double Az=-90; Az<=90.; Az+= 5.) 
		{
		Viewshed = atan(CliffHeight*cos((M_PI/180.)*Az)/X);
		F+= d_theta_phi*Po_10Be_Spal/(3.3)*pow(sin(Viewshed),3.3);
		}

		return (FMax-F)/FMax;
	}
	else return 0;
}

double RockyCoastCRN::GetGeoMagScalingFactor(double Time)
{
	/*
	Function to get a Geomagnetic and solar variability related scaling factor
	based on data from Lifton et al 2014 scaling model

	Martin Hurst
	October 2014
	*/

	double GeoMagScalingFactor = 1;
  
	//interpolate to get value
	int TimeCondition = 0;
	int ind = 0;
	
	//if too far back in the past return 1
	if (GeoMagTime.size() == 0) return GeoMagScalingFactor;
	
	while (TimeCondition == 0)
	{
		if (Time > GeoMagTime[ind]) ++ind;
		else TimeCondition = 1;
	}
	double Factor = (Time-GeoMagTime[ind-1])/(GeoMagTime[ind]-GeoMagTime[ind-1]);
	GeoMagScalingFactor = GeoMagScalingFactors[ind-1] + Factor*(GeoMagScalingFactors[ind]-GeoMagScalingFactors[ind-1]);

	return GeoMagScalingFactor;
}

double RockyCoastCRN::GetSeaLevelRise(double Time)
{
	/*
	Function to get rate of Sea Level Rise from Bradley Model for Sussex
		
	Martin Hurst
	December 2014
	*/

	if (RSLTime.size() == 0) return SLRRate;
	else
	{
	  	//interpolate to get value
		double Factor, Rate;
		int TimeCondition = 0;
		int ind = 0;
		while (TimeCondition == 0)
		{
			if (Time > RSLTime[RSLTime.size()-1])
			{
				TimeCondition = 1;
				ind = RSLTime.size()-1;
			}
			else if (Time > RSLTime[ind]) ++ind;
			else TimeCondition = 1;
		}
		Factor = (Time-RSLTime[ind-1])/(RSLTime[ind]-RSLTime[ind-1]);
		Rate = RSLRate[ind-1]+Factor*(RSLRate[ind]-RSLRate[ind-1]);

		return -Rate;
	}
}

void RockyCoastCRN::GetSinWaveBeachWidth(double Time)
{
  // Function to return beach width where beach width varies as a sinosoidal 
  // of time.
  // Martin Hurst
  // Feb 11th 2016
  
  double BeachWidthAmp = 30;
  double BeachWidthWaveLength = 100;
  
  BeachWidth = MeanBeachWidth+BeachWidthAmp*sin((2.*M_PI*Time)/BeachWidthWaveLength);
}

void RockyCoastCRN::GetThinningBeachWidth(double Time)
{
  // Function to return beach width where beach width reduces with time
  // Martin Hurst
  // Feb 11th 2016
  
  double ThinTime = 1000;
  
  if (Time < ThinTime) BeachWidth = (Time/ThinTime)*InitialBeachWidth;
  else BeachWidth = InitialBeachWidth;
}

void RockyCoastCRN::WriteProfile(string OutputFileName, double Time)
{
  //test if output file already exists
  int FileExists = 0;
	ifstream oftest(OutputFileName.c_str());
	if (oftest) FileExists = 1;
	oftest.close();

	//open the output filestream and write headers
	ofstream WritePlatform;
	if (FileExists == 0)
	{
		WritePlatform.open(OutputFileName.c_str());
		if (WritePlatform.is_open()) WritePlatform << NXNodes << " " << PlatformWidth << endl;
	}
	WritePlatform.close();

	//open output filestream again to  coastline data
	WritePlatform.open(OutputFileName.c_str(), fstream::app|fstream::out);

	//Check if file exists if not open a new one and write headers
	if (WritePlatform.is_open())
	{
		//write PlatformElevation
		WritePlatform << Time;
		for (int i=0; i<NXNodes; ++i) WritePlatform << setprecision(5) << " " << X[i];
		WritePlatform << endl;
		
		//write SurfaceElevation
		WritePlatform << Time;
		for (int i=0; i<NXNodes; ++i) WritePlatform << setprecision(5) << " " << SurfaceElevation[i];
		WritePlatform << endl;
	}
}

void RockyCoastCRN::WriteCRNProfile(string OutputFileName, double Time)
{
	//test if output file already exists
	int FileExists = 0;
	ifstream oftest(OutputFileName.c_str());
	if (oftest) FileExists = 1;
	oftest.close();

	//open the output filestream and write headers
	ofstream WritePlatform;
	if (FileExists == 0)
	{
		WritePlatform.open(OutputFileName.c_str());
		if (WritePlatform.is_open()) 
		{
			//write which nuclides we used
			for (int n=0; n<NoNuclides; ++n) WritePlatform << Nuclides[n] << " ";
			WritePlatform << endl;
			
			//dX as header line #2
			WritePlatform << dX << endl;
		}
	}
	WritePlatform.close();

	//open output filestream again to  coastline data
	WritePlatform.open(OutputFileName.c_str(), fstream::app|fstream::out);

	//Check if file exists if not open a new one and write headers
	if (WritePlatform.is_open())
	{
		//write concentrations
		for (int n=0; n<NoNuclides; ++n)
		{
			WritePlatform << Time;
			for (int j=0; j<NXNodes; ++j) WritePlatform << setprecision(5) << " " << SurfaceN[n][j];
			WritePlatform << endl;
		}
	}
}

void  RockyCoastCRN::WriteNuclideArray(string OutputFileName, double Time, int Nuclide)
{
  /* Writes a matrix of N concentrations to file
		File format is 	
		
		Time	
			N[0][0]     |    N[0][1]    |   N[0][2]     =====>    N[0][NXNodes]
			N[1][0]     |    N[1][1]    |   N[1][2]     =====>    N[0][NXNodes]
			N[2][0]     |    N[2][1]    |   N[2][2]     =====>    N[0][NXNodes]
		      ||               ||             ||                      ||
		      \/               \/             \/                      \/
		N[NZNodes][0]  | N[NZNodes][1] | N[NZNodes][2] =====> N[NZNodes][NXNodes] */
      
  	//open the output filestream and write headers
	ofstream WriteFile;
	WriteFile.open(OutputFileName.c_str());
	WriteFile << Time << " " << dZ << " " << dX << endl;
	
	//Check if file exists if not open a new one and write headers
	if (WriteFile.is_open())
	{
		//write Resistance
		for (int i=0; i<NZNodes; ++i)
		{
			for (int j=0;j<NXNodes; ++j)
			{
				WriteFile << setprecision(6) << N[Nuclide][j][i] << " ";
			}
			WriteFile << endl;
		}
	}
	else
	{
		//report errors
		cout << "RockyCoastCRN.WriteNuclideArray: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
}
#endif
