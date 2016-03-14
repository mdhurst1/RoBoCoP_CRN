/*------------------------------------------------------------------------

	RockyCoastCRN.cpp
	
	Object to evolve the CRN concentration across a coastal profile as a 
	function of tides, topographic shielding, relative sea level rise, 
	cliff retreat and platform downwear. For comparison with measured
	10Be CRN measurements across a coastal platform.
	
	Martin D. Hurst, British Geological Survey, December 2014

	Copyright (C) 2015, Martin Hurst
	
	Developer can be contacted:
  mhurst@bgs.ac.uk
  
  Martin D. Hurst
  British Geological Survey,
  Environmental Science Centre,
  Nicker Hill,
  Keyworth,
  Nottingham,
  UK,
  NG12 5GG
  
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
#include "CRN_global_variables.hpp"

using namespace std;

void RockyCoastCRN::Initialise()
{
  /* initialise an empty platform object
	retreatrate is the rate of cliff retreat (m/yr)
	beachwidth is the width of the beach (m), protecting the platform from CRN production
	JunctionElevation is the elevation of the platform/cliff junction
	platformgradient is the average slope of the platform surface (m/m) 
	cliffheight is the height of the adjacent cliff (m)
	tidalamplitude is the average tidal amplitude for diurnal tides. */
	
  double NDV = -9999;
  RetreatRate1 = NDV;
	RetreatRate2 = NDV;
	ChangeTime = NDV;
	PlatformGradient = NDV;
	CliffHeight = NDV;
	BeachWidth = NDV;
	BermHeight = NDV;
	JunctionElevation = NDV;
	TidalAmplitude = NDV;
}

void RockyCoastCRN::Initialise(double retreatrate, double beachwidth, double bermheight, double platformgradient, double cliffheight, double junctionelevation, double tidalamplitude, int steppedplatform, double stepsize)
{
	/* initialise a platform object for single retreat rate scenarios
	retreatrate is the rate of cliff retreat (m/yr)
	beachwidth is the width of the beach (m), protecting the platform from CRN production
	JunctionElevation is the elevation of the platform/cliff junction
	platformgradient is the average slope of the platform surface (m/m) 
	cliffheight is the height of the adjacent cliff (m)
	tidalamplitude is the average tidal amplitude for diurnal tides. */
  double changetime = 0;
	Initialise(retreatrate, retreatrate, changetime, beachwidth, bermheight, platformgradient, cliffheight, junctionelevation, tidalamplitude,steppedplatform, stepsize);
}
	
void RockyCoastCRN::Initialise(double retreatrate1, double retreatrate2, int retreattype, double changetime, double beachwidth, double bermheight, double platformgradient, double cliffheight, double junctionelevation, double tidalamplitude, int steppedplatform, double stepsize)
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
	BermHeight = bermheight;
	JunctionElevation = junctionelevation;
	TidalAmplitude = tidalamplitude;
	SteppedPlatform = steppedplatform;
	StepSize = stepsize;
	
	//constant for Bruun profile beaches
	A = 0.125;
		
	//initialise tides, geomag and RSL data
	InitialiseTides();
	InitialiseGeomagData();
	InitialiseRSLData();
	
	//Initialise Platform
	InitialisePlanarPlatformMorphology();
}

void RockyCoastCRN::Initialise(RoBoCoP RoBoCoPCoast)
{
  /* initialise a RockyCoastCRN object for use with a RoBoCoP object */
  printf("\nRockyCoastCRN.Initialise: Initialised a RockyCoastCRN object for use with a RoBoCoP object\n");

  //set tidal amplitude based on RoBoCoP
	TidalAmplitude = RoBoCoPCoast.TidalAmplitude;
	
	//get max extent of shoreface from RoBoCoP
	double XMin = RoBoCoPCoast.X[0];
	double XMax = RoBoCoPCoast.X[RoBoCoPCoast.NoNodes-1];
	double ZMin = RoBoCoPCoast.Z[0];
	double ZMax = RoBoCoPCoast.Z[RoBoCoPCoast.NoNodes-1];
	CliffHeight = ZMax;
		
	//Setup CRN domain based on RoBoCoP
	NXNodes = RoBoCoPCoast.NoNodes;
	NZNodes = RoBoCoPCoast.NoNodes;
	NDV = -9999;
	double dX = (XMax-XMin)/(NXNodes);
		
	//Setup CRN Arrays
	vector<double> EmptyX(NXNodes,0.0);
	vector<double> EmptyZ(NZNodes,0.0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	vector< vector<double> > EmptyVV(NXNodes,EmptyZ);
	X = EmptyX;
	Z = EmptyZ;
	PlatformElevation = EmptyXNDV;
	PlatformElevationOld = EmptyXNDV;
	SurfaceElevation = EmptyXNDV;
	SurfaceN = EmptyX;
	N = EmptyVV;
	
	//Get surface morphology from RoBoCoP
  int Ind = 0;
	for (int i=0; i<RoBoCoPCoast.NoNodes; ++i)
  {
    X[i] = dX*i;
    while (RoBoCoPCoast.X[Ind] <= X[i]) ++Ind;
    if (RoBoCoPCoast.X[Ind] == X[i]) Z[i] = RoBoCoPCoast.Z[Ind];
    else Z[i] = RoBoCoPCoast.Z[Ind-1] + (RoBoCoPCoast.Z[Ind]-RoBoCoPCoast.Z[Ind-1])*((X[i]-RoBoCoPCoast.X[Ind-1])/(RoBoCoPCoast.X[Ind]-RoBoCoPCoast.X[Ind-1]));
  }	
}	
void RockyCoastCRN::InitialiseTides()
{
	/// TIDES 
	TidalPeriod=12.;
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
	  exit(EXIT_SUCCESS);
	}
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
	  exit(EXIT_SUCCESS);
	}
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

void RockyCoastCRN::InitialisePlanarPlatformMorphology()
{
  vector<double> EmptyX(NXNodes,0.0);
	vector<double> EmptyZ(NZNodes,0.0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	vector< vector<double> > EmptyVV(NXNodes,EmptyZ);
	X = EmptyX;
	Z = EmptyZ;
	PlatformElevation = EmptyXNDV;
	PlatformElevationOld = EmptyXNDV;
	SurfaceElevation = EmptyXNDV;
	SurfaceN = EmptyX;
	N = EmptyVV;
	for (int j=0; j<NZNodes; ++j) Z[j] = ((PlatformDepth/2.)-j*(PlatformDepth/(NZNodes-1)));
	for (int i=0; i<NXNodes; ++i) X[i] = (i*(PlatformWidth/(NXNodes-1)));
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
  
  //Reset Concentrations
  //Setup Vectors
	vector<double> EmptyX(NXNodes,0.0);
	vector<double> EmptyZ(NZNodes,0.0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	vector< vector<double> > EmptyVV(NXNodes,EmptyZ);
	
	PlatformElevation = EmptyXNDV;
	PlatformElevationOld = EmptyXNDV;
	SurfaceElevation = EmptyXNDV;
	BeachThickness = EmptyX;
	
	SurfaceN = EmptyX;
	N = EmptyVV;
	
	//set Sea level parameters
	SLR = 0;          //Rate of relative sea level rise (m/y)  
	SeaLevel = 0;			//Sea Level Tracker	
	
	//Time control
	Time = 0;			//time in years
	dt = 1;		    //time step
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

  //limit start time to 7ka
	//if (MaxTime > 7000)	Time = 7000;
	//else Time = MaxTime;
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
	    else TempBeachThickness = PlatformElevation[i] + (JunctionElevation+BermHeight-A*pow((X[i]-CliffPositionX-BeachWidth),2./3.));
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
		//GetThinningBeachWidth(Time);
		
		//Get geomag scaling factor
		GeoMagScalingFactor = GetGeoMagScalingFactor(Time);
		
		//Get sea level rise from local record
		//SLR = GetSeaLevelRise(Time);
		
    //update CRN concentrations
    UpdateCRNs();
    
    //update morphology
    UpdateMorphology();
    
    //write output?
    if ((WriteResultsFlag != 0) && (Time <= WriteTime))
    {
      WriteProfile();
      WriteTime -= WriteInterval;
    }
    
		//update cliff position and time
		CliffPositionX -= RetreatRate*dt;
		Time -= dt;
		SeaLevel += SLR*dt;
		if (CliffPositionX < X[CliffPositionInd]) CliffPositionInd -= 1;
		//if (SeaLevel < Z[ZTrackInd]) ZTrackInd += 1;
    if (CliffPositionInd < 0) CliffPositionInd = 0;
    
	}
	
	WriteProfile();
	
  CliffPositionInd = 0;
    
	//Write result to file?
	if (WriteResultsFlag != 0)
	{
		ofstream write_results;
	  write_results.open((OutFileName+"_CRNs.txt").c_str());
	  write_results << "X, N\n";
	  for (int i=0; i<NXNodes; ++i) write_results << X[i] << " " << SurfaceN[i] << endl;
	  write_results.close();
  }
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
  vector<double> EmptyVector(NXNodes,0);
  P_Spal = EmptyVector;
  P_Muon = EmptyVector;
  
  // variation in production on the platform at depth as a function of a tidal period
	// First get a Tidal Sequence of Water Levels, clipping for negative depths
	for (int i=0, NN=WaterLevels.size(); i<NN;++i) WaterLevels[i] = SeaLevel+TideLevels[i];
		
  //LOOP ACROSS THE ACTIVE PART OF THE PLATFORM
	for (int i=CliffPositionInd; i<NXNodes; ++i)
	{
		//only work on active nodes
		if ((X[i] >= CliffPositionX) && (X[i] <= XMax))
		{
			//Get topographic shielding factor
			TopoShieldingFactor = GetTopographicShieldingFactor(X[i]-CliffPositionX, CliffHeight);

			//get water levels for this profile
			//reset production params
			P_Spal[i] = 0;
			P_Muon[i] = 0;

			for (int a=0, NN=WaterLevels.size(); a<NN; ++a)
			{
				if (WaterLevels[a] >= SurfaceElevation[i]) WaterDepths[a] = WaterLevels[a]-SurfaceElevation[i];
				else WaterDepths[a] = 0;

				//Calculate Production for this profile
				P_Spal[i] += GeoMagScalingFactor*TopoShieldingFactor*Po_Spal*exp(-WaterDepths[a]/z_ws);
				P_Muon[i] += TopoShieldingFactor*Po_Muon*exp(-WaterDepths[a]/z_wm);
			}

			//find mean production rate at surface
			P_Spal[i] /= NTidalValues;
			P_Muon[i] /= NTidalValues;
		
		  //loop through depths and update concentrations
		  int Top = 0;
		  for (int j=ZTrackInd; j<NZNodes; ++j)
		  {
			  //cout << Z[j] << " " << SurfaceElevation[i] << " " << i << " " << j << endl;
			  if ((Z[j] < PlatformElevation[i]) && (Z[j] > PlatformElevation[i]-20.))
			  {
				  if (Top == 0)
				  {	
					  //linearly interpolate to get concentration at the surface
					  if (PlatformElevationOld[i] == -9999) PlatformElevationOld[i] = SeaLevel+JunctionElevation;
					  else if (SurfaceN[i] > 0) SurfaceN[i] -= (((PlatformElevationOld[i]-PlatformElevation[i])/(PlatformElevationOld[i]-Z[j]))*(SurfaceN[i]-N[i][j]));
					  
					  //update concentration at the platform surface, accounting for beach cover
					  SurfaceN[i] += dt*P_Spal[i]*exp((0-((SurfaceElevation[i]-PlatformElevation[i])))/z_rs);
					  Top = 1;
				  }
				
				  //update concentrations at depth
				  //This is kept as SurfaceElevation not Platform Elevation for now
				  //NB This assumes that material density of the beach is the same as the bedrock!
				  N[i][j] += dt*P_Spal[i]*exp((Z[j]-SurfaceElevation[i])/z_rs);	//spallation
				  N[i][j] += dt*P_Muon[i]*exp((Z[j]-SurfaceElevation[i])/z_rm);	//muons
				  
				  //remove atoms due to radioactive decay
				  N[i][j] -= dt*Lambda;
			  }
		  }
		}
	}
}

void RockyCoastCRN::UpdateMorphology()
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
	for (int i=CliffPositionInd; i<NXNodes; ++i)
	{
		PlatformElevationOld[i] = PlatformElevation[i];
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
				  TempPlatformElevChange = dt*(RetreatRate*PlatformGradient-SLR);
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

double RockyCoastCRN::GetTopographicShieldingFactor(double X, double CliffHeight)
{
	/* 
	Function to calculate topographic shielding factor due to presence of a cliff.
	Cliff is assumed to be straight in profile, of fixed height, and vertical
	 
	Martin Hurst
	October 2014
	*/    

	if (X == 0) return 0.5;

	//##### Cliff Shielding paramters #####
	double d_theta_phi = (M_PI/180.)*5.0;			//azimuth and angle stepping
	double FMax = 2.0*M_PI*Po_Spal/(3.3);			//Maximum Intensity

	double Viewshed;
	double F = 0;

	for (double Az=-90; Az<=90.; Az+= 5.) 
	{
		Viewshed = atan(CliffHeight*cos((M_PI/180.)*Az)/X);
		F+= d_theta_phi*Po_Spal/(3.3)*pow(sin(Viewshed),3.3);
	}

	 return (FMax-F)/FMax;
		
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
	if (Time > GeoMagTime[GeoMagTime.size()-1]) return GeoMagScalingFactor;
	
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

  double Rate = 0.001;
  
	//interpolate to get value
	double Factor;
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

void RockyCoastCRN::GetSinWaveBeachWidth(double Time)
{
  // Function to return beach width where beach width varies as a sinosoidal 
  // of time.
  // Martin Hurst
  // Feb 11th 2016
  
  double MeanBeachWidth = 50;
  double BeachWidthAmp = 30;
  double BeachWidthWaveLength = 100;
  
  BeachWidth = MeanBeachWidth+BeachWidthAmp*sin((2.*M_PI*Time)/BeachWidthWaveLength);
}

void RockyCoastCRN::GetThinningBeachWidth(double Time)
{
  // Function to return beach width where beach width varies as a sinosoidal 
  // of time.
  // Martin Hurst
  // Feb 11th 2016
  
  double InitialBeachWidth = 100;
  double ThinTime = 1000;
  
  if (Time < ThinTime) BeachWidth = (Time/ThinTime)*InitialBeachWidth;
  else BeachWidth = InitialBeachWidth;
}

void RockyCoastCRN::WriteProfile()
{
  //Print to screen
	cout.flush();
	cout << "RockyCoastCRN: Time is " << setprecision(0) << fixed << Time << " years\r";
	
	//test if output file already exists
  string OutputFileName = OutFileName + "_Morphology.txt";
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
		for (int i=0; i<NXNodes; ++i) WritePlatform << setprecision(4) << " " << PlatformElevation[i];
		WritePlatform << endl;
		
		//write SurfaceElevation
		WritePlatform << Time;
		for (int i=0; i<NXNodes; ++i) WritePlatform << setprecision(4) << " " << SurfaceElevation[i];
		WritePlatform << endl;
	}
}
#endif
