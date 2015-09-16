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

#ifndef RockyCoastCRN_CPP
#define RockyCoastCRN_CPP

void RockyCoastCRN::Initialise()
{
  /* initialise an empty platform object
	retreatrate is the rate of cliff retreat (m/yr)
	beachwidth is the width of the beach (m), protecting the platform from CRN production
	elevinit is the elevation of the platform at the beach edge 
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
	ElevInit = NDV;
	TidalAmplitude = NDV;
}
void RockyCoastCRN::Initialise(double retreatrate, double beachwidth, double elevinit, double platformgradient, double cliffheight, double tidalamplitude)
{
	/* initialise a platform object for single retreat rate scenarios
	retreatrate is the rate of cliff retreat (m/yr)
	beachwidth is the width of the beach (m), protecting the platform from CRN production
	elevinit is the elevation of the platform at the beach edge 
	platformgradient is the average slope of the platform surface (m/m) 
	cliffheight is the height of the adjacent cliff (m)
	tidalamplitude is the average tidal amplitude for diurnal tides. */

	//Geometric paramters
	NXNodes = 201;											//Number of nodes in cross shore
	NZNodes = 201;											//Number of nodes in profile

	PlatformWidth = 1000.;						//width of model domain (m)
	PlatformDepth = 40.;							//Depth to which CRN will be tracked (needs to be large enough that sea level rise is ok)
	NDV = -9999;											//Place holder for no data
	
	//setup vectors
	vector<double> EmptyX(NXNodes,0.0);
	vector<double> EmptyZ(NZNodes,0.0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	vector< vector<double> > EmptyVV(NXNodes,EmptyZ);
	X = EmptyX;
	Z = EmptyZ;
	SurfaceElevation = EmptyXNDV;
	SurfaceElevationOld = EmptyXNDV;
	SurfaceN = EmptyX;
	N = EmptyVV;
	for (int j=0; j<NZNodes; ++j) Z[j] = (20.-j*(PlatformDepth/(NZNodes-1)));
	for (int i=0; i<NXNodes; ++i) X[i] = (i*(PlatformWidth/(NXNodes-1)));
	
	//Assign parameters
	RetreatRate1 = retreatrate;
	RetreatRate2 = retreatrate;
	ChangeTime = 0;
	PlatformGradient = platformgradient;
	CliffHeight = cliffheight;
	BeachWidth = beachwidth;
	ElevInit = elevinit;
	TidalAmplitude = tidalamplitude;
	
	/// TIDES 
	TidalPeriod=12.;
	for (double TT = 0; TT <= TidalPeriod; TT += 0.2) TideLevels.push_back(-TidalAmplitude*sin((2.*M_PI*TT)/(TidalPeriod)));
	NTidalValues = (double)TideLevels.size();
	WaterDepths.resize(NTidalValues);
	WaterLevels.resize(NTidalValues);
	
	/// GEOMAG VARIATION AND RELATIVE SEALEVEL
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

void RockyCoastCRN::Initialise(double retreatrate1, double retreatrate2, double changetime, double beachwidth, double platformgradient, double cliffheight, double elevinit, double tidalamplitude)
{
  /* initialise a platform object for two retreat rate scenarios
  retreatrate1 runs from 7.5ka until changetime, after which retreatrate2 continues to present
	retreatrate1 is the first rate of cliff retreat (m/yr)
	retreatrate2 is the 2nd rate of cliff retreat (m/yr)
	changetime is the time to switch from retreatrate1 to retreatrate2 (years BP)
	beachwidth is the width of the beach (m), protecting the platform from CRN production
	elevinit is the elevation of the platform at the beach edge 
	platformgradient is the average slope of the platform surface (m/m) 
	cliffheight is the height of the adjacent cliff (m)
	tidalamplitude is the average tidal amplitude for diurnal tides. */
	
	//Geometric parameters
	NXNodes = 201;											//Number of nodes in cross shore
	NZNodes = 201;											//Number of nodes in profile

	PlatformWidth = 1000.;						//width of model domain (m)
	PlatformDepth = 40.;							//Depth to which CRN will be tracked (needs to be large enough that sea level rise is ok)
	NDV = -9999;											//Place holder for no data
	
	//setup vectors
	vector<double> EmptyX(NXNodes,0.0);
	vector<double> EmptyZ(NZNodes,0.0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	vector< vector<double> > EmptyVV(NXNodes,EmptyZ);
	X = EmptyX;
	Z = EmptyZ;
	SurfaceElevation = EmptyXNDV;
	SurfaceElevationOld = EmptyXNDV;
	SurfaceN = EmptyX;
	N = EmptyVV;
	for (int j=0; j<NZNodes; ++j) Z[j] = (20.-j*(PlatformDepth/(NZNodes-1)));
	for (int i=0; i<NXNodes; ++i) X[i] = (i*(PlatformWidth/(NXNodes-1)));
	
	//Assign Parameters
	RetreatRate1 = retreatrate1;
	RetreatRate2 = retreatrate2;
	ChangeTime = changetime;
	PlatformGradient = platformgradient;
	CliffHeight = cliffheight;
	BeachWidth = beachwidth;
	ElevInit = elevinit;
	TidalAmplitude = tidalamplitude;
	
	/// TIDES
	TidalPeriod=12.;
	for (double TT = 0; TT <= TidalPeriod; TT += 0.5) TideLevels.push_back(-TidalAmplitude*sin((2.*M_PI*TT)/(TidalPeriod)));
	NTidalValues = (double)TideLevels.size();
	WaterDepths.resize(NTidalValues);
	WaterLevels.resize(NTidalValues);
	
	/// GEOMAG VARIATION AND RELATIVE SEALEVEL
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
	GeoMagTime.pop_back();
	GeoMagScalingFactors.pop_back();
	GeoMagIn.close();
	
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
	RSLTime.pop_back();
	RSLRate.pop_back();
	GIAIn.close();
}

void RockyCoastCRN::UpdateParameters( double RetreatRate1_Test, double RetreatRate2_Test, 
                                    double ChangeTime_Test, double BeachWidth_Test, 
                                    double ElevInit_Test)
{
  /*  Function to update model parameters for a new model run. 
  This is designed for use with the MCMC_RockyCoast object in order to save a 
  bit of time on reinitialising each time in the Markov chain 
  RetreatRate1_Test is the new RetreatRate1
  RetreatRate2_Test is the new RetreatRate2
  ChangeTime_Test is the new ChangeTime
  BeachWidth_Test is the new BeachWidth
  ElevInit_Test is the new ElevInit    */
  
  RetreatRate1 = RetreatRate1_Test;
  RetreatRate2 = RetreatRate2_Test;
  ChangeTime = ChangeTime_Test;
  BeachWidth = BeachWidth_Test;
  ElevInit = ElevInit_Test;

}

void RockyCoastCRN::RunModel(int RetreatType, int WriteResultsFlag)
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
  
  //Reset Concentrations!
  //setup vectors
	vector<double> EmptyX(NXNodes,0.0);
	vector<double> EmptyZ(NZNodes,0.0);
	vector<double> EmptyXNDV(NXNodes,NDV);
	vector< vector<double> > EmptyVV(NXNodes,EmptyZ);
	
	SurfaceElevation = EmptyXNDV;
	SurfaceElevationOld = EmptyXNDV;
	SurfaceN = EmptyX;
	N = EmptyVV;
	  
	double SLR = 0;					//Rate of relative sea level rise (m/y)
		
	//Time control
	double Time = 0;					//time in years
	double dt = 4;						//time step
	double SeaLevel = 0;			//Sea Level Tracker
	double XMax = 1000.; 	    //Distance from cliff

	//Work out start time from assumed cliff retreat rates.
	double MaxTime;
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
	if (MaxTime > 7000)	Time = 7000;
	else Time = MaxTime;
	MaxTime = Time;
	
	double XTrack;
	if (RetreatType == 0) XTrack = MaxTime*RetreatRate1;
	else if (RetreatType == 1) XTrack = ChangeTime*RetreatRate2 + (Time-ChangeTime)*RetreatRate1; 		      //for tracking cliff posiiton in X
	else if (RetreatType == 2) XTrack = MaxTime*(RetreatRate1+RetreatRate2)/2.;
	else 
	{
	  cout << "Retreat Type Unknown" << endl;
	  exit(EXIT_SUCCESS);
	}
	
	XMax = XTrack;
	
	//Get surface elevations
	int XTrackInd = 0;
	int ZTrackInd = 0;
	for (int i=0; i<NXNodes; ++i)
	{
	  if (X[i] > XTrack) SurfaceElevation[i] = ElevInit-(X[i]-(XTrack+BeachWidth))*PlatformGradient;
	  else XTrackInd = i;
	}
  SurfaceElevationOld = SurfaceElevation;
  
  //Set Z tracker
  for (int j=0; j<NZNodes; ++j) if (Z[j] > ElevInit) ZTrackInd = j;

	double P_Spal, P_Muon;
	double GeoMagScalingFactor, TopoShieldingFactor;
	double SurfaceElevChange;
	double RetreatRate;
	double dTime,Factor,dRR;

	//////////////////////////////////////////////////////////
	//MAIN MODEL LOOP
	//////////////////////////////////////////////////////////

	while (Time > 0)
	{
		//Get retreat rate depending on style of retreat
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
		
		//Get shielding and scaling factors
		GeoMagScalingFactor = GetGeoMagScalingFactor(Time);
		SLR = GetSeaLevelRise(Time);
		
		// variation in production on the platform at depth as a function of a tidal period
		// First get a Tidal Sequence of Water Levels, clipping for negative depths
		for (int i=0, NN=WaterLevels.size(); i<NN;++i) WaterLevels[i] = SeaLevel+TideLevels[i];
		
		//LOOP ACROSS THE ACTIVE PART OF THE PLATFORM
		for (int i=XTrackInd; i<NXNodes; ++i)
		{
			//only work on active nodes beyond the beach
			if ((X[i] > XTrack+BeachWidth) && (X[i] <= XMax))
			{
				//Get topographic shielding factor
				TopoShieldingFactor = GetTopographicShieldingFactor(X[i]-XTrack, CliffHeight);

				//get water levels for this profile
				//reset production params
				P_Spal = 0;
				P_Muon = 0;

				for (int a=0, NN= WaterLevels.size(); a<NN; ++a)
				{
					if (WaterLevels[a] >= SurfaceElevation[i]) WaterDepths[a] = WaterLevels[a]-SurfaceElevation[i];
					else WaterDepths[a] = 0;

					//Calculate Production for this profile
					P_Spal += GeoMagScalingFactor*TopoShieldingFactor*Po_Spal*exp(-WaterDepths[a]/z_ws);
					P_Muon += TopoShieldingFactor*Po_Muon*exp(-WaterDepths[a]/z_wm);
				}

				//find mean production rate at surface
				P_Spal /= NTidalValues;
				P_Muon /= NTidalValues;
			
			  //loop through depths and update concentrations
			  int Top = 0;
			  for (int j=ZTrackInd; j<NZNodes; ++j)
			  {
				  //cout << Z[j] << " " << SurfaceElevation[i] << " " << i << " " << j << endl;
				  if ((Z[j] <= SurfaceElevation[i]) && (Z[j] > SurfaceElevation[i]-20.))
				  {
					  if (Top == 0)
					  {	
						  //linearly interpolate to get concentration at the surface
						  if (SurfaceElevationOld[i] == -9999) SurfaceElevationOld[i] = SeaLevel+ElevInit;
						  else if (SurfaceN[i] > 0) SurfaceN[i] -= (((SurfaceElevationOld[i]-SurfaceElevation[i])/(SurfaceElevationOld[i]-Z[j]))*(SurfaceN[i]-N[i][j]));
						  
						  //update concentration at the surface
						  SurfaceN[i] += dt*P_Spal*exp((0-((SurfaceElevationOld[i]-SurfaceElevation[i])))/z_rs);
						  Top = 1;
					  }
					
					  //update concentrations at depth
					  N[i][j] += dt*P_Spal*exp((Z[j]-SurfaceElevation[i])/z_rs);	//spallation
					  N[i][j] += dt*P_Muon*exp((Z[j]-SurfaceElevation[i])/z_rm);	//muons
					  
					  //remove atoms due to radioactive decay
					  N[i][j] -= dt*Lambda;
				  }
			  }
			}
		}
		
		//update position and time
		XTrack -= RetreatRate*dt;
		Time -= dt;
		SeaLevel += SLR*dt;
		if (XTrack < X[XTrackInd]) XTrackInd -= 1;
		if (SeaLevel < Z[ZTrackInd]) ZTrackInd += 1;
    if (XTrackInd < 0) XTrackInd = 0;

		for (int i=XTrackInd; i<NXNodes; ++i)
		{
			SurfaceElevationOld[i] = SurfaceElevation[i];
			if ((X[i] >= XTrack) && (XTrack <= XMax))
			{
				if (SurfaceElevation[i] == NDV) SurfaceElevation[i] = SeaLevel+ElevInit-(X[i]-XTrack)*PlatformGradient;
				else
				{
					SurfaceElevChange = dt*(RetreatRate*PlatformGradient-SLR);
					if (SurfaceElevChange < 0) SurfaceElevChange = 0;
					SurfaceElevation[i] -= SurfaceElevChange;
				}
			}
		}
	}
	
  XTrackInd = 0;
    
	//Write result to file?
	if (WriteResultsFlag != 0)
	{
		char FileName[128];
	  sprintf(FileName, "CRN_Model_Results_R1%4.4f_R2%4.4f_t%4.4f_Beta%4.4f_T%1.1f.txt", RetreatRate1, RetreatRate2, ChangeTime, PlatformGradient, TidalAmplitude);
	  ofstream write_results;
	  write_results.open(FileName);
	  write_results << "X, N\n";
	  for (int i=0; i<NXNodes; ++i) write_results << X[i] << " " << SurfaceN[i] << endl;
	  write_results.close();
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

#endif
