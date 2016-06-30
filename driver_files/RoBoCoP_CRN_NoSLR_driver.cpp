/*------------------------------------------------------------------------

	RoBoCop_CRN_Driver.cpp
	
	ROck and BOttom COastal Profile (RoBoCop) Model
	
	Martin D. Hurst, British Geological Survey, February 2016

	Copyright (C) 2016, Martin Hurst
	
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
#include <unistd.h>
#include "../RoBoCoP.hpp"
#include "../RockyCoastCRN.hpp"

using namespace std;

int main()
{
  //initialisation parameters
  double dZ = 0.1;
  double PlatformGradient = 1./10.;
  double CliffPositionX = 0.;
  
  //Rate of sea level rise
  double SLR = 0.0;
  
  	// Time Control
	double EndTime = 10000.;
	double Time = 0.;
	double TimeInterval = 1.;

  //initialise RoBoCoP
  RoBoCoP PlatformModel = RoBoCoP(dZ,PlatformGradient,CliffPositionX,TimeInterval);

  //Initialise RockyCoastCRN
	RockyCoastCRN PlatformCRN = RockyCoastCRN(PlatformModel);
	
	//Initialise Tides
	double TidalAmplitude = 1.0;
	double TidalPeriod = 12.42;
	PlatformModel.InitialiseTides(TidalAmplitude, TidalPeriod);
	PlatformCRN.InitialiseTides(TidalAmplitude, TidalPeriod);
	
	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double MeanWaveHeight = 1.0;
	double StdWaveHeight = 0.0;
	double MeanWavePeriod = 6.;
	double StdWavePeriod = 0.;
	PlatformModel.InitialiseWaves(MeanWaveHeight, StdWaveHeight, MeanWavePeriod, StdWavePeriod);
	
	//Print Control
	double PrintInterval = 100.;
	double PrintTime = PrintInterval;
	string OutputMorphologyFileName = "ShoreProfile.xz";
	string OutputConcentrationFileName = "CRNConcentrations.xn";
	
	//print initial conditions to file first
	PlatformModel.WriteProfile(OutputMorphologyFileName, Time);
  PlatformCRN.WriteCRNProfile(OutputConcentrationFileName, Time);
  
	//Loop through time
	while (Time < EndTime)
	{
	  //Update Sea Level
	  PlatformModel.UpdateSeaLevel(SLR);
	  
	  //Get a new wave
	  PlatformModel.GetWave();
	  
	  //Evolve the coast
	  PlatformModel.EvolveCoast();
	  
	  //Update the morphology inside RockyCoastCRN
	  PlatformCRN.UpdateMorphology(PlatformModel);
	  
	  //Update the CRN concentrations
	  PlatformCRN.UpdateCRNs();
    
    //update time
	  Time += TimeInterval;
	      
	  //print?
	  if (Time >= PrintTime)
	  {
	    PlatformModel.WriteProfile(OutputMorphologyFileName, Time);
	    PlatformCRN.WriteCRNProfile(OutputConcentrationFileName, Time);
	    PrintTime += PrintInterval;
	  }
	  //cout << "Time is " << Time << endl;
	}
	//write a gap at the end
	cout << endl << endl;
}


