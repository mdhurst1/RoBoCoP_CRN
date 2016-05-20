/*------------------------------------------------------------------------

	RoBoCop_Driver.cpp
	
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

using namespace std;

int main()
{
  //initialisation parameters
  double dZ = 0.1;
  //double PlatformGradient = 1./50.;
  //double CliffPositionX = 0.;
  
  //Rate of sea level rise
  double SLR = 0.0002;
  
  //initialise RoBoCoP
  RoBoCoP PlatformModel = RoBoCoP(dZ);
	
	//Initialise Tides
	double TidalAmplitude = 1.;
	double TidalPeriod = 12.42;
	PlatformModel.InitialiseTides(TidalAmplitude, TidalPeriod);
	
	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight = 1.2;
	double WavePeriod = 6.;
	PlatformModel.InitialiseWaves(WaveHeight,WavePeriod);
	
	// Time Control
	double EndTime = 10000.;
	double Time = 0.;
	double TimeInterval = 1;
	
	//Print Control
	double PrintInterval = 100.;
	double PrintTime = Time;
	string OutputFileName = "ShoreProfile.xz";
	string ErosionFileName = "Erosion.ez";
	
	//Loop through time
	while (Time < EndTime)
	{
	  //Update Sea Level
	  PlatformModel.UpdateSeaLevel(SLR,TimeInterval);
	  
	  //Evolve the coast
	  PlatformModel.EvolveCoast(TimeInterval);
	  
	  //print?
	  if (Time >= PrintTime)
	  {
	    PlatformModel.WriteProfile(OutputFileName, Time);
	    PlatformModel.WriteErosion(ErosionFileName, Time);
	    PrintTime += PrintInterval;
	    cout << "Time is " << Time << endl;
	  }
	  
	  //update time
	  Time += TimeInterval;
	}
		
	//a few blank lines to finish
	cout << endl << endl;
	
}


