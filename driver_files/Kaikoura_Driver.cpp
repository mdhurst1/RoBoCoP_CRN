/*------------------------------------------------------------------------

	Kaikoura_Driver.cpp
	
	DC++ implementation of Hiro Matsumoto's Shore Platform Model

	Matsumoto, H., Dickson, M. E., & Kench, P. S. (2016a)
	An exploratory numerical model of rocky shore profile evolution. 
	Geomorphology, 268, 98-109. http://doi.org/10.1016/j.geomorph.2016.05.017
	
	Matsumoto, H., Dickson, M.E., and Kench, P.S. (2016b)
	Modelling the Development of Varied Shore Profile Geometry on Rocky Coasts.
	Journal of Coastal Research, 75, 597–601, http://dx.doi.org/10.2112/SI75-120.1

	Martin D. Hurst, University of Glasgow
	Hironori Matsumoto, University of Auckland
	
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
#include <Hiro.hpp>

using namespace std;

int main()
{
	//initialisation parameters
	double dZ = 0.1;
	double dX = 0.1;
	double Gradient = 0.5;
	double CliffHeight = 10.;
	
	//Time control parameters
	double EndTime = 10000;
	double Time = 0.;
	double TimeInterval = 1;

//	// Do an earthquake
	double EarthquakeTime = 5000.;
	double EarthquakeSeaLevel = -1.;

	//Print Control
	double PrintInterval = 1.;
	double PrintTime = Time;
	string OutputFileName = "Kaikoura_ShoreProfiles_1m.xz";
	
	//initialise Hiro Model
	Hiro PlatformModel = Hiro(dZ, dX, Gradient, CliffHeight);
	
	//Initialise Tides
	double TidalRange = 2.;
	PlatformModel.InitialiseTides(TidalRange);
		
	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight_Mean = 1.;
	double WaveHeight_StD = 0.;
	double WavePeriod_Mean = 6.;
	double WavePeriod_StD = 0;
	PlatformModel.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);

	//Sea level rise?
	double SLR = 0;
	PlatformModel.InitialiseSeaLevel(SLR);

	// Wave coefficient constant
	double StandingCoefficient = 0.01;
	double BreakingCoefficient = 1.;
	double BrokenCoefficient = 1.;
	PlatformModel.Set_WaveCoefficients(StandingCoefficient, BreakingCoefficient, BrokenCoefficient);

	//reset the geology
	double CliffFailureDepth = 1.;
	double Resistance = 0.001;
	double WeatheringRate = 0.0001;
	PlatformModel.InitialiseGeology(CliffHeight, CliffFailureDepth, Resistance, WeatheringRate);
				

	//Loop through time
	while (Time <= EndTime)
	{
		//Update Sea Level
		if (Time > EarthquakeTime)
		{
			PlatformModel.UpdateSeaLevel(EarthquakeSeaLevel);
			PlatformModel.UpdateMorphology();
			EarthquakeTime = 0.;
		}

		//Get the wave conditions
		PlatformModel.GetWave();

		//Calculate forces acting on the platform
		PlatformModel.CalculateBackwearing();
		//PlatformModel.CalculateDownwearing();

		//Do erosion
		PlatformModel.ErodeBackwearing();
		//PlatformModel.ErodeDownwearing();

		//Update the Morphology 
		PlatformModel.UpdateMorphology();	  
		
		//Implement Weathering
		PlatformModel.IntertidalWeathering();
		
		//Update the Morphology 
		PlatformModel.UpdateMorphology();

		//Check for Mass Failure
		PlatformModel.MassFailure();
		
		//Update the Morphology 
		PlatformModel.UpdateMorphology();
				
		//print?
		if (Time >= PrintTime)
		{
			PlatformModel.WriteProfile(OutputFileName, Time);
			PrintTime += PrintInterval;
			//cout << endl;
		}
		
		//update time
		Time += TimeInterval;
		
	}
	
	// Need to implement a WRITE METADATA FUNCTION!
	// GIVE EACH MODEL RUN A UNIQUE ID NUMBER?
	// OR INPUT FILES TO CONTROL PARAMETERS?
	
	//a few blank lines to finish
	cout << endl << endl;
	
}


