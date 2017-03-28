/*------------------------------------------------------------------------

	Hiro_Driver.cpp
	
	Driver file for running the shore platform model of Matsumoto et al. (2016)
	
	Matsumoto, H., Dickson, M. E., & Kench, P. S. (2016)
	An exploratory numerical model of rocky shore profile evolution. 
	Geomorphology, 268, 98-109. http://doi.org/10.1016/j.geomorph.2016.05.017
	
	Martin D. Hurst, University of Glasgow
	Hironori Matsumoto, University of Auckland
	
	February 2017
	
	Copyright (C) 2017, Martin Hurst
	
	Developer can be contacted:
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
#include "../Hiro.hpp"

using namespace std;

int main()
{
	//initialisation parameters
	double dZ = 0.1;
	double dX = 0.1;

	//Time control parameters
	double EndTime = 10000.;
	double Time = 0.;
	double TimeInterval = 1;

	//Print Control
	double PrintInterval = 100.;
	double PrintTime = Time;
	string OutputFileName = "ShoreProfile.xz";
	
	//initialise Hiro Model
	Hiro PlatformModel = Hiro(dZ, dX);
	
//	//Initialise Tides
// //This will become setting up the erosion shape function
//	double TidalAmplitude = 1.;
//	double TidalPeriod = 12.42;
//	PlatformModel.InitialiseTides(TidalAmplitude, TidalPeriod);
	double TidalRange = 1;
	PlatformModel.InitialiseTides(TidalRange);

	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight_Mean = 2.;
	double WaveHeight_StD = 0;
	double WavePeriod_Mean = 6.;
	double WavePeriod_StD = 0;
	PlatformModel.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);

	//Sea level rise?
	double SLR = 0;
	PlatformModel.InitialiseSeaLevel(SLR);
	

	//Loop through time
	while (Time <= EndTime)
	{
		//Update Sea Level
		PlatformModel.UpdateSeaLevel();

		//Get the wave conditions
		PlatformModel.GetWave();

		//Calculate forces acting on the platform
		PlatformModel.CalculateBackwearing();
		PlatformModel.CalculateDownwearing();

		//Do erosion
		PlatformModel.ErodeBackwearing();
		PlatformModel.ErodeDownwearing();

		//Implement Weathering
		PlatformModel.IntertidalWeathering();

		//Update the Morphology 
		PlatformModel.UpdateMorphology();	  
		
		//Check for Mass Failure
		PlatformModel.MassFailure();
		
		//print?
		if (Time >= PrintTime)
		{
			PlatformModel.WriteProfile(OutputFileName, Time);
			PrintTime += PrintInterval;
			cout << endl;
		}

		//update time
		Time += TimeInterval;
	}
	
	string ResistanceFileName = "ResistanceArray.xz";
	string MorphologyFileName = "MorphologyArray.xz";
	PlatformModel.WriteResistanceArray(ResistanceFileName, Time);
	PlatformModel.WriteMorphologyArray(MorphologyFileName, Time);
	
	//a few blank lines to finish
	cout << endl << endl;
	
}


