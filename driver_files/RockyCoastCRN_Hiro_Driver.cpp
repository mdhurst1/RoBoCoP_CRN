/*==============================================================

RockyCoastCRN_Hiro_Driver.hpp

A driver function to the concentrations of 10Be in a shore platform
with morphological evolution read from a text file. Developed to allow
morphological evolution to be simulated by soft-coupling to the model
of Hiro Matsumoto et al. (2016).

Developed by:
Martin D. Hurst

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

==============================================================*/
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
#include "../RockyCoastCRN.hpp"

using namespace std;

int main(int nNumberofArgs,char *argv[])
{
	//print heading to screen
	cout << endl << "RockyCoastCRN: Hiro Driver" << endl;
	
	//Check input arguments
	if (nNumberofArgs!=3)
	{
		cout << "The program requires two inputs" << endl;
		cout << "\t- The input morphology file name" << endl;
		cout << "\t- The tidal range (in metres)" << endl;
		exit(EXIT_SUCCESS);
	}
	
	//Parse arguments
	string FileName = argv[1];
	float TidalRange = argv[2];
	float TidalAmp = 0.5*TidalRange;
	
	//Input parameters
	float CliffHeight = 0.;             //Cliff height for shielding
	
	//Initialise Tides
	PlatformCRN.InitialiseTides(TidalAmplitude, TidalPeriod);
	float TidalPeriod = 12.42;
	
	// Time Control
	float EndTime = 6000.;
	float Time = 0.;
	float TimeInterval = 1.;
	
	// Declare X and Z
	vector<double> X;
	vector<double> Z;
	for (float i=0; i<=750; ++i) Z[i] = 37.5-i*0.1;
	
	// Read time series of coastal change from file
	// File format is each row represents 1 year, there are 750 columns each representing 
	// an elevation from 37.5 to -37.5 m in 0.1 m increments. The value in each represents
	// the horizontal X position of the coast
	ifstream ReadProfileIn;
	ReadProfileIn.open(FileName.c_str());
	// declare a string for the line and individual value
	string Line;
	// declare a float to hold the newly read value
	float Value;
	// declare a string stream for converting the line into values
	stringstream PositionString;
	
	//Loop through time
	while (Time < EndTime)
	{
		//Get a new morphology from input file
		X.clear();
		getline(ReadProfileIn, Line);	
		PositionString(Line);
		while (PositionString >> Value) X.pushback(Value);
		
		//Update the morphology inside RockyCoastCRN
		PlatformCRN.UpdateMorphology(X,Z);

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
}


