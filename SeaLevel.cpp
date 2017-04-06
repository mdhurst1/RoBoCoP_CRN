/*------------------------------------------------------------------------

	SeaLevel.cpp
	
	Sea Level Object
	
	Martin D. Hurst, University of Glasgow, April 2017

	Copyright (C) 2017, Martin Hurst
	
	Developer can be contacted:
	martin.hurst@glasgow.ac.uk

	Martin D. Hurst
	
  
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

#ifndef SeaLevel_CPP
#define SeaLevel_CPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "SeaLevel.hpp"

using namespace std;

void SeaLevel::Initialise()
{
	// Initialising SeaLevel using pretend Milankovitch Cycles
	// Parameterisation of this is a complete fabrication to get
	// the right sort of pattern
	
	// Set the amplitude and frequency of each cyclicity
	double Obliquity_Amp = 5.;
	double Obliquity_Period = 41000.;
	double Obliquity_Offset = (double)rand()/RAND_MAX;
	
	double Precession_Amp = 8.;
	double Precession_Period = 26000.;
	double Precession_Offset = (double)rand()/RAND_MAX;
	
	double Eccentricity_Amp = 90.;
	double Eccentricity_Period = 100000.;
	double Eccentricity_Offset = 10000.; //set to roughly match end of Holocene
	
	
	//Set up Times vector
	// Calculate the relative contributions
	vector<double> Empty(710000,0.);
	Times = Empty;
	SeaLevels = Empty;
	for (int t=0, T=Times.size(); t<T; ++t) 
	{
		Times[t] = -700000.+t;
		// Obliquity signal
		SeaLevels[t] += Obliquity_Amp*cos(2.*M_PI*(Times[t]+Obliquity_Offset*Obliquity_Period)/Obliquity_Period);
		// Precession signal
		SeaLevels[t] += Precession_Amp*cos(2.*M_PI*(Times[t]+Precession_Offset*Precession_Period)/Precession_Period);
		// Eccentricity signal
		SeaLevels[t] += Eccentricity_Amp*cos(2.*M_PI*(Times[t]+Eccentricity_Offset*Eccentricity_Period)/Eccentricity_Period);
	}
}
void SeaLevel::Initialise(double SLR)
{
	// Initialising SeaLevel using a rate of Sea Level Rise
	// Set to zero for constant Sea Level
	
}
void SeaLevel::Initialise(string SeaLevelDataFile)
{
	// Initialising Sea Level from a datafile
	// Use this for SeaLevel output from other models or other RSL records
	double indata;
	char Dummy[32];
	
	//read in RSL data from Bradley et al. (2011) model
	ifstream SeaLevelIn(SeaLevelDataFile.c_str());
	if (!SeaLevelIn)
	{
	  printf("SeaLevel::%s: line %d Sea Level data file \"%s\" doesn't exist\n\n", __func__, __LINE__, SeaLevelFilename.c_str());
	  printf("Setting Realtive Sea Level to zero");
	}
	else
	{
		SeaLevelIn >> Dummy;
		SeaLevelIn >> Dummy;
		while (!SeaLevelIn.eof())
		{
		  SeaLevelIn >> indata;
		  Times.push_back(indata);
		  SeaLevelIn >> indata;
		  SeaLevels.push_back(indata);
		}
		SeaLevelIn.close();
	}
}
  		
void SeaLevel::get_SeaLevel(double Time)
{
	//interpolate to get value
	double Factor, Rate;
	int TimeCondition = 0;
	int ind = 0;
	while (TimeCondition == 0)
	{
		if (Time >= Times[Times.size()-1])
		{
			TimeCondition = 1;
			ind = RSLTime.size()-1;
		}
		else if (Time >= RSLTime[ind]) ++ind;
		else TimeCondition = 1;
	}
	Factor = (Time-Times[ind-1])/(Times[ind]-Times[ind-1]);
	SeaLevel = SeaLevels[ind-1]+Factor*(SeaLevels[ind]-SeaLevels[ind-1]);

	return SeaLevel;
}
