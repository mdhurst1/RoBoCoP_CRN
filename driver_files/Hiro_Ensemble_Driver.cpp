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

template <typename T> string tostr(const T& t)
{ 
   ostringstream os; 
   os<<t; 
   return os.str(); 
}

int main()
{
	//initialisation parameters
	double dZ = 0.1;
	double dX = 0.1;

	//initialise Hiro Model
	Hiro PlatformModel = Hiro(dZ, dX);
	
	//Initialise Tidal Ranges to check
	vector<double> TidalRanges;
	TidalRanges.push_back(1.);
	TidalRanges.push_back(8.);
	
	//Initialise WaveHeights to check
	vector<double> WaveHeights;
	WaveHeights.push_back(1.);
	WaveHeights.push_back(3.);
	
	//Initialise weathering efficacy to check
	vector<double> WeatheringRates;
	WeatheringRates.push_back(0.005);
	WeatheringRates.push_back(0.5);
	
	//Geology
	double CliffHeight = 10;
	double CliffFailureDepth = 10;
	double RockResistance = 1; 
	double WeatheringConst = 0.005;

	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight_Mean = 2.;
	double WaveHeight_StD = 0;
	double WavePeriod_Mean = 6.;
	double WavePeriod_StD = 0;
	
	
	//Loop across parameter space
	for (int a=0, A=TidalRanges.size(); a<A; ++a)
	{
		for (int b=0, B=WaveHeights.size(); b<B; ++b)
		{
			for (int c=0, C=WeatheringRates.size(); c<C; ++c)
			{
				//setup the output file
				PlatformModel.OutputFileName = "ShoreProfile_T"+tostr(TidalRanges[a])+"_H"+tostr(WaveHeights[b])+"_W"+tostr(WeatheringRates[c])+".xz";
		
				//reset the model
				PlatformModel.ResetModel();
		
				//Reset tidal range
				double TidalRange = TidalRanges[a];
				PlatformModel.InitialiseTides(TidalRange);
		
				//reset wave height
				WaveHeight_Mean = WaveHeights[b];
				PlatformModel.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);
		
				//reset the geology
				WeatheringConst = WeatheringRates[c];
				PlatformModel.InitialiseGeology(CliffHeight, CliffFailureDepth, RockResistance, WeatheringConst);
				
				//run the model!
				cout << setprecision(1);
				cout << "Running Model with..." << endl;
				cout << "\tTidal Range " << TidalRange << " m" << endl;
				cout << "\tWave Height " << WaveHeight_Mean << " m" << endl;
				cout << "\tMax Weathering Rate " << WeatheringConst << " m/yr?" << endl;
				PlatformModel.EvolveCoast();
				cout << "\nDone\n";
			}
		}
	}	
	//a few blank lines to finish
	cout << endl << endl;
	
}


