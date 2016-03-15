/*------------------------------------------------------------------------

	RoBoCop.cpp
	
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

#ifndef RoBoCoP_CPP
#define RoBoCoP_CPP

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
#include "RoBoCoP.hpp"

using namespace std;

void RoBoCoP::Initialise()
{
  /* initialise an empty RoBoCoP object as default */
  printf("\nRoBoCoP.Initialise: Initialised an empty RoBoCoP object\n");
}

void RoBoCoP::Initialise(double dZ)
{
  /* initialise a verictal cliff RoBoCoP object */
  printf("\nRoBoCoP.Initialise: Initialised a RoBoCoP as a vertical cliff\n");
  
  //Declare stuff
  NoNodes = round(20./dZ)+1;
  NDV = -9999;
  
  //declare an array of zeros for assignment of x and z vectors
  vector<double> Empty(NoNodes,0.0);
	X = Empty;
	Z = Empty;
	
	//Loop through array and calculate X and Z
	for (int i=0; i<NoNodes; ++i) 
	{
	  Z[i] = 10.-i*dZ;
	  X[i] = 0;
  }
}

void RoBoCoP::Initialise(double dZ, double PlatformGradient)
{
	/* initialise a sloped platform RoBoCoP object */
	printf("\nRoBoCoP.Initialise: Initialised a RoBoCoP as a planar, sloped platform\n");
	
  //Declare stuff
  NoNodes = 1+round(20./dZ);
  NDV = -9999;
  
  //declare an array of zeros for assignment of x and z vectors
  vector<double> Empty(NoNodes,0.0);
	X = Empty;
	Z = Empty;
	
	//Loop through array and calculate X and Z
	for (int i=0; i<NoNodes; ++i) 
	{
	  Z[i] = 10.-i*dZ;
	  X[i] = -Z[i]/PlatformGradient;
  }
}

void RoBoCoP::Initialise(double dZ, double PlatformGradient, double CliffPositionX)
{
	/* initialise a sloped platform backed by a vertical cliff RoBoCoP object */
  printf("\nRoBoCoP.Initialise: Initialised a RoBoCoP as planar, sloped platform backed by a cliff\n");
  
  //Declare stuff
  NoNodes = 1+round(20./dZ);
  NDV = -9999;
  
  //declare an array of zeros for assignment of x and z vectors
  vector<double> Empty(NoNodes,0.0);
	X = Empty;
	Z = Empty;
	
	//Loop through array and calculate X and Z
	for (int i=0; i<NoNodes; ++i) 
	{
	  Z[i] = 10.-i*dZ;
	  if (Z[i] < 0) X[i] = CliffPositionX-(Z[i]/PlatformGradient);
	  else X[i] = CliffPositionX;
  }
}	

void RoBoCoP::Initialise(string FileNameIn)
{
  /* read coastal profile from file to RoBoCoP object */
  
  //Declare stuff
  NoNodes = 1+round(20./dZ);
  NDV = -9999;
}

void RoBoCoP::InitialiseTides(double TidalAmplitude, double TidalPeriod)
{
  /* intialise the tides as a simple cosine wave (i.e. simple diurnal/semidiurnal) */

  //Tidal Timestep (in hours)
  double dT = 0.01;    
  double TideTime;
  
  //Setup tide water levels array
  NTideValues = (int)(TidalPeriod/(2.*dT));
  vector<double> Empty(NTideValues,-9999);
  WaterLevels = Empty;
  
  //loop through and calculate water levels
  for (int i=0, N=WaterLevels.size(); i<N; ++i)
  {
    TideTime = i*dT;
    WaterLevels[i] = -TidalAmplitude*cos(2.*M_PI*TideTime/TidalPeriod);
  }  
}

void RoBoCoP::InitialiseWaves(double OffshoreWaveHeight, double WavePeriod)
{
  /* intialise waves as a single wave */

  WavePeriod = WavePeriod;
  
  //Breaking Wave Height calculated following Komar and Gaughan (1972)
  BreakingWaveHeight = 0.39*pow(g,0.2)*pow(WavePeriod,0.4)*pow(OffshoreWaveHeight,2.4);
  
  //Water depth of breaking wave
  BreakingWaveWaterDepth = BreakingWaveHeight/0.78;
}

void RoBoCoP::EvolveCoast(double TimeInterval)
{
  /* Function to evolve the coastal profile through time following
      Trenhaile (2000)  */
  
  //declarations
  double SurfZoneTopZ, SurfZoneBottomZ, SurfZoneTopX, SurfZoneBottomX, SurfZoneWidth, SurfForce, WaterDepth;
  
  //vector to sum erosion
  vector<double> Empty(NoNodes,0);
  Erosion = Empty;
  
  //Determine surf zone width
  for (int t=0, T=WaterLevels.size(); t<T; ++t)
  {
    //Elevation of surf zone limits
    SurfZoneBottomZ = SeaLevel+WaterLevels[t]-BreakingWaveWaterDepth;
    SurfZoneTopZ = SeaLevel+WaterLevels[t];
    
    //find horizontal extent of surfzone
    bool SurfZone=false;
    int i=NoNodes-1;
    while (SurfZone == false)
    {
      if (Z[i] > SurfZoneBottomZ)
      {
        //Find Position X of Surfzone bottom by interpolating
        SurfZoneBottomX = X[i] + (X[i+1]-X[i])*(Z[i]-SurfZoneBottomZ)/(Z[i+1]-Z[i]);
        SurfZone = true;
      }
      --i;
    }
    while (SurfZone == true)
    {
      if (Z[i] > SurfZoneTopZ)
      {
        //Find Position X of Surfzone bottom by interpolating
        SurfZoneTopX = X[i] + (X[i+1]-X[i])*(Z[i]-SurfZoneTopZ)/(Z[i+1]-Z[i]);
        SurfZone = false;
        ++i;
      }
      --i;
    }
    
    //Calculate Surf Zone width
    SurfZoneWidth = SurfZoneBottomX-SurfZoneTopX;
    
    //Calculate Surf Force reaching the water line
    SurfForce = 0.5*rho_w*(BreakingWaveWaterDepth)*exp(-k*SurfZoneWidth);
    
    //Loop offshore and calculate erosion
    while (i < NoNodes)
    {
      WaterDepth = SeaLevel+WaterLevels[t]-Z[i];
      if (fabs(WaterDepth)<0.0001) WaterDepth=0;
      else if (WaterDepth < 0) WaterDepth = 100.;

      Erosion[i] += M*SurfForce*exp(-WaterDepth);
      ++i;
    }
  }
  //Do the erosion/update the profile
  double XCliff = X[NoNodes-1];
  for (int i=NoNodes-1; i>-1; --i)
  {
    X[i] -= Erosion[i]*TimeInterval;
    
    //check for overhangs and remove
    if (X[i] > XCliff) X[i] = XCliff;
    XCliff = X[i];
  }
}

void RoBoCoP::WriteProfile(string OutputFileName, double Time)
{
  /* Writes a coastline object X and Z coordinates to file for a given time
		If the file already exists the data will be ed else a new file is
		created.

		File format is 	Time | X[0] | X[1] | X[2] =====> X[NoNodes]
								Time | Z[0] | Z[1] | Z[2] =====> Z[NoNodes]   */
      
  
	//Print to screen
	cout.flush();
	cout << "RoBoCoP: Time is " << setprecision(2) << fixed << Time << " years\r";

	//test if output file already exists
	int FileExists = 0;
	ifstream oftest(OutputFileName.c_str());
	if (oftest) FileExists = 1;
	oftest.close();

	//open the output filestream and write headers
	ofstream WriteCoastFile;
	if (FileExists == 0)
	{
		WriteCoastFile.open(OutputFileName.c_str());
		if (WriteCoastFile.is_open()) WriteCoastFile << "Header :)" << endl;
	}
	WriteCoastFile.close();

	//open output filestream again to  coastline data
	WriteCoastFile.open(OutputFileName.c_str(), fstream::app|fstream::out);

	//Check if file exists if not open a new one and write headers
	if (WriteCoastFile.is_open())
	{
		//write X
		WriteCoastFile << setprecision(4) << Time;
		for (int i=0; i<NoNodes; ++i) WriteCoastFile << setprecision(10) << " " << X[i];
		WriteCoastFile << endl;

		//write Y
		WriteCoastFile << setprecision(4) << Time;
		for (int i=0; i< NoNodes; ++i) WriteCoastFile << setprecision(10) << " " << Z[i];
		WriteCoastFile << endl;
	}
	
	else
	{
		//report errors
		cout << "RoBoCoP.WriteCoast: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
}

void RoBoCoP::WriteErosion(string OutputFileName, double Time)
{
  /* Writes a coastline object X and Z coordinates to file for a given time
		If the file already exists the data will be ed else a new file is
		created.

		File format is 	Time | X[0] | X[1] | X[2] =====> X[NoNodes]
								Time | Z[0] | Z[1] | Z[2] =====> Z[NoNodes]   */
      
  
	//Print to screen
	cout.flush();
	cout << "RoBoCoP: Time is " << setprecision(2) << fixed << Time << " years\r";

	//test if output file already exists
	int FileExists = 0;
	ifstream oftest(OutputFileName.c_str());
	if (oftest) FileExists = 1;
	oftest.close();

	//open the output filestream and write headers
	ofstream WriteCoastFile;
	if (FileExists == 0)
	{
		WriteCoastFile.open(OutputFileName.c_str());
		if (WriteCoastFile.is_open()) WriteCoastFile << "Header :)" << endl;
	}
	WriteCoastFile.close();

	//open output filestream again to  coastline data
	WriteCoastFile.open(OutputFileName.c_str(), fstream::app|fstream::out);

	//Check if file exists if not open a new one and write headers
	if (WriteCoastFile.is_open())
	{
		//write Erosion
		WriteCoastFile << setprecision(4) << Time;
		for (int i=0; i<NoNodes; ++i) WriteCoastFile << setprecision(10) << " " << Erosion[i];
		WriteCoastFile << endl;

		//write Y
		WriteCoastFile << setprecision(4) << Time;
		for (int i=0; i< NoNodes; ++i) WriteCoastFile << setprecision(10) << " " << Z[i];
		WriteCoastFile << endl;
	}
	
	else
	{
		//report errors
		cout << "RoBoCoP.WriteCoast: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
}
#endif
