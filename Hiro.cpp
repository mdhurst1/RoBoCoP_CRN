/*------------------------------------------------------------------------

	Hiro.cpp
	
	C++ implementation of Hiro Matsumoto's Shore Platform Model

	Matsumoto, H., Dickson, M. E., & Kench, P. S. (2016)
	An exploratory numerical model of rocky shore profile evolution. 
	Geomorphology, 268, 98-109. http://doi.org/10.1016/j.geomorph.2016.05.017
	
	Martin D. Hurst, University of Glasgow
	Hironori Matsumoto, University of Auckland
	
	February 2017
	
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

#ifndef Hiro_CPP
#define Hiro_CPP

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
#include "Hiro.hpp"

using namespace std;

void Hiro::Initialise()
{
  /* initialise an empty Hiro object as default */
  printf("\nHiro.Initialise: Initialised an empty Hiro object\n");
}

void Hiro::Initialise(double dZ_in, double dX_in)
{
	/* initialise a verictal cliff Hiro object */
	printf("\nHiro.Initialise: Initialised a Hiro as a vertical cliff\n");
	
	//Declare spatial stuff
	dZ = dZ_in;
	dX = dZ_in;
	NXNodes = round(500./dX);
	NZNodes = round(75./dZ);	

	//declare an array of zeros for assignment of Z vector
	Z = vector<double>(NZNodes,0.0);
	for (int i=0; i<NZNodes; ++i) Z[i] = dZ*(0.5*NZNodes-i);

	//declare arrays of zeros to initalise various other vectors
	vector<double> ZZeros(NZNodes,0);
	vector<double> XZeros(NXNodes,0);
	X = ZZeros;
	Bw_Erosion = ZZeros;
	Dw_Erosion = XZeros;
	Weathering = ZZeros;

	//declare an array of ones for assignment of the Morphology Array
	MorphologyArray = vector< vector<int> >(NZNodes,vector<int>(NXNodes,1));
	ResistanceArray = vector< vector<double> >(NZNodes,vector<double>(NXNodes,0));
		
	//default time interval
	dt = 1.;
	
	//Set sea level to zero to begin with
	SeaLevel = 0;	
}

//void Hiro::InitialiseTides(double TidalAmplitude, double TidalPeriod)
//{
//  /* intialise the tides as a simple cosine wave (i.e. simple diurnal/semidiurnal) */

//  //Tidal Timestep (in hours)
//  double dT = 0.01;    
//  double TideTime;
//  
//  //Setup tide water levels array
//  NTideValues = (int)(TidalPeriod/(2.*dT));
//  vector<double> Empty(NTideValues,-9999);
//  WaterLevels = Empty;
//  
//  //loop through and calculate water levels
//  for (int i=0, N=WaterLevels.size(); i<N; ++i)
//  {
//    TideTime = i*dT;
//    WaterLevels[i] = -TidalAmplitude*cos(2.*M_PI*TideTime/TidalPeriod);
//  }  
//}

void Hiro::InitialiseWaves(double WaveHeight_Mean, double WaveHeight_StD, double WavePeriod_Mean, double WavePeriod_StD)
{
  /* intialise waves as a single wave */

  MeanWavePeriod = WavePeriod_Mean;
  StdWavePeriod = WavePeriod_StD;
  MeanWaveHeight = WaveHeight_Mean;
  StdWaveHeight = WaveHeight_StD;
}

void Hiro::GetWave()
{
	//declare temp variables
	double OffshoreWaveHeight;
	double rand1, rand2;

	//Get two random numbers and generate wave data
	rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	WavePeriod = MeanWavePeriod + StdWavePeriod*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	OffshoreWaveHeight = MeanWaveHeight + StdWaveHeight*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));

	//Breaking Wave Height calculated following Komar and Gaughan (1972)
	BreakingWaveHeight = 0.39*pow(g,0.2)*pow(WavePeriod,0.4)*pow(OffshoreWaveHeight,2.4);
	BreakingWaveDist = BreakingWaveHeight/2.;

	//Water depth of breaking wave
	BreakingWaveWaterDepth = BreakingWaveHeight/0.78;
	
}

void Hiro::UpdateSeaLevel(double SLRRate)
{
	/*Update sea level based on a constant sea level rise rate*/
	SeaLevel += SLRRate*dt;
}

Hiro::Backwear()
{
	//Reset backwear vector
	vector<double> ZZeros(NZNodes,0);
	Bw_Erosion = ZZeros;
	
	//Loop across all intertidal elevations
	LowTideInd = SeaLeveli-0.5*TidalRange/dZ;
	HighTideInd = SeaLeveli+0.5*TidalRange/dZ;
	
	for (int i=LowTideInd; i<=HighTideInd; ++i)
	{
		//Estimate horizontal breaking point
		//Elevation of breaker point
		SurfZoneBottomZ = Z[i]-BreakingWaveWaterDepth;
		//BreakingPointInd = i-round(BreakingWaveWaterDepth*dZ);?
		
		bool SurfZone=false;
		int ii=i;
		while (SurfZone == false)
		{
			if (Z[ii] < SurfZoneBottomZ)
			{
				//Find Position X of Surfzone 
				BreakingPointX = Xz[ii];
				BreakingPointInd = ii;
				SurfZone = true;
			}
			--ii;
		}
		
		//Set Wave Type
		if (X[i] == 1) WaveType = 1;
		else if (X[i]-BreakingPointX<=0) WaveType = 1;
		else if ((X[i]-BreakingPointX)<BreakingWaveDist) WaveType = 2;
		else WaveType = 3;
		
		//Determine Surfzone Width
		//Get surf zone mean platform gradient
		//Set it to super steep if vertical.
		if (X[i-1] != X[BreakingPointInd+1]) SurfZoneGradient = abs((Zx[i-1]-Zx[BreakingPointInd+1])/(X[i-1]-X[BreakingPointInd+1]));
		else SurfZoneGradient = 100000;
		
		SurfZoneWidth = WaveHeight/SurfZoneGradient;
		
		//Set the wave attenuation constant #2
		WaveAttenuationConst = -log(SurfZoneWidth)/SurfZoneWidth;
		
		//Determine Backwear erosion
		//Standing wave
		if (WaveType == 1)
		{
			//Loop across pressure distribution function
			//This may have some problems!
			for (int ii=i-PressureDistZ1/dZ+1; ii>i+PressureDistZ2/dZ-1; ++ii)
			{
				// Calculate wave force and update backwear at each elevation
				WaveForce = StandingWaveConst*WaveHeight*ErosionShapeFunction[i]*UnbrokenWavePressure[ii];
				Bw_Erosion[ii] += WaveForce;
			}
		}
		//Breaking wave
		//For a breaking wave, first deal with backwear for the standing wave part, 
		// then the breaking part
		else if (WaveType == 2)
		{
			//Loop across pressure distribution function 
			//This may have some problems!
			for (int ii=i-PressureDistZ1/dZ+1; ii>i+PressureDistZ2/dZ-1; ++ii)
			{
				//need to add for condition where changes to broken wave above water level in pressure distribution function
				if (X[ii] < BreakingPointX) WaveForce = StandingWaveConst*WaveHeight*ErosionShapeFunction[i]*UnbrokenWavePressure[ii];
				else WaveForce = BrokenWaveConst*WaveHeight*BreakingWaveDecay*ErosionShapeFinction[i]*BreakingWavePressure[ii]*exp(-WaveAttenuationConst*(X[ii]-BreakingWaveDist));
				Bw_Erosion[ii] += WaveForce;
			}			
		}
		//Broken wave
		else if (WaveType == 3)
		{
			//Loop across pressure distribution function 
			//This may have some problems!
			for (int ii=i-PressureDistZ1/dZ+1; ii>i+PressureDistZ2/dZ-1; ++ii)
			{
				if (X[ii] < BreakingPointX) WaveForce = StandingWaveConst*WaveHeight*ErosionShapeFunction[i]*UnbrokenWavePressure[ii];
				else if (X[ii]<(BreakingPointX+BreakingWaveDist)) WaveForce = BreakingWaveConst*WaveHeight*BreakingWaveDecay*ErosionShapeFinction[i]*BrokenWavePressure[ii]*exp(-WaveAttenuationConst*(X[ii]-BreakingWaveDist));
				else force = BrokenWaveConst*WaveHeight*ErosionShapeFunction[i]*BrokenWavePressure*exp(-WaveAttenuationConst*(X[i]-(BreakingPointX+BreakingWaveDist)));
				Bw_Erosion[ii] += WaveForce;
			}
		}
	}
}

Hiro::Downwear()
{
	//Reset downwear vector
	vector<double> ZZeros(NZNodes,0);
	Dw_Erosion = ZZeros;
	
	//Find surficial cells in intertidal zone
	//Loop across all intertidal elevations
	LowTideInd = SeaLeveli-0.5*TidalRange/dZ;
	HighTideInd = SeaLeveli+0.5*TidalRange/dZ;
	
	for (int i=LowTideInd; i<=HighTideInd; ++i)
	{
		//Get wave function needs to calculate a bunch of stuff?
		GetWave();
		
		//Loop from water level down and determine force
		for (int ii=i; ii<Something; --ii)
		{
			//Standing Waves
			if (X[i]<BreakingPointX)
			{
				WaveForce = StandingWaveConst*WaveHeight*ErosionShapeFunction[i]*UnbrokenWavePressure;
				DepthDecay = -log(SubmarineDecayConst)/WaveHeight;
			}	
			//Breaking Waves
			else if (X[i]<(BreakingPointX+BreakingWaveDist))
			{
				WaveForce = BreakingWaveConst*WaveHeight*ErosionShapeFunction[i]*BreakingWavePressure*exp(-BreakingWaveDecay*(X[i]-BreakingPointX));
				DepthDecay = -log(SubmarineDecayConst)/(WaveHeight*exp(-BreakingWaveDecay*(X[i]-BreakingPointX)));
			}
			//Broken Waves
			else
			{
				WaveForce = BrokenWaveConst*WaveHeight*ErosionShapeFunction[i]*BrokenWavePressure*exp(-BrokenWaveDecay*(X[i]-(BreakingPointX+BreakingWaveDist));
				DepthDecay = -log(SubmarineDecayConst)/(WaveHeight*exp(-BrokenWaveDecay*(X[i]-(BreakingPointX+BreakingWaveDist))));
			}
			WaterDepth = Z[i]-Z[ii];
			Dw_Erosion += WaveForce*exp(-DepthDecay*WaterDepth);
		}
	}
}

Hiro::Weathering()
{
	//Reset weathering vector
	vector<double> ZZeros(NZNodes,0);
	Weathering = ZZeros;
	
	//Loop across the tidal range
	
	//Find surficial cells in intertidal zone
	//Loop across all intertidal elevations
	LowTideInd = SeaLeveli-0.5*TidalRange/dZ;
	HighTideInd = SeaLeveli+0.5*TidalRange/dZ;
	
	for (int i=LowTideInd; i<=HighTideInd; ++i)
	{
		//Calculate Weathering
		WeatheringForce = WeatheringConst*WeatheringEfficacy[i];
		
		//How are we going to get j? i.e. x-position in the array?
		//Need a loop in here moving from bottom to top of tidal range in x-position
		RemainingResistance = ResistanceArray[i][j];
		ResistanceArray[i][j] -= WeatheringForce; 
		
		//If resistance is less than zero then cell is lost to weathering erosion
		//excess weathering force is applied to the block behind
		if (ResistanceArray[i][j] < 0)
		{
			MorphologyArray[i][j] = 0;
			ResistanceArray[i][j] = 0;
			WeatheringForce -= RemainingResistance;
		}
	
	}
	ResistanceArray;
	
}
//void Hiro::EvolveCoast()
//{
//  /* Function to evolve the coastal profile through time following
//      Trenhaile (2000)  */
//  
//  //declarations
//  double SurfZoneTopZ, SurfZoneBottomZ, SurfZoneTopX, SurfZoneBottomX, SurfZoneWidth, SurfForce, WaterDepth;
//  
//  //vector to sum erosion
//  vector<double> Empty(NoNodes,0);
//  Erosion = Empty;
//  
//  //Determine surf zone width
//  for (int t=0, T=WaterLevels.size(); t<T; ++t)
//  {
//    //Elevation of surf zone limits
//    SurfZoneBottomZ = SeaLevel+WaterLevels[t]-BreakingWaveWaterDepth;
//    SurfZoneTopZ = SeaLevel+WaterLevels[t];
//    
//    //find horizontal extent of surfzone
//    bool SurfZone=false;
//    int i=NoNodes-1;
//    while (SurfZone == false)
//    {
//      if (Z[i] > SurfZoneBottomZ)
//      {
//        //Find Position X of Surfzone bottom by interpolating
//        SurfZoneBottomX = X[i] + (X[i+1]-X[i])*(Z[i]-SurfZoneBottomZ)/(Z[i+1]-Z[i]);
//        SurfZone = true;
//      }
//      --i;
//    }
//    while (SurfZone == true)
//    {
//      if (Z[i] > SurfZoneTopZ)
//      {
//        //Find Position X of Surfzone bottom by interpolating
//        SurfZoneTopX = X[i] + (X[i+1]-X[i])*(Z[i]-SurfZoneTopZ)/(Z[i+1]-Z[i]);
//        SurfZone = false;
//        ++i;
//      }
//      --i;
//    }
//    
//    //Calculate Surf Zone width
//    SurfZoneWidth = SurfZoneBottomX-SurfZoneTopX;
//    
//    //Calculate Surf Force reaching the water line
//    SurfForce = 0.5*rho_w*(BreakingWaveWaterDepth)*exp(-k*SurfZoneWidth);
//    
//    //Loop offshore and calculate erosion
//    while (i < NoNodes)
//    {
//      WaterDepth = SeaLevel+WaterLevels[t]-Z[i];
//      if (fabs(WaterDepth)<0.0001) WaterDepth=0;
//      else if (WaterDepth < 0) WaterDepth = 100.;

//      Erosion[i] += M*SurfForce*exp(-0.5*WaterDepth);
//      ++i;
//    }
//  }
//  //Do the erosion/update the profile
//  double XCliff = X[NoNodes-1];
//  for (int i=NoNodes-1; i>-1; --i)
//  {
//    X[i] -= Erosion[i]*dt;
//    
//    //check for overhangs and remove
//    if (X[i] > XCliff) X[i] = XCliff;
//    XCliff = X[i];
//  }
//}

void Hiro::WriteProfile(string OutputFileName, double Time)
{
  /* Writes a Hiro object X coordinates to file, each value spans dZ in elevation
		File format is 	
		
		StartZ dZ
		Time | X[0] | X[1] | X[2] =====> X[NoNodes] */
      
  
	//Print to screen
	cout.flush();
	cout << "Hiro: Writing output at Time " << setprecision(0) << fixed << Time << " years\r";

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
		if (WriteCoastFile.is_open()) WriteCoastFile << Z[0] << " " << dZ << endl;
	}
	WriteCoastFile.close();

	//open output filestream again to  coastline data
	WriteCoastFile.open(OutputFileName.c_str(), fstream::app|fstream::out);

	//Check if file exists if not open a new one and write headers
	if (WriteCoastFile.is_open())
	{
		//write X
		WriteCoastFile << setprecision(4) << Time;
		for (int i=0; i<NZNodes; ++i) WriteCoastFile << setprecision(10) << " " << X[i];
		WriteCoastFile << endl;
	}
	else
	{
		//report errors
		cout << "Hiro.WriteCoast: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
}
#endif
