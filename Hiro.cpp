/*------------------------------------------------------------------------

	Hiro.cpp
	
	C++ implementation of Hiro Matsumoto's Shore Platform Model

	Matsumoto, H., Dickson, M. E., & Kench, P. S. (2016)
	An exploratory numerical model of rocky shore profile evolution. 
	Geomorphology, 268, 98-109. http://doi.org/10.1016/j.geomorph.2016.05.017
	
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
	
	//PHYSICAL CONSTANTS
	rho_w = 1025.;
	g = 9.81;

	//Define these in an input parameter file?
	SubmarineDecayConst = 0.1;
	StandingWaveConst = 0.01;
	BreakingWaveConst = 10.;
	BrokenWaveConst = 0.1;
	BreakingWaveDecay = 0.1;
	BrokenWaveDecay = 0.01;
	WeatheringConst = 0.005;
	RockResistance = 0.5;
	
	//Wave pressure parameters, check these with Hiro at some point
	StandingWavePressure_Bw = 0.1;
	BreakingWavePressure_Bw = 0.1;
	BrokenWavePressure_Bw = 0.1;
	StandingWavePressure_Dw = 0.1;
	BreakingWavePressure_Dw = 0.1;
	BrokenWavePressure_Dw = 0.1;
	
	//Cliff control params
	CliffHeight = 10.;
	CliffFailureDepth = 1.;
	
	//Declare spatial stuff
	dZ = dZ_in;
	dX = dX_in;
	NXNodes = 1000;
	NZNodes = round(2.*CliffHeight/dZ)+1;	

	//declare an array of zeros for assignment of Z vector
	Z = vector<double>(NZNodes,0.0);
	for (int i=0; i<NZNodes; ++i) Z[i] = dZ*(0.5*(NZNodes-1)-i);

	//declare arrays of zeros to initalise various other vectors
	vector<double> ZZeros(NZNodes,0);
	vector<double> XZeros(NXNodes,0);
	X = XZeros;
	for (int j=0; j<NXNodes; ++j) X[j] = dX*j;
	Zx = XZeros;
	Xz = ZZeros;
	
	Bw_Erosion = ZZeros;
	Dw_Erosion = XZeros;
	Weathering = ZZeros;
	
	//Indices trackers
	XInd = vector<int>(NZNodes,0);
	ZInd = vector<int>(NXNodes,0);
	
	//declare an array of ones for assignment of the Morphology Array
	MorphologyArray = vector< vector<int> >(NZNodes,vector<int>(NXNodes,1));
	ResistanceArray = vector< vector<double> >(NZNodes,vector<double>(NXNodes,RockResistance));
		
	//default time interval
	Time = 0;
	TimeInterval = 1;
	EndTime = 10000;
	
	//default print conditions to print every timestep
	PrintTime = 0;
	PrintInterval = TimeInterval;
	
	//Set sea level to zero to begin with, and the ind, this will get updated later
	SeaLevelRise = 0;
	SeaLevel = 0;
	SeaLevelInd = 0;	
}

void Hiro::InitialiseTides(double TideRange)
{
	//setup tides
	
	//declare temporary variables
	double Total = 0;
	double Max = 0;
	
	TidalRange = TideRange;
	UpdateMorphology();
	
	// Make erosion shape function based on tidal duration
	NTideValues = (int)(TidalRange/dZ)+1;
	vector<double> EmptyTideVec(NTideValues,TidalRange/2.);	
	ErosionShapeFunction = EmptyTideVec;
	
	// Loop over tidal range and assign weights
	// If narrow tidal range just use two separate sin functions
	// If wide use two overlapping tidal functions
	if (NTideValues < 20)
	{
		for (int i=0; i<NTideValues; ++i)
		{
			ErosionShapeFunction[i] = sin(i*dZ*M_PI/(0.5*TidalRange));
			if (i == (NTideValues-1)/2) ErosionShapeFunction[i] += ErosionShapeFunction[i-1];
			//else 
			if (i*dZ>0.5*TidalRange) ErosionShapeFunction[i] *= -1;
			Total += ErosionShapeFunction[i];
			if (ErosionShapeFunction[i] > Max) Max = ErosionShapeFunction[i];
		}
	}
	else
	{
		for (int i=0; i<0.55*NTideValues; ++i)
		{
			ErosionShapeFunction[i] = sin(i*dZ*M_PI/(0.55*TidalRange));
			Total += ErosionShapeFunction[i];
			if (ErosionShapeFunction[i] > Max) Max = ErosionShapeFunction[i];
		}
		for (int i=0.45*NTideValues; i<NTideValues; ++i)
		{
			ErosionShapeFunction[i] += sin((i*dZ-0.45*TidalRange)*M_PI/(0.55*TidalRange));
			Total += ErosionShapeFunction[i];
			if (ErosionShapeFunction[i] > Max) Max = ErosionShapeFunction[i];
		}
	}
	
	//Normalise values to total
	//Or should I be normalising to Max Value!?
	for (int i=0; i<NTideValues; ++i) ErosionShapeFunction[i] /= Total;
	//for (int i=0; i<NTideValues; ++i) ErosionShapeFunction[i] /= Max;
	
	//Initialise weathering shape function
	InitialiseWeathering();	
}

void Hiro::InitialiseGeology(double CliffHeightNew, double CliffFailureDepthNew, double RockResistanceNew, double WeatheringConstNew)
{
	/* Function to set the cliff height, failure depth, rock resistance and 
		weathering rate constant */
		
	CliffHeight = CliffHeightNew;
	CliffFailureDepth = CliffFailureDepthNew;
	RockResistance = RockResistanceNew;
	WeatheringConst = WeatheringConstNew;
	
	//Loop across the resistance array and reset all values
	for (int i=0;i<NZNodes; ++i)
	{
		for (int j=0; j<NXNodes; ++j)
		{
			ResistanceArray[i][j] = RockResistance;
		}
	}
}

void Hiro::InitialiseWeathering()
{
	/* Weathering Efficacy Function guided by Trenhaile and Kanayay (2005)
	using a log-normal distribution across the tidal range */
	
	// declare control parameters for distribution
	//double sigma = 0.5;
	//double Theta = 0;
	double MaxEfficacy = 0;
	
	// This m value is tailored to cause a distribution peak at 3/4 of the tidal range
	// as in Matsumoto et al. (2016)
	//double m = 1.1665; 
	
	// Make weathering shape function based on tidal duration
	NTideValues = (int)(TidalRange/dZ)+1;
	vector<double> EmptyTideVec(NTideValues,0);	
	WeatheringEfficacy = EmptyTideVec;
	
	//Create a vector ranging from 0 to 10 with NTideValues
	vector<double> LogNormalDistX(NTideValues,0);
	for (int i=0; i<NTideValues; ++i) LogNormalDistX[i] = 10.*i/(NTideValues-1);

	//Normalise so that Max value is 1
	//MaxEfficacy = exp(-(pow(log(2.5-Theta)-m,2.)/(2.*pow(sigma,2.)))) / ((2.5-Theta)*sigma*sqrt(2.*M_PI)); 
	
	//Create log normal weathering efficacy shape function	
	//for (int i=1; i<NTideValues ;++i) 
	//{
	//	WeatheringEfficacy[i] = (exp(-((pow(log(LogNormalDistX[i]-Theta)-m,2.))/(2*pow(sigma,2.)))) / ((LogNormalDistX[i]-Theta)*sigma*sqrt(2.*M_PI)));
	//	if (WeatheringEfficacy[i] > MaxEfficacy) MaxEfficacy = WeatheringEfficacy[i];
	//}
	
	//Create log normal weathering efficacy shape function	
	for (int i=1; i<NTideValues ;++i) 
	{
		if (i<((NTideValues-1)/4)) WeatheringEfficacy[i] = exp(-(pow(i-(NTideValues/4.),2.)/(NTideValues/2.)));
      else WeatheringEfficacy[i] = exp(-(pow(i-(NTideValues/4.),2.))/(NTideValues*NTideValues/10.));
		if (WeatheringEfficacy[i] > MaxEfficacy) MaxEfficacy = WeatheringEfficacy[i];
	}
	
	for (int i=1; i<NTideValues ;++i) WeatheringEfficacy[i] /= 	MaxEfficacy;
}

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
	//double OffshoreWaveHeight;
	//double rand1, rand2;

	//Get two random numbers and generate wave data
	//rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	//WavePeriod = MeanWavePeriod + StdWavePeriod*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	//rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	//OffshoreWaveHeight = MeanWaveHeight + StdWaveHeight*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));

	//Breaking Wave Height calculated following Komar and Gaughan (1972)
	//BreakingWaveHeight = 0.39*pow(g,0.2)*pow(WavePeriod,0.4)*pow(OffshoreWaveHeight,2.4);
	WaveHeight = MeanWaveHeight;
	BreakingWaveHeight = MeanWaveHeight;
	BreakingWaveDist = BreakingWaveHeight/2.;

	//Water depth of breaking wave
	BreakingWaveWaterDepth = BreakingWaveHeight/0.78;
	
	//Get Wave Pressure Distribution as Uniform
	PressureDistMaxInd = 0.5*BreakingWaveHeight/dZ; 
	PressureDistMinInd = -0.5*BreakingWaveHeight/dZ;
}

void Hiro::InitialiseSeaLevel(double SLR)
{
	SeaLevelRise = SLR;
}

void Hiro::UpdateSeaLevel()
{
	/*Update sea level based on a constant sea level rise rate*/
	SeaLevel += SeaLevelRise*TimeInterval;
}

void Hiro::CalculateBackwearing()
{
	//Declare temporary variables
	double WaveForce, SurfZoneBottomZ; //, SurfZoneBottomX;
	int WaveType;
	int BreakingPointZInd = 0;
	
	//Reset backwear vector
	vector<double> ZZeros(NZNodes,0);
	Bw_Erosion = ZZeros;
	
	//Loop across all intertidal elevations
	for (int i=MaxTideZInd+1; i<MinTideZInd; ++i)
	{
		//Estimate horizontal breaking point
		//Elevation of breaker point
		SurfZoneBottomZ = Z[i]-BreakingWaveWaterDepth;
		for (int ii=i; ii<NZNodes; ++ii) 
		{
			if (Z[ii] < SurfZoneBottomZ)
			{
				BreakingPointZInd = ii;
				break;
			}
		}
		
		//If waves are breaking at the seaward edge, find the top of the seaward edge
		for (int ii=BreakingPointZInd; ii>i; --ii) 
		{
			if (Xz[ii] == Xz[BreakingPointZInd]) BreakingPointZInd = ii;
			else break;
		}

		//Set Wave Type
		if (Xz[i] == 0) WaveType = 1;
		else if (Xz[i]-Xz[BreakingPointZInd]<=0) WaveType = 1;
		else if ((Xz[i]-Xz[BreakingPointZInd])<BreakingWaveDist) WaveType = 2;
		else WaveType = 3;
		
		//Find the location where the broken wave starts
		//int BrokenWaveXInd = BreakingPointXInd + round(BreakingWaveDist/dX);
		
		//Determine Surfzone Width
		//Get surf zone mean platform gradient
		//This is critical!
//		double Slope;
//		double SumSlopes = 0;
//		int NSlopes = 0;
		if ((Xz[MaxXZInd] != Xz[BreakingPointZInd]))
		{
			SurfZoneGradient = abs((Z[MaxXZInd]-Z[BreakingPointZInd])/(Xz[MaxXZInd]-Xz[BreakingPointZInd]));
//			for (int jj=BrokenWaveXInd-1; X[jj]>Xz[i]; --jj)
//			{
//				if ((X[jj] != X[jj+1]) && (Zx[jj+1] < 0.5*CliffHeight))	
//				{
//					Slope = fabs((Zx[jj+1]-Zx[jj])/(X[jj+1]-X[jj]));
//					SumSlopes += Slope;
//					NSlopes++;
//				}
//			}
//			SurfZoneGradient = SumSlopes/NSlopes;
		}
		else 
		{
			SurfZoneGradient = 1.;
			//cout << endl << "setting SurfZoneGradient = 1" << endl;
		}
		
		//Limt SurfZoneGradient to 45 degrees!
		if (SurfZoneGradient > 1.) SurfZoneGradient = 1.;
		
		SurfZoneWidth = WaveHeight/SurfZoneGradient;
		
		//Set the wave attenuation constant #2
		BreakingWaveAttenuation = -log(BreakingWaveDecay)/BreakingWaveDist;
		BrokenWaveAttenuation = -log(BrokenWaveDecay)/SurfZoneWidth;
		
		//Determine Backwear erosion
		//Standing wave
		if (WaveType == 1)
		{
			//Loop across pressure distribution function, currently a const
			for (int ii=i+PressureDistMinInd; ii<=i+PressureDistMaxInd; ++ii)
			{
				// Calculate wave force and update backwear at each elevation
				WaveForce = StandingWaveConst*WaveHeight*ErosionShapeFunction[i-MaxTideZInd]*StandingWavePressure_Bw;
				Bw_Erosion[ii] += WaveForce;
			}
		}
		//Breaking wave
		//For a breaking wave, first deal with backwear for the standing wave part, 
		// then the breaking part
		else if (WaveType == 2 || WaveType == 3)
		{
			//Loop across pressure distribution function 
			//This may have some problems!
			for (int ii=i+PressureDistMinInd; ii<=i+PressureDistMaxInd; ++ii)
			{
				//need to add for condition where changes to broken wave above water level in pressure distribution function
				if (Xz[ii] < Xz[BreakingPointZInd]) 
				{
					WaveForce = StandingWaveConst*WaveHeight*ErosionShapeFunction[i-MaxTideZInd]*StandingWavePressure_Bw;
				}
				else if (Xz[ii] <= (Xz[BreakingPointZInd]+BreakingWaveDist))
				{
					BreakingWaveHeight = WaveHeight*exp(-BreakingWaveAttenuation*(Xz[ii]-Xz[BreakingPointZInd]));
					WaveForce = BreakingWaveConst*BreakingWaveHeight*ErosionShapeFunction[i-MaxTideZInd]*BreakingWavePressure_Bw;
				}
				else
				{
					BrokenWaveHeight = BreakingWaveDecay*WaveHeight*exp(-BrokenWaveAttenuation*(Xz[ii]-(Xz[BreakingPointZInd]+BreakingWaveDist)));
					WaveForce = BrokenWaveConst*BrokenWaveHeight*ErosionShapeFunction[i-MaxTideZInd]*BrokenWavePressure_Bw;
				}
				Bw_Erosion[ii] += WaveForce;
			}			
		}
	}
}

void Hiro::CalculateDownwearing()
{
	//Declare temporary variables
	double WaveForce, WaterDepth;
	
	//Reset downwear vector
	vector<double> ZZeros(NZNodes,0);
	Dw_Erosion = ZZeros;
	
	for (int i=MinTideZInd; i<=MaxTideZInd; ++i)
	{
		//Get wave function needs to calculate a bunch of stuff?
		GetWave();
		
		//Loop from water level down and determine force
		for (int ii=WaterLevelYInd; ii>0; --ii)
		{
			//Standing Waves
			if (X[i]<BreakingPointX)
			{
				WaveForce = StandingWaveConst*WaveHeight*ErosionShapeFunction[i-MaxTideZInd]*StandingWavePressure_Dw;
				DepthDecay = -log(SubmarineDecayConst)/WaveHeight;
			}	
			//Breaking Waves
			else if (X[i]<(BreakingPointX+BreakingWaveDist))
			{
				WaveForce = BreakingWaveConst*WaveHeight*ErosionShapeFunction[i-MaxTideZInd]*BreakingWavePressure_Dw*exp(-BreakingWaveDecay*(X[i]-BreakingPointX));
				DepthDecay = -log(SubmarineDecayConst)/(WaveHeight*exp(-BreakingWaveDecay*(X[i]-BreakingPointX)));
			}
			//Broken Waves
			else
			{
				WaveForce = BrokenWaveConst*WaveHeight*ErosionShapeFunction[i-MaxTideZInd]*BrokenWavePressure_Dw*exp(-BrokenWaveDecay*(X[i]-(BreakingPointX+BreakingWaveDist)));
				DepthDecay = -log(SubmarineDecayConst)/(WaveHeight*exp(-BrokenWaveDecay*(X[i]-(BreakingPointX+BreakingWaveDist))));
			}
			WaterDepth = Z[i]-Z[ii];
			Dw_Erosion[i] += WaveForce*exp(-DepthDecay*WaterDepth);
		}
	}
}

void Hiro::IntertidalWeathering()
{
	//Declare temporay variables
	double RemainingResistance, WeatheringForce;
	
	//Reset weathering vector
	vector<double> ZZeros(NZNodes,0);
	Weathering = ZZeros;
	
	//Loop across the tidal range
	for (int i=MaxTideZInd+1; i<MinTideZInd; ++i)
	{
		//Calculate Weathering
		WeatheringForce = WeatheringConst*WeatheringEfficacy[i-MaxTideZInd];
		
		//How are we going to get j? i.e. x-position in the array?
		//Need a loop in here moving from bottom to top of tidal range in x-position
		for (int j=0; j<=MaxXXInd; ++j)
		{
			//Check we're at a a surface cell
			if ((MorphologyArray[i][j] == 1) && (j == 0))
			{
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
			else if ((MorphologyArray[i][j] == 1) && ((MorphologyArray[i-1][j] == 0) || (MorphologyArray[i][j-1] == 0)))
			{
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
		}
	}
}


void Hiro::ErodeBackwearing()
{
	//Loop over all wet cells
	int j=0;
	double RemainingForce;
	
	for (int i=MinTideZInd+PressureDistMaxInd; i>=MaxTideZInd+PressureDistMinInd; --i)
	{
		//Find j ind somehow
		//loop across the active shoreface
		j=0;
		while (j < MaxXXInd)
		{
			if (MorphologyArray[i][j] == 1) break;
			++j;
		}
		
		//Check Backwear Force vs Resistance
		RemainingForce = Bw_Erosion[i];
		
		while (RemainingForce > ResistanceArray[i][j])
		{
			//Update remaining force
			RemainingForce -= ResistanceArray[i][j];
			 
			//For now assume that only one block can be removed at a time
			//Hiro has code that allows multiple blocks to be removed
			MorphologyArray[i][j] = 0;
			ResistanceArray[i][j] = 0;
			
			//iterate landward
			++j;
			if (j > MaxXXInd) MaxXXInd = j;
		
			//May also need soemthing to move ix_max landward by 1
			//If we run out of land dynamically add more columns
			if (j == NXNodes-2) 
			{
				//Grow them
				X.push_back(NXNodes*dX);
				Zx.push_back(Zx[NXNodes-1]);
				ZInd.push_back(ZInd[NXNodes-1]);
				for (int i=0; i<NZNodes; ++i)
				{
					MorphologyArray[i].push_back(MorphologyArray[i][NXNodes-1]);
					ResistanceArray[i].push_back(ResistanceArray[i][NXNodes-1]);
				}
				++NXNodes;
			}
		}
	}
}

void Hiro::ErodeDownwearing()
{
	//is there any reason why this needs to be a separate function?
	//Loop over all cells that get wet
	for (int j=0; j<MaxTideXInd; ++j)
	{
		// loop over the tidal range?
		for (int i=MaxTideZInd; i<=MinTideZInd; ++i)
		{
			// Check Downwear Force vs Resistance
			if (Dw_Erosion[j] >= ResistanceArray[i][j])
			{
				//For now assume that only one block can be removed at a time
				//Hiro has code that allows multiple blocks to be removed
				//I doubt this happens very often with downwear
				MorphologyArray[i][j] = 0;
				ResistanceArray[i][j] = 0;
			
				//Hiro then has some code to count the number of blocks removed
				//but not clear why this is needed
			}
		}
	}
}

void Hiro::SupratidalWeathering()
{
	//add this later
}

void Hiro::UpdateMorphology()
{
	//function to update morphology vectors and indices
	 
	// Find Sea Level in vertical
	// Only need to do this once if sea level isnt changing
	if ((SeaLevelInd == 0) || (SeaLevelRise != 0))
	{
		for (int i=0; i<NZNodes; ++i)
		{
			if (Z[i] < SeaLevel)
			{
				SeaLevelInd = i-1;
				break;
			}
		}
	}
	
	//Loop across all intertidal elevations
	MinTideZInd = SeaLevelInd+0.5*TidalRange/dZ;
	MaxTideZInd = SeaLevelInd-0.5*TidalRange/dZ;
	
	//Populate vector of X values in Z 
	MaxXXInd = 0;
	MaxXZInd = 0;
	for (int i=0; i<NZNodes; ++i)
	{
		for (int j=XInd[i]; j<NXNodes; ++j)
		{
			if (MorphologyArray[i][j] == 1) 
			{
				Xz[i] = X[j];
				XInd[i] = j;
				if (XInd[i] > MaxXXInd) 
				{
					MaxXXInd = XInd[i];
					MaxXZInd = i;
				}
				break;
			}
		}
	}
	
	//Grow the X direction arrays dynamically
	if (MaxXXInd > NXNodes-10)
	{
		//Grow them
		X.push_back(NXNodes*dX);
		Zx.push_back(Zx[NXNodes-1]);
		ZInd.push_back(ZInd[NXNodes-1]);
		for (int i=0; i<NZNodes; ++i)
		{
			MorphologyArray[i].push_back(MorphologyArray[i][NXNodes-1]);
			ResistanceArray[i].push_back(ResistanceArray[i][NXNodes-1]);
		}
		++NXNodes;
	}
	
	//Populate vector of Z values in X 
	for (int j=0; j<NXNodes; ++j)
	{
		for (int i=ZInd[j]; i<NZNodes; ++i)
		{
			if (MorphologyArray[i][j] == 1)
			{
				Zx[j] = Z[i]; 
				ZInd[j] = i;
				break;
			}
		}
	}
	
	//Determine intertidal range indices in X-direction
	bool LowTideFlag = false;
	bool HighTideFlag = false;
	for (int j=0; j<NXNodes; ++j)
	{
		if ((MorphologyArray[MaxTideZInd][j] == 1) && (HighTideFlag == false))
		{
			MaxTideXInd = j;
			HighTideFlag = true;
		}
		if ((MorphologyArray[MinTideZInd][j] == 1) && (LowTideFlag == false))
		{
			 MinTideXInd = j;
			 LowTideFlag = true;
		}
		if ((HighTideFlag == 1) && (LowTideFlag == 1)) break;
	}
}

void Hiro::MassFailure()
{
	//simple implementation for now, talk to Hiro about this
	//Cliff position taken from Highest elevation.
	double XCliff = Xz[0];
	
	//Find X position of notch
	double XMax = 0;
	int XMaxZInd = 0;
	
	for (int i=0;i<NZNodes; ++i)
	{
		if (Xz[i] > XMax)
		{
			XMax = Xz[i];
			XMaxZInd = i;
		}
	}
	
	//if big enough then delete
	if ((XMax-XCliff) > CliffFailureDepth)
	{
		for (int i=0; i<=XMaxZInd; ++i)
		{
			for (int j=0; j<XMax/dX; ++j)
			{
				MorphologyArray[i][j] = 0;
				ResistanceArray[i][j] = 0;
			}
		}
	}
}


void Hiro::EvolveCoast()
{
	/* 
		Function to evolve the coastal profile through time following
		Matsumoto et al. 2016 
     
		This is the main model loop. Can be executed here or from a driver file
  
		Martin Hurst 30/3/2017 
	*/
  
 	//Loop through time
	while (Time <= EndTime)
	{
		//Update Sea Level
		UpdateSeaLevel();

		//Get the wave conditions
		GetWave();

		//Calculate forces acting on the platform
		CalculateBackwearing();
		CalculateDownwearing();

		//Do erosion
		ErodeBackwearing();
		ErodeDownwearing();

		//Update the Morphology 
		UpdateMorphology();	  
		
		//Implement Weathering
		IntertidalWeathering();
		
		//Update the Morphology 
		UpdateMorphology();

		//Check for Mass Failure
		MassFailure();
		
		//Update the Morphology 
		UpdateMorphology();
		
		//print?
		if (Time >= PrintTime)
		{
			WriteProfile(OutputFileName, Time);
			PrintTime += PrintInterval;
		}

		//update time
		Time += TimeInterval;
	}
  
}

void Hiro::WriteProfile(string OutputFileName, double Time)
{
  /* Writes a Hiro object X coordinates to file, each value spans dZ in elevation
		File format is 	
		
		StartZ dZ
		Time | X[0] | X[1] | X[2] =====> X[NoNodes] */
      
  
	//Print to screen
	cout.flush();
	cout << "Hiro: Writing output at Time " << setprecision(2) << fixed << Time << " years\r";

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
		for (int i=0; i<NZNodes; ++i) WriteCoastFile << setprecision(10) << " " << Xz[i];
		WriteCoastFile << endl;
	}
	else
	{
		//report errors
		cout << "Hiro.WriteCoast: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
	WriteCoastFile.close();
}

void  Hiro::WriteResistanceArray(string OutputFileName, double Time)
{
  /* Writes a Hiro object Resistance matrix coordinates to file
		File format is 	
		
		Time	
			X[0][0]     |    X[0][1]    |   X[0][2]     =====>    X[0][NXNodes]
			X[1][0]     |    X[1][1]    |   X[1][2]     =====>    X[0][NXNodes]
			X[2][0]     |    X[2][1]    |   X[2][2]     =====>    X[0][NXNodes]
		      ||               ||             ||                      ||
		      \/               \/             \/                      \/
		X[NZNodes][0]  | X[NZNodes][1] | X[NZNodes][2] =====> X[NZNodes][NXNodes] */
      
  	//open the output filestream and write headers
	ofstream WriteFile;
	WriteFile.open(OutputFileName.c_str());
	WriteFile << Time << " " << dZ << " " << dX << endl;
	
	//Check if file exists if not open a new one and write headers
	if (WriteFile.is_open())
	{
		//write Resistance
		for (int i=0; i<NZNodes; ++i)
		{
			for (int j=0;j<NXNodes; ++j)
			{
				WriteFile << setprecision(5) << ResistanceArray[i][j] << " ";
			}
			WriteFile << endl;
		}
	}
	else
	{
		//report errors
		cout << "Hiro.WriteResistance: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
}

void  Hiro::WriteMorphologyArray(string OutputFileName, double Time)
{
  /* Writes a Hiro object Resistance matrix coordinates to file
		File format is 	
		
		Time	
			X[0][0]     |    X[0][1]    |   X[0][2]     =====>    X[0][NXNodes]
			X[1][0]     |    X[1][1]    |   X[1][2]     =====>    X[0][NXNodes]
			X[2][0]     |    X[2][1]    |   X[2][2]     =====>    X[0][NXNodes]
		      ||               ||             ||                      ||
		      \/               \/             \/                      \/
		X[NZNodes][0]  | X[NZNodes][1] | X[NZNodes][2] =====> X[NZNodes][NXNodes] */
      
  	//open the output filestream and write headers
	ofstream WriteFile;
	WriteFile.open(OutputFileName.c_str());
	WriteFile << Time << " " << dZ << " " << dX << endl;
	
	//Check if file exists if not open a new one and write headers
	if (WriteFile.is_open())
	{
		//write Resistance
		for (int i=0; i<NZNodes; ++i)
		{
			for (int j=0;j<NXNodes; ++j)
			{
				WriteFile << setprecision(5) << MorphologyArray[i][j] << " ";
			}
			WriteFile << endl;
		}
	}
	else
	{
		//report errors
		cout << "Hiro.WriteMorphology: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
}


#endif
