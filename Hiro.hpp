/*------------------------------------------------------------------------

	Hiro.hpp
	
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

/** @file Hiro.hpp
@author Martin D. Hurst, University of Glasgow

@version Version 0.0.1
@brief Hiro object for simulating the evolution of rocky coastal profile
@details This object contains a simple geomorphic model for the evolution of a rocky
  coastal profile 
*/

/**
@mainpage
This is the documentation for the Hiro model
These pages describe the software.
@author Martin D. Hurst, University of Glasgow

------------------------------------------------------------------------*/

#ifndef Hiro_HPP
#define Hiro_HPP

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

using namespace std;

/*/////////////////////////////////////////////////////////////////////////////////////////
//TEMPLATES
/////////////////////////////////////////////////////////////////////////////////////////*/

//Dummy class for RockyCoastCRN
class RockyCoastCRN;

///@brief Main coastal platform object.
class Hiro
{
  friend class RockyCoastCRN;
  
	private:
	
	  /* MEMBER DECLARATIONS */
		
		// SPATIAL DOMAIN DECLARATIONS
		int NXNodes, NZNodes;               // Number of nodes across the coastline
		double dZ, dX;                      // Vertical and horizontal spacing of nodes
		
		vector<double> Z;			// elevation (m)
		vector<double> Xz;	   // cross shore distance at each elevation (m)
		
		vector <double> X;		// cross shore distance (m)
		vector <double> Zx;		// elevation at each cross shore distance(m)
		
		vector< vector<int> > MorphologyArray;		// array to store morphology
		vector< vector<double> > ResistanceArray;	// array to store resistance
		
		//Positions of min and max tide in X and Y
		int MinTideXInd, MaxTideXInd, MinTideYInd, MaxTideYInd;
		
		//
		double SurfZoneGradient;
		double SurfZoneWidth;
		
		// PROCSES DOMAIN DECLARTIONS
		vector<double> Bw_Erosion;		//Back wear erosion
		vector<double> Dw_Erosion;		//Down wear erosion
		vector<double> Weathering;		//Wearthering erosion
		
		// SEA LEVEL DECLARATIONS
		//vector<double> RSLTime;             //Times for relative sea level elevations
		//vector<double> RSLRate;             //Relative sea level elevations, will be length[1] if constant?
		double SeaLevelRise;
		double SeaLevel;
		int SeaLevelInd;
		
		// TIDES DECLARATIONS
		double TidalRange;              		//Tidal Range in metres
//		double TidalPeriod;                 //Tidal Period
//		vector<double> TideLevels;          //Vector of tide levels
//		vector<double> WaterLevels;         //vector containing water levels
//		vector<double> WaterDepths;         //vector containing water depths
		int NTideValues;                //number of tidal values (.size() of tidelevels vector)
		int WaterLevelYInd, WaterLevelXInd;
		
		// WAVE DECLARATIONS
		double MeanWavePeriod;
		double StdWavePeriod;
	  	double MeanWaveHeight;
  		double StdWaveHeight;
  		double WavePeriod;
		double WaveHeight;
		double BreakingWaveHeight;
		double BreakingWaveDist;
		double BreakingWaveWaterDepth;
		double WaveAttenuationConst;
		double BreakingPointX, BreakingPointY;
		int BreakingPointXInd, BreakingPointYInd;
		
		
		//TIME CONTROL PARAMETERS
		double Time, MaxTime, dt;
		
		//PHYSICAL CONSTANTS
		double rho_w;
		double g;
		
		//Constants and controlling parameters
		//These should be read from an input file?
		double SubmarineDecayConst, StandingWaveConst;
		double BreakingWaveConst, BrokenWaveConst;
		double BreakingWaveDecay, BrokenWaveDecay;
		double WeatheringConst;
		
		//Wave pressure parameters, check these with Hiro at some point
		double StandingWavePressure_Bw, BreakingWavePressure_Bw, BrokenWavePressure_Bw;
		double StandingWavePressure_Dw, BreakingWavePressure_Dw, BrokenWavePressure_Dw;
		
		//downwear decay const as a function of water depth
		//depends on wave height so defined inline
		double DepthDecay;		
		double CliffWeatheringRate;
		
		//This will need to be populated in the initialise tides function
		vector<double> WeatheringEfficacy;
		vector<double> ErosionShapeFunction;
		
		int PressureDistMinInd, PressureDistMaxInd;
		
		
		double NDV;    // No data value
	
    /* FUNCTION DECLARATIONS */

		//Initialise Functions
		void Initialise();
		void Initialise(double dZ, double dX);
		
	protected:

	public:
	
		/* PlatformCRN Initialisation functions */

		/// @brief Empty initialisation function for Hiro
		/// @author Martin D. Hurst 
		/// @date 27/02/2017
		Hiro()	
		{
			Initialise();
		}
		
		Hiro(double dZ, double dX)
		{
			Initialise(dZ, dX);
		}
		
    		
		//Initialise Tides
		void InitialiseTides(double TideRange);
		
		//Initialise Waves
		void InitialiseWaves(double WaveHeight_Mean, double WaveHeight_StD, double WavePeriod_Mean, double WavePeriod_StD);
		
		//Update Sea Level
		void UpdateSeaLevel(double SLRRate);
		
		//Sample a wave
		void GetWave();
		
		// Function to initialise weathering shape function 
		void InitialiseWeathering();
		
		/// @brief Launch the main program loop to evolve Hiro coast
		/// @details This function evolves a rocky coastal platform through time.
		///	@author Martin D. Hurst 
		/// @date 27/02/2017
		void EvolveCoast();
		
		void CalculateBackwearing();
		void CalculateDownwearing();
		
		void ErodeBackwearing();
		void ErodeDownwearing();
		void MassFailure();
		
		void IntertidalWeathering();
		void SupratidalWeathering();
		
		void UpdateMorphology();
		
		/// @brief Writes the platform morphology to file
		/// @details This function writes the elevations of the platform surface at the current time to
		///   a file. If the file exists, this is appended.
		///	@author Martin D. Hurst 
		/// @date 27/02/2017
		void WriteProfile(string OutputFileName, double Time);
		
		/// @brief Writes the ResistanceArray to file
		/// @details This function writes the rock Resistance at the current time to
		///   a file.
		/// @author Martin D. Hurst 
		/// @date 09/03/2017
		void WriteResistance(string OutputFileName, double Time);
		
		/// @brief Get X coordinates
		/// @return X coordinates
		///	@author Martin D. Hurst 
		/// @date 27/02/2017
		vector<double> get_X() const { return X; }
		
		/// @brief Get surface CRN concentration
		/// @return Surface CRN concentration
		///	@author Martin D. Hurst 
		/// @date 27/02/2017
		vector<double> get_Z() const { return Z; }
		
};

#endif
