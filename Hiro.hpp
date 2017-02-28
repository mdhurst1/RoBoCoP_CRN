/*------------------------------------------------------------------------

	Hiro.hpp
	
	C++ implementation of Hiro Matsumoto's Shore Platform Model

	Matsumoto, H., Dickson, M. E., & Kench, P. S. (2016)
	An exploratory numerical model of rocky shore profile evolution. 
	Geomorphology, 268, 98-109. http://doi.org/10.1016/j.geomorph.2016.05.017
	
	Martin D. Hurst, University of Glasgow
	Hironori Matsumoto, University of Auckland

	Copyright (C) 2016, Martin Hurst
	
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
		double dZ;                          // Vertical spacing of nodes
		
		vector<double> X;		          		// cross shore distance (m)
		vector<double> Z;							// elevation (m)
		vector< vector<int> > MorphologyArray;		// array to store morphology
		vector< vector<double> > ResistanceArray;		// array to store resistance
		
		// SEA LEVEL DECLARATIONS
		vector<double> RSLTime;             //Times for relative sea level elevations
		vector<double> RSLRate;             //Relative sea level elevations, will be length[1] if constant?
		double SeaLevel;
		
		// TIDES DECLARATIONS
		double TidalPeriod;                 //Tidal Period
		double TidalRange;              //Tidal Amplitude
		vector<double> TideLevels;          //Vector of tide levels
		vector<double> WaterLevels;         //vector containing water levels
		vector<double> WaterDepths;         //vector containing water depths
		int NTideValues;                //number of tidal values (.size() of tidelevels vector)
		
		// WAVE DECLARATIONS
		double WaveHeight;
		double BreakingWaveHeight;
		double BreakingWaveWaterDepth;
		
		//TIME CONTROL PARAMETERS
		double Time, MaxTime, dt;
		
		//PHYSICAL CONSTANTS
		static const double rho_w = 1025.;
		static const double g = 9.81;
		static const double k = 0.02;
		static const double M = 0.0001;
		
		double WeatheringRate;
		double CliffWeatheringRate;
		
		static const double NDV = -9999;    // No data value
	
    /* FUNCTION DECLARATIONS */

		//Initialise Functions
		void Initialise();
		void Initialise(double dZ);
		
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
		
		Hiro(double dZ)
		{
			Initialise(dZ);
		}
		
    		
		//Initialise Tides
		void InitialiseTides(double TidalAmplitude, double TidalPeriod);
		
		//Initialise Waves
		void InitialiseWaves(double WaveHeight_Mean, double WaveHeight_StD, double WavePeriod_Mean, double WavePeriod_StD);
		
		//Update Sea Level
		void UpdateSeaLevel(double SLRRate);
		
		//Sample a wave
		void GetWave();
		
		/// @brief Launch the main program loop to evolve Hiro coast
		/// @details This function evolves a rocky coastal platform through time.
		///	@author Martin D. Hurst 
		/// @date 27/02/2017
		void EvolveCoast();
		
		void CalculateBackWearingForce();
		void CalculateDownWearingForce();
		void ErodeBackWearing();
		void ErodeDownWearing();
		void MassFailure();
		void Weathering();
		
		/// @brief Writes the platform morphology to file
		/// @details This function writes the elevations of the platform surface at the current time to
		///   a file. If the file exists, this is appended.
		///	@author Martin D. Hurst 
		/// @date 27/02/2017
		void WriteProfile(string OutputFileName, double Time);
		void WriteErosion(string OutputFileName, double Time);
		
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
