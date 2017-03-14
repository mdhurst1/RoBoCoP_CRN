/*------------------------------------------------------------------------

	RoBoCoP.hpp
	
	ROck and BOttom COastal Profile (RoBoCop) Model
	
	Martin D. Hurst, British Geological Survey, February 2016

	Copyright (C) 2015, Martin Hurst

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

/** @file RoBoCoP.hpp
@author Martin D. Hurst, British Geological Survey

@version Version 0.0.1
@brief RoBoCoP object for simulating the evolution of rocky coastal profile
@details This object contains a simple geomorphic model for the evolution of a rocky
  coastal profile 
*/

/**
@mainpage
This is the documentation for the RoBoCoP model
These pages describe the software.
@author Martin D. Hurst, British Geological Survey

------------------------------------------------------------------------*/

#ifndef RoBoCoP_HPP
#define RoBoCoP_HPP

#include <vector>

using namespace std;

/*/////////////////////////////////////////////////////////////////////////////////////////
//TEMPLATES
/////////////////////////////////////////////////////////////////////////////////////////*/

//Dummy class for RockyCoastCRN
class RockyCoastCRN;

///@brief Main coastal platform object.
class RoBoCoP
{
  friend class RockyCoastCRN;
  
	private:
	
	  /* MEMBER DECLARATIONS */
		
		// SPATIAL DOMAIN DECLARATIONS
		int NoNodes;	                      // Number of nodes across the coastline
		double NDV;                         // No data value
		double dZ;                          // Vertical spacing of nodes
		vector<double> X;		          			// cross shore distance (m)
		vector<double> Z;										// elevation (m)
		vector<double> Erosion;             //Total erosion at each elevation;
		
		// SEA LEVEL DECLARATIONS
		vector<double> RSLTime;             //Times for relative sea level elevations
		vector<double> RSLRate;             //Relative sea level elevations, will be length[1] if constant?
		double SeaLevel;
		
		// TIDES DECLARATIONS
		double TidalPeriod;                 //Tidal Period
		double TidalAmplitude;              //Tidal Amplitude
		vector<double> TideLevels;          //Vector of tide levels
		vector<double> WaterLevels;         //vector containing water levels
		vector<double> WaterDepths;         //vector containing water depths
		int NTideValues;                //number of tidal values (.size() of tidelevels vector)
		
		// WAVE DECLARATIONS
		double MeanWavePeriod;
		double StdWavePeriod;
		double MeanWaveHeight;
		double StdWaveHeight;
		double WavePeriod;
		double BreakingWaveHeight;
		double BreakingWaveWaterDepth;
		
		//CLIFF POSITION DECLARATIONS
		double CliffPositionX;        //tracks the cliff position in X
		int CliffPositionInd;         //tracks the index of the cliff position in X
				
		//TIME CONTROL PARAMETERS
		double Time, MaxTime, dt;
					
    /* FUNCTION DECLARATIONS */

		//Initialise Functions
		void Initialise();
		void Initialise(double dZ);
		void Initialise(double dZ, double PlatformGradient);
		void Initialise(double dZ, double PlatformGradient, double CliffPositionX, double TimeInterval=1.);
		void Initialise(string FileNameIn);
		
	protected:

	public:
	
		/* PlatformCRN Initialisation functions */

		/// @brief Empty initialisation function for RoBoCoP
	  ///	@author Martin D. Hurst 
    /// @date 25/02/2016
		RoBoCoP()	{
			Initialise();
		}
		
		RoBoCoP(double dZ)
		{
			Initialise(dZ);
		}
		
    RoBoCoP(double dZ, double PlatformGradient)
		{
			Initialise(dZ, PlatformGradient);
		}
		
		RoBoCoP(double dZ, double PlatformGradient, double CliffPositionX, double TimeInterval=1.)
		{
			Initialise(dZ, PlatformGradient, CliffPositionX, TimeInterval);
		}
		
		RoBoCoP(string FileNameIn)
		{
			Initialise(FileNameIn);
		}
		
		//Initialise Tides
		void InitialiseTides(double TidalAmplitude, double TidalPeriod);
		
		//Initialise Waves
		void InitialiseWaves(double WaveHeight_Mean, double WaveHeight_StD, double WavePeriod_Mean, double WavePeriod_StD);
		
		//Update Sea Level
		void UpdateSeaLevel(double SLRRate);
		
		//Sample a wave
		void GetWave();
		
		/// @brief Launch the main program loop to evolve RoBoCoP coast
		/// @details This function evolves a rocky coastal platform through time.
		///	@author Martin D. Hurst 
		/// @date 25/02/2016
		void EvolveCoast();
		
		/// @brief Writes the platform morphology to file
		/// @details This function writes the elevations of the platform surface at the current time to
		///   a file. If the file exists, this is appended.
		///	@author Martin D. Hurst 
		/// @date 25/02/2016
		void WriteProfile(string OutputFileName, double Time);
		void WriteErosion(string OutputFileName, double Time);
		
		/// @brief Get X coordinates
		/// @return X coordinates
		///	@author Martin D. Hurst 
		/// @date 25/02/2016
		vector<double> get_X() const { return X; }
		
		/// @brief Get surface CRN concentration
		/// @return Surface CRN concentration
		///	@author Martin D. Hurst 
		/// @date 25/02/2016
		vector<double> get_Z() const { return Z; }
		
};

#endif
