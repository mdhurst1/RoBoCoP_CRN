/*------------------------------------------------------------------------

	PlatformCRN.hpp
	
	Object to evolve the CRN concentration across a coastal profile as a 
	function of tides, topographic shielding, relative sea level rise, 
	cliff retreat and platform downwear. For comparison with measured
	10Be CRN measurements across a coastal platform.
	
	Martin D. Hurst, British Geological Survey, December 2014

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

/** @file RockyCoastCRN.hpp
@author Martin D. Hurst, British Geological Survey

@version Version 0.0.1
@brief RockyCoastCRN object for predicting CRN concentrations on a coastal platform
@details This object contains a simple geomorphic model for the evolution of a rocky
  coastal profile and routines to predict the accumulation of cosmogenic radionuclides
  (10Be only at the moment) in the coastal platform. Model is intended to allow inversion
  to determine rates of cliff retreat from the distribution of CRNs by linkage with the
  MCMC_RockyCoast object which calls it.
*/

/**
@mainpage
This is the documentation for the RockyCoastCRN model
These pages describe the software.
@author Martin D. Hurst, British Geological Survey

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

using namespace std;

#ifndef RockyCoastCRN_HPP
#define RockyCoastCRN_HPP

/*/////////////////////////////////////////////////////////////////////////////////////////
//TEMPLATES
/////////////////////////////////////////////////////////////////////////////////////////*/

///@brief Main coastal platform object.
class RockyCoastCRN
{
	private:
		int NXNodes;	//Number of nodes across the coastline
		int NZNodes;	//Number of nodes depth
		
		vector<double> X;										//cross shore distance (m)
		vector<double> Z;										//elevation (m)
		vector<double> SurfaceN;							//Surface concentrations	(a/g)
		vector<double> SurfaceElevation;			//Surface Elevations
		vector<double> SurfaceElevationOld;	//Surface Elevations
		vector< vector<double> > N;					//concentration of nuclides (a/g)
		
		vector<double> GeoMagTime;
		vector<double> GeoMagScalingFactors;  //Holds data from Lifton et al. (2014) Geomag model
		vector<double> RSLTime;
		vector<double> RSLRate;             //Holds data from Bradley et al. (2011) GIA model
	
		double PlatformWidth;								//Wdith of model domain
		double PlatformDepth;								//depth of model domain
		double NDV;													//No Data Value placeholder
		
		//CRN
		double P_spal, P_muon;		//local production rates for spallation and muons
		
		//TIDES
		double TidalPeriod;
		double TidalAmplitude;
		vector<double> TideLevels;
		vector<double> WaterLevels;
		vector<double> WaterDepths;
		double NTidalValues;
		
		//RetreatRates
		double RetreatRate1,RetreatRate2;
		double ChangeTime;
		double PlatformGradient;
		double CliffHeight;
		double BeachWidth;
		double ElevInit;
	
		//Initialise Function
		void Initialise();
		void Initialise(double retreatrate, double platformgradient, double cliffheight, double beachwidth, double elevinit, double tidalamplitude);
		void Initialise(double retreatrate1, double retreatrate2, double changetime, double platformgradient, double cliffheight, double beachwidth, double elevinit, double tidalamplitude);
		
		//function to retrieve topographic shielding factor
		double GetTopographicShieldingFactor(double X, double CliffHeight);
		
		//function ot retrieve geomagnetic scaling factor		
		double GetGeoMagScalingFactor(double Time);
		
		//factor to retrieve a rate of sea level rise
		double GetSeaLevelRise(double Time);
			
	protected:

	public:
	
		/* PlatformCRN Initialisation functions */

		/// @brief Initialisation function for two retreat rate scenario.
		/// @param retreatrate1 First rate of cliff retreat (m/yr)
		/// @param retreatrate2 Second rate of cliff retreat (m/yr)
		/// @param changetime Time to switch from retreatrate1 to retreatrate2 (years BP)
		/// @param beachwidth Width of the beach (constant) blocking CRN production
		/// @param platformgradient Gradient (dz/dx) of coastal platform
		/// @param cliffheight Height of retreating cliff (constant)
		/// @param elevinit Elevation of platform at the cliff
		/// @param tidalamplitude Amplitude of diurnal tides
	  ///	@author Martin D. Hurst 
    /// @date 14/09/2015
		PlatformCRN(double retreatrate1, double retreatrate2, double changetime, double beachwidth, double platformgradient, double cliffheight, double elevinit, double tidalamplitude)
		{
			Initialise(retreatrate1, retreatrate2, changetime, beachwidth, platformgradient, cliffheight, elevinit, tidalamplitude);
		}
		
		/// @brief Initialisation function for a single retreat rate scenario.
		/// @param retreatrate Rate of cliff retreat (m/yr)
		/// @param beachwidth Width of the beach (constant) blocking CRN production
		/// @param platformgradient Gradient (dz/dx) of coastal platform
		/// @param cliffheight Height of retreating cliff (constant)
		/// @param elevinit Elevation of platform at the cliff
		/// @param tidalamplitude Amplitude of diurnal tides
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		PlatformCRN(double retreatrate, double beachwidth, double platformgradient, double cliffheight, double elevinit, double tidalamplitude)
		{
			Initialise(retreatrate, beachwidth, platformgradient, cliffheight, elevinit, tidalamplitude);
		}
		
		/// @brief Empty initialisation function, will throw an error.
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		PlatformCRN()
		{
			Initialise();
		}
		
		/// @brief Update RockyCoastCRN object parameters for a new model run
		/// @param RetreatRate1_Test New first rate of cliff retreat (m/yr)
		/// @param RetreatRate2_Test second rate of cliff retreat (m/yr)
		/// @param ChangeTime_Test New time to switch from retreatrate1 to retreatrate2 (years BP)
		/// @param BeachWidth_Test New width of the beach (constant) blocking CRN production
		/// @param ElevInit_Test New elevation of platform at the cliff
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		void UpdateParameters(double RetreatRate1_Test, double RetreatRate2_Test, double ChangeTime_Test, double BeachWidth_Test, double ElevInit_Test);
		
		/// @brief Launch the main program loop to evolve coast and predict CRN concentrations
		/// @details This function evolves a rocky coastal platform through time assumming gradual
		///   cliff retreat and platform downwear proportional to cliff retreat (i.e. equillibrium
		///   cliff retreat; translating a constant crossshore profile landward through time). The
		///   model takes geomagnetic scaling as a text file (generated by Lifton et al. (2014) code)
		///   and a relative sea level history text file (generated from Bradley et al. (2011) model).
		/// @param RetreatType Single (=0), step change (=1), or gradual change (=2) retreat rate scenario
		/// @param WriteResultsFlag flag to write results to file (=1) or not (=0), default is on.
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		void RunModel(int RetreatType,int WriteResultsFlag=1);
		
		/// @brief Get X coordinates
		/// @return X coordinates
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		vector<double> get_X() const { return X; }
		
		/// @brief Get surface CRN concentration
		/// @return Surface CRN concentration
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		vector<double> get_SurfaceN() const { return SurfaceN; }
		
};

#endif
