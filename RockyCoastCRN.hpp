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

#ifndef RockyCoastCRN_HPP
#define RockyCoastCRN_HPP

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

/*/////////////////////////////////////////////////////////////////////////////////////////
//TEMPLATES
/////////////////////////////////////////////////////////////////////////////////////////*/

///@brief Main coastal platform object.
class RockyCoastCRN
{
	friend class RoBoCoP;
	
	private:
		vector<double> X;	  //cross shore distance (m)
		vector<double> Z;		//elevation (m)
		
		int NXNodes;	//Number of nodes across the coastline
		int NZNodes;	//Number of nodes depth
		
		double dX;    //Nodes spacing in cross shore (m)
		double dY;    //Node spacing in vertical (m)
				
		vector<double> SurfaceN;						//CRN surface concentrations	(a/g)
		
		vector<double> PlatformElevation;		  //Platform Surface Elevations (m)
		vector<double> PlatformElevationOld;	//Platform Surface Elevations Old (m)
		vector<double> SurfaceElevation;		  //Surface Elevations (including beach cover) (m)
		
		vector< vector<double> > N;					//concentration of nuclides (a/g) as a function of position and depth
		
		vector<double> GeoMagTime;
		vector<double> GeoMagScalingFactors;  //Holds data from Lifton et al. (2014) Geomag model
		vector<double> RSLTime;
		vector<double> RSLRate;               //Holds data from Bradley et al. (2011) GIA model
	
		double PlatformWidth;								//Width of model domain
		double PlatformDepth;								//Depth of model domain
		double NDV;													//No Data Value placeholder
		
		//CRN
		//double P_Spal, P_Muon;		//local production rates for spallation and muons
		vector<double> P_Spal, P_Muon;
		
		//TIDES
		double TidalPeriod;
		double TidalAmplitude;
		vector<double> TideLevels;
		vector<double> WaterLevels;
		vector<double> WaterDepths;
		double NTidalValues;
		
		//RetreatRates
		double RetreatRate1,RetreatRate2;
		double RetreatRate;
		double RetreatType;
		double ChangeTime;
		double PlatformGradient;
		double CliffHeight;
		double JunctionElevation;
		double CliffPositionX;        //tracks the cliff position in X (m)
		int CliffPositionInd;         //tracks the index of the cliff position in X
		double XMin, XMax;            //Extent of the model domain in X (m)
		int ZTrackInd;
				
		//sea level parameters
		double SeaLevel, SLR;
		
		//Scaling parameters
		double GeoMagScalingFactor, TopoShieldingFactor;
		
		//Run time control parameters
		double Time, MaxTime, dt;
		
		//RetreatSyle
		int SteppedPlatform;
		double StepSize;
	
	  //Beach profile stuff
	  vector<double> BeachThickness;
	  double InitialBeachWidth; //Initial value of beach width when using a thinning beach width
	  double BeachWidth;      //Width of beach at BermHeight
	  double MeanBeachWidth;  //Mean Width of Beach when using a variable beach width
		double BermHeight;      //Height of Berm
		int BeachType;          //Style of beach evolution, 0 = fixed beach, 1 = sinusoidal beach width, 2 = thinning beach width
		double A;               //Sediment Scale Parameter in Bruun Profile (m^1/3)

	  string OutFileName;
	  
		//Initialise Function
		void Initialise();
		void Initialise(double retreatrate, double beachwidth, int beachtype, double bermheight, double platformgradient, double cliffheight, double junctionelevation, double tidalamplitude,int steppedplatform=0, double stepsize=0);
		void Initialise(double retreatrate1, double retreatrate2, int retreattype, double changetime, double beachwidth, int beachtype, double bermheight, double platformgradient, double cliffheight, double junctionelevation, double tidalamplitude,int steppedplatform=0, double stepsize=0);
		void Initialise(RoBoCoP RoBoCoPCoast);
		
		//functions to initialise platform morphology
		void InitialisePlanarPlatformMorphology();
				
		//function to retrieve topographic shielding factor
		double GetTopographicShieldingFactor(double X, double CliffHeight);
		
		//functions to read and retrieve geomagnetic scaling factor
		void InitialiseGeomagData();		
		double GetGeoMagScalingFactor(double Time);
		
		//factor to retrieve a rate of sea level rise
		void InitialiseRSLData();
		double GetSeaLevelRise(double Time);
		
		//function to get sinusoid beachwidth through time
		void GetSinWaveBeachWidth(double Time);
		void GetThinningBeachWidth(double Time);
		
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
		/// @param junctionelevation Elevation of platform at the cliff
		/// @param tidalamplitude Amplitude of diurnal tides
	  ///	@author Martin D. Hurst 
    /// @date 14/09/2015
		RockyCoastCRN(double retreatrate1, double retreatrate2, int retreattype, double changetime, double beachwidth, int beachtype, double bermheight, double platformgradient, double cliffheight, double junctionelevation, double tidalamplitude, int steppedplatform=0, double stepsize=0)
		{
			Initialise(retreatrate1, retreatrate2, retreattype, changetime, beachwidth, beachtype, bermheight, platformgradient, cliffheight, junctionelevation, tidalamplitude, steppedplatform, stepsize);
		}
		
		/// @brief Initialisation function for a single retreat rate scenario.
		/// @param retreatrate Rate of cliff retreat (m/yr)
		/// @param beachwidth Width of the beach (constant) blocking CRN production
		/// @param platformgradient Gradient (dz/dx) of coastal platform
		/// @param cliffheight Height of retreating cliff (constant)
		/// @param junctionelevation Elevation of platform at the cliff
		/// @param tidalamplitude Amplitude of diurnal tides
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		RockyCoastCRN(double retreatrate, double beachwidth, int beachtype, double platformgradient, double cliffheight, double junctionelevation, double tidalamplitude, int steppedplatform=0, double stepsize=0)
		{
			Initialise(retreatrate, beachwidth, beachtype, platformgradient, cliffheight, junctionelevation, tidalamplitude, steppedplatform, stepsize);
		}
		
		/// @brief Initialisation function with friend class RoBoCoP as the morphological model
		/// @param RoBoCoP RoBoCoPCoast a RoBoCoP coastal morphology object 
		///	@author Martin D. Hurst 
    /// @date 8/3/2015
		RockyCoastCRN(RoBoCoP RoBoCoPCoast)
		{
			Initialise(RoBoCoPCoast);
		}
		
		/// @brief Empty initialisation function, will throw an error.
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		RockyCoastCRN()
		{
			Initialise();
		}
		
		/// @brief function to initialise the tides
		/// @details This function initialises tides as a single cosine wave with fixed amplitude and
		///   period. 
		/// @param A Tidal Amplitude (metres)
		/// @param T Tidal Period (hours)
		///	@author Martin D. Hurst 
    /// @date 15/03/2016
		void InitialiseTides(double A, double T);
				
		/// @brief Update RockyCoastCRN object parameters for a new model run
		/// @param RetreatRate1_Test New first rate of cliff retreat (m/yr)
		/// @param RetreatRate2_Test second rate of cliff retreat (m/yr)
		/// @param ChangeTime_Test New time to switch from retreatrate1 to retreatrate2 (years BP)
		/// @param BeachWidth_Test New width of the beach (constant) blocking CRN production
		/// @param JunctionElevation_Test New elevation of platform at the cliff
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		void UpdateParameters(double RetreatRate1_Test, double RetreatRate2_Test, double ChangeTime_Test, double BeachWidth_Test, double JunctionElevation_Test);
		
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
		void RunModel(string OutFileName, int WriteResultsFlag=1);
		
		/// @brief Determine current retreat rate
		/// @details This function determines the cliff retreat rate for the current timestep based on
		///   the selected scenario/mode of cliff retreat. Retreat type = 0 is constant retreat rate
		///   1 = step change in retreat rate at time ChangeTime and 2 = gradual change in reteat rate
		///   through time.
		///	@author Martin D. Hurst 
    /// @date 09/02/2016
		void GetRetreatRate();
		
		/// @brief Updates the CRN concentrations at the platform surface and at depth
		/// @details This function calculates the accumulation of 10Be in the platform surface and
		///   at depth by both spallation and muogenic production. 
		/// @param TimeInterval the time step in years
		///	@author Martin D. Hurst 
    /// @date 09/02/2016
		void UpdateCRNs();

		/// @brief Updates the platform morphology
		/// @details This function calculates the amount of platform downwear and updates the elevations
		///   of the platform surface. Currently just does gradual uniform downwear or step-retreat.
		///	@author Martin D. Hurst 
    /// @date 09/02/2016
		void UpdateEquillibriumMorphology();
		
		/// @brief Updates the platform morphology
		/// @details This function calculates the amount of platform downwear and cliff retre updates the elevations
		///   of the platform surface based on an iteration of the RoBoCoP coast object.
		///	@author Martin D. Hurst
		/// @param RoBoCoPCoast A RoBoCoP Coastal model object
    /// @date 14/03/2016
		void UpdateMorphology(RoBoCoP RoBoCoPCoast);
		
		/// @brief Writes the platform morphology to file
		/// @details This function writes the elevations of the platform surface at the current time to
		///   a file. If the file exists, this is appended.
		///	@author Martin D. Hurst 
    /// @date 09/02/2016
		void WriteProfile(string OutputFileName, double Time);
		
		/// @brief Writes the platform morphology to file
		/// @details This function writes the CRN concnetrations on the platform surface at the current time to
		///   a file. If the file exists, this is appended.
		///	@author Martin D. Hurst 
    /// @date 15/03/2016
		void WriteCRNProfile(string OutputFileName, double Time);
		
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
