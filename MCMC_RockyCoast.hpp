/*------------------------------------------------------------------------

	MCMC_RockyCoast.hpp
	
	MCMC Object to evolve a coastal profile many times in in order to find 
	the most likely history of coastal retreat in keeping with 10Be CRN 
	measurements across a coastal platform.
	
	Built around code developed by Hilley, Mudd and Hurst to invert 
	hillslope topography for most likely boundary conditions at the 
	Dragon's Back Pressure Ridge (see Hurst et al 2013).

	Martin D. Hurst, British Geological Survey, January 2015

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

  /** @file MCMC_RockyCoast.hpp
  @author Martin D. Hurst, British Geological Survey

  @version Version 0.0.1
  @brief MCMC_RockyCoast object for running a Markov Chain Monte Carlo inversion to get cliff
    retreat rates from CRN data on a coastal platform
  @details This object contains routines to carry out a Markov Chain Monte Carlo analysis 
    comparing measured and modelled concentrations of cosmogenic radionuclides in order to
    determine the most likely parameters and cliff retreat rates given a set of CRN concentration
    measurements. Calls the RockyCoastCRN object.
  */

  /**
  @mainpage
  This is the documentation for the MCMC_RockyCoast algorithms
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
#include "RockyCoastCRN.hpp"

using namespace std;

#ifndef MCMC_RockyCoast_HPP
#define MCMC_RockyCoast_HPP

/*/////////////////////////////////////////////////////////////////////////////////////////
//TEMPLATES
/////////////////////////////////////////////////////////////////////////////////////////*/

///@brief Main Markov Chain Monte Carlo object.
class MCMC_RockyCoast
{
	private:
		//Initialise Function
		void Initialise();
		void Initialise(char* CRNDatafile);
		void Initialise(char* CRNDatafile, char* PlatformXSectionDatafile);

		//Vectors to hold CRN data
		int NData;
		vector<double> XData;
		vector<double> CRNConcData;
		vector<double> CRNConcErrorData;
		
		//Vectors to hold Model results
		int NModelData;
		vector<double> XDataModel;
		vector<double> CRNConcModel;
		
		//Vectors to hold topographic data
		int NTopoData;
		vector<double> TopoXData;
		vector<double> TopoZData;
		
		RockyCoastCRN MCMCPlatformCRN;
		
		// calcualtes the likelihood using measured and modelled data
		long double CalculateLikelihood();
		
		// runs a single iteration of the coastal model, then reports the likelihood of the parameters
		long double RunCoastIteration(double RetreatRate1_Test, double RetreatRate2_Test, double ChangeTime_Test, double BeachWidth_Test, double ElevInit_Test, int RetreatType);

	public:
	
		/* MCMC_RockyCoast Initialisation functions */

		/// @brief Empty initialisation function, will throw an error.
		///	@author Martin D. Hurst 
    /// @date 16/09/2015
		MCMC_RockyCoast()
		{
			Initialise();
		}
		
		/// @brief Initialisation function where only 1 argument provided throws an error.
		/// @param CRNDatafile
		///	@author Martin D. Hurst 
    /// @date 16/09/2015
		MCMC_RockyCoast(char* CRNDatafile)	
		{
		  Initialise(CRNDatafile);
		}
		
		/// @brief Initialisation function for MCMC_RockyCoast objects
		/// @param CRNDatafile Data file containing observed CRN concentrations
		/// @param PlatformXSectionDatafile Data file containing platform cross section (currently not used).
		///	@author Martin D. Hurst 
    /// @date 16/09/2015
		MCMC_RockyCoast(char* CRNDatafile, char* PlatformXSectionDatafile)	
		{
		  Initialise(CRNDatafile, PlatformXSectionDatafile);
		}
		
		// this runs the metropolis algorithm along a chain with NIterations it prints to the file 'OutFilename'
		/// @brief Launch the main MCMC program loop to search for most likely parameters
		/// @details This function runs a Markov Chain Monte Carlo analysis to find the most likely
		///   combination of parameters to fit observed concentrations of 10Be from a coastal platform.
		///   The MCMC algorithm runs the RockyCoastCRN model many times (NIterations c. 200k), checks the 
		///   likelihood betweeen modelled and measured cncentrations and either accepts or rejects the 
		///   new parameters. The model sometimes accepts less likely results to allow exploration of the
		///   parameter space. 
		/// @param NIterations Number of times to run the RockyCoastCRN model in the Markov Chain
		/// @param ParamFilename Filename of file containing parameters for the MCMC analysis (see example file)
		/// @param OutFilename File to write the results of each iteration of the chain to.
		///	@author Martin D. Hurst 
    /// @date 16/09/2015
    void RunMetropolisChain(int NIterations, char* ParamFilename, char* OutFilename);

};

#endif
