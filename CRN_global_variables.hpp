/*-----------------------------------------------------------------------------------------

Global parameters for Rocky Coast Cosmogenic Radionuclide Model

Martin Hurst
British Geological Survey

December 2014

-----------------------------------------------------------------------------------------*/

/**
@file   CRN_global_variables.hpp
@author Martin D. Hurst (mhurst@bgs.ac.uk)
@date   14/09/2015
@brief  Global variables for modelling cosmogenic nuclide accumulation. 
*/
#include <cmath>
#include <vector>


using namespace std;

#ifndef CRN_GLOBALS_HPP
#define CRN_GLOBALS_HPP

//Sea level high latitude production rates
///@brief   10Be spallation reference production rate.
///@details Spallation (a/g/yr) calibrated 10Be production rate (Lifton et al. 2014).
float Po_Spal = 4.007166;

///@brief   10Be muogneic reference produciton rate.
///@details Total muogenic production rate (a/g/yr) following Braucher et al. (2013).
float Po_Muon = 0.028;

// Attenuation Lengths
///@brief Spallogenic attenuation length (kg/m^2).
float Lamb_Spal = 1600.0;
///@brief Muogenic attenuation length (kg/m^2) following Braucher et al. (2013).
float Lamb_Muon = 42000.0;	

//Densities
///@brief Rock (chalk) density (kg/m^3) from Mortimore et al. (2004).
float rho_r = 1800.0;
///@brief	Density of sea water (kg/m^3)	
float rho_w = 1035.0;

//Decay length scales
float z_rs = Lamb_Spal/rho_r;		//Decay length scale chalk spallation
float z_ws = Lamb_Spal/rho_w;		//Decay length scale sea water spallation
float z_rm = Lamb_Muon/rho_r;		//Decay length scale chalk muons
float z_wm = Lamb_Muon/rho_w;		//Decay length scale sea water muons

///@brief Half life of 10Be (Korschineck et al. 2010).
float Lambda = 4.99*pow(10.,-7);

//float Po_Spal = 4.34;				//Spallation (a/g/yr) to match Regard et al (2012)
//float Po_Fast = 0.039;			//Fast Muons (a/g/yr) to match Regard et al (2012)
//float Po_Slow = 0.012;			//Slow Muons (a/g/yr) to match Regard et al (2012)
//float Lamb_Fast = 15000.0;	//Slow Muons (kg/m^2) to match Regard et al (2012)
//float Lamb_Slow = 43200.0;	//Fast Muons (kg/m^2) to match Regard et al (2012)
//float z_rfast = Lamb_Fast/rho_r;	//Decay length scale chalk fast muons , to match Regard et al (2012)
//float z_wfast = Lamb_Fast/rho_w;	//Decay length scale sea water fast muons, to match Regard et al (2012)
//float z_rslow = Lamb_Slow/rho_r;	//Decay length scale chalk slow muons, to match Regard et al (2012)
//float z_wslow = Lamb_Slow/rho_w;	//Decay length scale sea water slow muons, to match Regard et al (2012)


#endif
