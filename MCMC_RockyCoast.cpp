/*------------------------------------------------------------------------

	MCMC_RockyCoast.cpp
	
	MCMC Object to evolve a coastal profile many times in in order to find 
	the most likely history of coastal retreat in keeping with 10Be CRN 
	measurements across a coastal platform.
	
	Built around code developed by Simon M. Mudd, George Hilley and Hurst 
	to invert hillslope topography for most likely boundary conditions at the 
	Dragon's Back Pressure Ridge (see Hurst et al. 2013).

	Martin D. Hurst, British Geological Survey, March 2015

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
#include "MCMC_RockyCoast.hpp"

using namespace std;

#ifndef MCMC_RockyCoast_CPP
#define MCMC_RockyCoast_CPP


void MCMC_RockyCoast::Initialise()
{
  /* initialise an empty Markov Chain Monte Carlo object
      Martin Hurst, March 2015 */
  cout << "unable to initialise, TestObject is an empty object" << endl;
  exit(EXIT_SUCCESS);  
}

void MCMC_RockyCoast::Initialise(char* CRNDatafile)
{
  /* initialise an empty Markov Chain Monte Carlo object
 	   Martin Hurst, March 2015 */
 	   
  cout << "unable to initialise with char*, TestObject is an empty object" << endl;
  exit(EXIT_SUCCESS);  
}

void MCMC_RockyCoast::Initialise(char* CRNDatafile, char* PlatformXSectionDatafile)
{
  /* initialise a Markov Chain Monte Carlo object with CRN data and the platform XSection
      Martin Hurst, March 2015 */
  
  //Declare temp variables
	char Dummy[32];
	float TempXData, TempCRNConcData, TempCRNConcErrorData, TempTopoXData, TempTopoZData;
  
  //Generate input filestream and read data into vectors
	ifstream ReadCRNDataFile(CRNDatafile);
	if (!ReadCRNDataFile)
	{
	  printf("MCMC_Coast::%s line %d: Input CRN data file \"%s\" doesn't exist\n\n", __func__, __LINE__, CRNDatafile);
	  exit(EXIT_SUCCESS);
	}
	
	ReadCRNDataFile >> Dummy >> Dummy >> Dummy >> Dummy;
	while(ReadCRNDataFile >> Dummy >> TempXData >> TempCRNConcData >> TempCRNConcErrorData)
	{
    XData.push_back(TempXData);
    CRNConcData.push_back(TempCRNConcData);
    CRNConcErrorData.push_back(TempCRNConcErrorData);
  }

  
  //Generate input filestream and read topography into vectors
  ifstream ReadTopoDataFile(PlatformXSectionDatafile);
	if (!ReadTopoDataFile)
	{
	  printf("MCMC_Coast::%s: line %d Input topographic data file \"%s\" doesn't exist\n\n", __func__, __LINE__, PlatformXSectionDatafile);
	  exit(EXIT_SUCCESS);
	}
	ReadTopoDataFile >> Dummy >> Dummy;
  while(ReadTopoDataFile >> TempTopoXData >> TempTopoZData)
  {
    TopoXData.push_back(TempTopoXData);
    TopoZData.push_back(TempTopoZData);
  }
    
  NTopoData = TopoXData.size();
  NData = XData.size();
  
  RockyCoastCRN MCMCRockyCoastCRN = RockyCoastCRN();  
}

long double MCMC_RockyCoast::CalculateLikelihood()
{
	/* Function to calculate the likelihood by comparing measured and modelled data
	   Martin Hurst, March 2015 */
	   
	//declarations
	double DiffX, Scale;
	long double Likelihood = 1.L;

  //Work out the resulting model predictions
	XDataModel = MCMCPlatformCRN.get_X();
	CRNConcModel = MCMCPlatformCRN.get_SurfaceN();
	vector<double> NModel(NData);
	vector<double> Residuals(NData);
	
	//Interpolate to sample locations
	for (int i=0; i<NData; ++i)
	{
	  //Take X value of sample and interpolate to get model results at this point
	  int j=0;
	  while ((XDataModel[j]-XData[i]) < 0) ++j;
	  DiffX = XDataModel[j]-XData[i];
    Scale = DiffX/(XDataModel[j]-XDataModel[j-1]);
  
    //Get Interpolated N value
    NModel[i] = CRNConcModel[j]-Scale*(CRNConcModel[j]-CRNConcModel[j-1]);
	}
	
	//Calculate likelihood
	for (int i=0; i<NData; i++)
	{
		Residuals[i] = (CRNConcData[i]-NModel[i])*(CRNConcData[i]-NModel[i]);
		Likelihood *= exp(-(fabs(Residuals[i]))/(CRNConcErrorData[i]*CRNConcErrorData[i]));
		
	}
	return Likelihood;
}



long double MCMC_RockyCoast::RunCoastIteration(double RetreatRate1_Test, double RetreatRate2_Test, double ChangeTime_Test, double BeachWidth_Test, double ElevInit_Test, int RetreatType)
{
	/* runs a single instance of the RockyCostCRN model, then report 
	    the likelihood of the parameters 
	    Martin Hurst, March 2015 */
	    
	//Run a coastal iteration
	int WriteResultsFlag = 0;
	MCMCPlatformCRN.UpdateParameters(RetreatRate1_Test, RetreatRate2_Test, ChangeTime_Test, BeachWidth_Test, ElevInit_Test);
	MCMCPlatformCRN.RunModel(RetreatType,WriteResultsFlag);
	
	//Calculate likelihood
	return CalculateLikelihood();
}

void MCMC_RockyCoast::RunMetropolisChain(int NIterations, char* ParameterFilename, char* OutFilename)
{
  /* Run the metropolis algorithm along a chain with NIterations
     
     ParameterFilename is the name of the parameter file containing 
      Retreat Type, Minimum, Maximum, Standard Deviation (for Metropolis search) and 
      Initial Values for each parameter, the Platform Gradient, 
      Cliff Height and Tidal Amplitude (see example)
     
     Prints to the chain results to file 'OutFilename'
     
     Martin Hurst, January 2015 */
  
	//Declarations
	long double LastLikelihood = 0.L;			//Last accepted likelihood
	long double NewLikelihood = 0.L;				//New likelihood
	long double LikelihoodRatio = 0.L;			//Ratio between last and new likelihoods
	double AcceptanceProbability; //New iteration is accepted if likelihood ratio exceeds
	
	int RetreatType;        //Style of cliff retreat 0 = single rate, 1 = step change in rates, 2 = linear change in rates
	
	double ElevInit = 0;    //initial elev init value set to zero
	
	int NAccepted = 0;      //count accepted parameters
	int NRejected = 0;      //count rejected parameters
	
	double Rand1, Rand2;    //For generating random numbers
	
	//Holders to define parameter space	
	double  RetreatRate1_New, RetreatRate1_Old, RetreatRate1_Min, RetreatRate1_Max, RetreatRate1_Std, RetreatRate1_Init,
	        RetreatRate2_New, RetreatRate2_Old, RetreatRate2_Min, RetreatRate2_Max, RetreatRate2_Std, RetreatRate2_Init,
	        ChangeTime_New, ChangeTime_Old, ChangeTime_Min, ChangeTime_Max, ChangeTime_Std, ChangeTime_Init,
	        BeachWidth_New, BeachWidth_Old, BeachWidth_Min, BeachWidth_Max, BeachWidth_Std, BeachWidth_Init,
	        PlatformGradient, CliffHeight, TidalAmplitude;
	double ElevInit_New;
	double dRR1, dRR2, dCT, dBW;  // change in parameter values for RetreatRate1, RetreatRate2, ChangeTime and BeachWidth
	double MeanChange = 0.;       //Change in parameter values centred on zero to allow changes in both directions (pos and neg)
	
	char Dummy[32];
	
	
	//Initialise seed for random number generation
	int RandomSeed = 1;
	srand(RandomSeed);
	
	//Create datafile out and write ParameterFilename
	ofstream ChainFileOut(OutFilename);
	ChainFileOut << "ParameterFile: " << ParameterFilename << endl;
	ChainFileOut  << "i RetreatRate1_New RetreatRate2_New ChangeTime_New BeachWidth_New ElevInit_New NewLikelihood LastLikelihood NAccepted NRejected" << endl;
	
	//Read in parameters for monte carlo run from parameter file
	//Min and max values for paramters from param file
	ifstream ParamFileIn(ParameterFilename);
	if (!ParameterFilename)
	{
	  printf("MCMC_Coast::%s: line %d Input parameter data file \"%s\" doesn't exist\n\n", __func__, __LINE__, ParameterFilename);
	  exit(EXIT_SUCCESS);
	}
	
	ParamFileIn >> Dummy >> RetreatType
	            >> Dummy >> RetreatRate1_Min >> Dummy >> RetreatRate1_Max >> Dummy >> RetreatRate1_Std >> Dummy >> RetreatRate1_Init
	            >> Dummy >> RetreatRate2_Min >> Dummy >> RetreatRate2_Max >> Dummy >> RetreatRate2_Std >> Dummy >> RetreatRate2_Init
	            >> Dummy >> ChangeTime_Min   >> Dummy >> ChangeTime_Max >> Dummy >> ChangeTime_Std >> Dummy >> ChangeTime_Init
	            >> Dummy >> BeachWidth_Min   >> Dummy >> BeachWidth_Max >> Dummy >> BeachWidth_Std >> Dummy >> BeachWidth_Init
	            >> Dummy >> PlatformGradient >> Dummy >> CliffHeight >> Dummy >> TidalAmplitude;
	
	ParamFileIn.close();
	
	//Initialise RockyCoastCRN object
	MCMCPlatformCRN = RockyCoastCRN(RetreatRate1_Init, RetreatRate2_Init, ChangeTime_Init, BeachWidth_Init, PlatformGradient, CliffHeight, ElevInit, TidalAmplitude);

	/*  start the chain with a guess this guess is a very coarse approximation of what the 'real' values 
	    might be. The Metropolis algorithm will sample around this */
	RetreatRate1_New = RetreatRate1_Init;
	RetreatRate2_New = RetreatRate2_Init;
	ChangeTime_New = ChangeTime_Init;
	BeachWidth_New = BeachWidth_Init;

	//Get ElevInit from Beach Width
  	//Find index that minimises TopoXData and BeachWidth
  	int i=0;
  	double DiffX, Scale;
  
  	while ((TopoXData[i]-BeachWidth_New) < 0) ++i;
  	DiffX = TopoXData[i]-BeachWidth_New;
  	Scale = DiffX/(TopoXData[i]-TopoXData[i-1]);
  
  	//Get Interpolated Z value
  	ElevInit_New = TopoZData[i]+Scale*(TopoZData[i]-TopoZData[i-1]);
   	
	//Run a single coastal iteration to get the initial Likelihood for the initial parameters
	LastLikelihood = RunCoastIteration(RetreatRate1_New, RetreatRate2_New, ChangeTime_New, BeachWidth_New, ElevInit_New, RetreatType);
	
	//set old parameters for comparison and updating
	RetreatRate1_Old = RetreatRate1_New;
	RetreatRate2_Old = RetreatRate2_New;
	ChangeTime_Old = ChangeTime_New;
	BeachWidth_Old = BeachWidth_New;
	
	int Accept = 0;

	//Do the metropolis algorithm
	for (int j=0; j<NIterations; ++j)
	{
		fflush(stdout);
 	  	printf("Iteration %d\r",j+1);  
		
		//Update the variables following a normal distribution
		
		//First Update RetreatRate1
		Accept = 0;
		while (Accept == 0)
		{	
			Rand1 = (double)rand()/RAND_MAX; Rand2 = (double)rand()/RAND_MAX;
			dRR1 = MeanChange + RetreatRate1_Std*sqrt(-2.*log(Rand1))*cos(2.*M_PI*(Rand2));
			RetreatRate1_New = RetreatRate1_Old + dRR1;
			if ((RetreatRate1_New < RetreatRate1_Min) || (RetreatRate1_New > RetreatRate1_Max)) continue;
	  		else Accept = 1;
		}
	  	// Update RetreatRate2
	 	Accept = 0;
		while (Accept == 0)
		{	
			Rand1 = (double)rand()/RAND_MAX; Rand2 = (double)rand()/RAND_MAX;
	  		dRR2 = MeanChange + RetreatRate2_Std*sqrt(-2.*log(Rand1))*cos(2.*M_PI*(Rand2));
	  		RetreatRate2_New = RetreatRate2_Old + dRR2;
			if ((RetreatRate2_New < RetreatRate2_Min) || (RetreatRate2_New > RetreatRate2_Max)) continue;
			else Accept = 1;
	  	}
	  	// Update ChangeTime
		Accept = 0;
		while (Accept == 0)
		{	
			Rand1 = (double)rand()/RAND_MAX; Rand2 = (double)rand()/RAND_MAX;
		  	dCT = MeanChange + ChangeTime_Std*sqrt(-2.*log(Rand1))*cos(2.*M_PI*(Rand2));
		  	ChangeTime_New = ChangeTime_Old + dCT;
		  	if ((ChangeTime_New < ChangeTime_Min) || (ChangeTime_New > ChangeTime_Max)) continue;
			else Accept = 1;
		}
	  
	  	// Update BeachWidth
		Accept = 0;
		while (Accept == 0)
		{	
			Rand1 = (double)rand()/RAND_MAX; Rand2 = (double)rand()/RAND_MAX;
			dBW = MeanChange + BeachWidth_Std*sqrt(-2.*log(Rand1))*cos(2.*M_PI*(Rand2));
	  		BeachWidth_New = BeachWidth_Old+dBW;
	  		if ((BeachWidth_New < BeachWidth_Min) || (BeachWidth_New > BeachWidth_Max)) continue;
			else Accept = 1;
		}
		
		//Get ElevInit from Beach Width
	  	//Find index that minimises TopoXData and BeachWidth
	  	i=0;
	  	
		while ((TopoXData[i]-BeachWidth_New) < 0) ++i;
	  	DiffX = TopoXData[i]-BeachWidth_New;
	  	Scale = DiffX/(TopoXData[i]-TopoXData[i-1]);
	  
	  	//Get Interpolated Z value
	  	ElevInit_New = TopoZData[i]+Scale*(TopoZData[i]-TopoZData[i-1]);

		//Run a model iteration with new parameters
		NewLikelihood = RunCoastIteration(RetreatRate1_New, RetreatRate2_New, ChangeTime_New, BeachWidth_New, ElevInit_New, RetreatType);
		
		//Get the likelihood ratio
		LikelihoodRatio = NewLikelihood/LastLikelihood;
		
		//Get acceptance probability (from uniform distribution between 0 and 1)
		//This allows some false results to be accepted in order to explore the parameter space.
		AcceptanceProbability = (double)rand()/RAND_MAX;
		
		//Test for acceptance
		if (LikelihoodRatio > AcceptanceProbability)
		{
			LastLikelihood = NewLikelihood;
			++NAccepted;
			
			//Last variables become equal to new variables
			RetreatRate1_Old = RetreatRate1_New;
			RetreatRate2_Old = RetreatRate2_New;
			ChangeTime_Old = ChangeTime_New;
			BeachWidth_Old = BeachWidth_New;
		}
		else ++NRejected;

		//write the result to the output file
		ChainFileOut  << j << " " 
		              << RetreatRate1_New << " " << RetreatRate2_New << " " 
		              << ChangeTime_New << " " << BeachWidth_New << " " << ElevInit_New << " "
		              << NewLikelihood << " " << LastLikelihood << " " 
		              << NAccepted << " " << NRejected << endl;
    }
    
	ChainFileOut.close();
	}
	
#endif
