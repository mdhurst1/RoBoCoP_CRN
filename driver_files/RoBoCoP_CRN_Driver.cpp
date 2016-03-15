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
#include <unistd.h>
#include "../RoBoCoP.hpp"
#include "../RockyCoastCRN.hpp"

using namespace std;

int main()
{
  //initialisation parameters
  double dZ = 0.1;
  double PlatformGradient = 1./10.;
  double CliffPositionX = 0.;
  
  //initialise RoBoCoP
  RoBoCoP PlatformModel = RoBoCoP(dZ,PlatformGradient,CliffPositionX);

	//Initialise Tides
	double TidalAmplitude = 1.;
	double TidalPeriod = 12.42;
	PlatformModel.InitialiseTides(TidalAmplitude, TidalPeriod);
	
	//Initialise RockyCoastCRN
	RockyCoastCRN PlatformCRN = RockyCoastCRN(PlatformModel);
	
	
	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight = 1.2;
	double WavePeriod = 6.;
	PlatformModel.InitialiseWaves(WaveHeight,WavePeriod);
	
	// Time Control
	double EndTime = 10000.;
	double Time = 0.;
	double TimeInterval = 1;
	
	//Print Control
	double PrintInterval = 100.;
	double PrintTime = PrintInterval;
	string OutputMorphologyFileName = "ShoreProfile.xz";
	string OutputCRNMorphologyFileName = "ShoreProfile2.xz";
	string OutputConcentrationFileName = "CRNConcentrations.xn";
	
	//print initial conditions to file first
	PlatformModel.WriteProfile(OutputMorphologyFileName, Time);
  PlatformCRN.WriteProfile(OutputCRNMorphologyFileName, Time);
  PlatformCRN.WriteCRNProfile(OutputConcentrationFileName, Time);
  
	//Loop through time
	while (Time < EndTime)
	{
	  //Evolve the coast
	  PlatformModel.EvolveCoast(TimeInterval);
	  
	  //Update the morphology inside RockyCoastCRN
	  PlatformCRN.UpdateMorphology(PlatformModel);
	  
	  //Update the CRN concentrations
	  PlatformCRN.UpdateCRNs();
    
    //update time
	  Time += TimeInterval;
	      
	  //print?
	  if (Time >= PrintTime)
	  {
	    PlatformModel.WriteProfile(OutputMorphologyFileName, Time);
	    PlatformCRN.WriteProfile(OutputCRNMorphologyFileName, Time);
	    PlatformCRN.WriteCRNProfile(OutputConcentrationFileName, Time);
	    PrintTime += PrintInterval;
	    //cout << "Time is " << Time << endl;
	  }
	  
	  
	}
	//write a gap at the end
	cout << endl << endl;
}


