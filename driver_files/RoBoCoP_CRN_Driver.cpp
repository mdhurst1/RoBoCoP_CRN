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
	double PrintTime = Time;
	string OutputFileName = "ShoreProfile.xz";
	
	//Loop through time
	while (Time < EndTime)
	{
	  //Evolve the coast
	  PlatformModel.EvolveCoast(TimeInterval);
	  
	  //Update the morphology inside RockyCoastCRN
	  PlatformCRN.UpdateMorphology(PlatformModel);
	  
	  //Update the CRN concentrations
	  PlatformCRN.UpdateCRNs();
        
	  //print?
	  if (Time >= PrintTime)
	  {
	    PlatformModel.WriteProfile(OutputFileName, Time);
	    PlatformCRN.WriteProfile();
	    PrintTime += PrintInterval;
	    //cout << "Time is " << Time << endl;
	  }
	  
	  //update time
	  Time += TimeInterval;
	}
}


