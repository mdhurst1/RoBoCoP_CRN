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

using namespace std;

int main()
{
  //initialisation parameters
  double dZ = 0.1;
  //double PlatformGradient = 1./50.;
  //double CliffPositionX = 0.;
  
  //initialise RoBoCoP
  RoBoCoP PlatformModel = RoBoCoP(dZ);
	
	//Initialise Tides
	double TidalAmplitude = 1.;
	double TidalPeriod = 12.42;
	PlatformModel.InitialiseTides(TidalAmplitude, TidalPeriod);
	
	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight = 1.2;
	double WavePeriod = 6.;
	PlatformModel.InitialiseWaves(WaveHeight,WavePeriod);
	
	//a few blank lines to finish
	cout << endl << endl;
	
	// Time Control
	double EndTime = 10000.;
	double Time = 0.;
	double TimeInterval = 1;
	
	//Print Control
	double PrintInterval = 100.;
	double PrintTime = Time;
	string OutputFileName = "ShoreProfile.xz";
	string ErosionFileName = "Erosion.ez";
	
	//Loop through time
	while (Time < EndTime)
	{
	  //Evolve the coast
	  PlatformModel.EvolveCoast(TimeInterval);
	  
	  //print?
	  if (Time >= PrintTime)
	  {
	    PlatformModel.WriteProfile(OutputFileName, Time);
	    PlatformModel.WriteErosion(ErosionFileName, Time);
	    PrintTime += PrintInterval;
	    //cout << "Time is " << Time << endl;
	  }
	  
	  //update time
	  Time += TimeInterval;
	}
}


