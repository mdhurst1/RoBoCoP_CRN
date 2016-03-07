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
  double PlatformGradient = 1./50.;
  double CliffPositionX = 0.;
  
  //initialise RoBoCoP
  RoBoCoP RoBoCoP = RoBoCoP(dZ);
	
	//Initialise Tides
	double TidalAmplitude = 2.5;
	double TidalPeriod = 12.42;
	RoBoCop.InitialiseTides(TidalAmplitude, TidalPeriod);
	
	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight = 1.2;
	double WavePeriod = 6.;
	
	//a few blank lines to finish
	cout << endl << endl;
	
	// Time Control
	int EndTime = 10000.;
	int Time = 0.;
	int TimeInterval = 0.1;
	
	//Print Control
	int PrintInterval = 100.;
	int PrintTime = Time;
	
	//Loop through time
	while (Time < EndTime)
	{
	  //Evolve the coast
	  EvolveCoast(TimeInterval);
	  
	  //print?
	  if (Time >= PrintTime)
	  {
	    WriteCoastFile();
	    PrintTime += PrintInterval;
	  }
	  
	  //update time
	  Time += TimeInterval;
	}
}


