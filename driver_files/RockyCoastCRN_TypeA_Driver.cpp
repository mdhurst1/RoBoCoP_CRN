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
#include "../RockyCoastCRN.hpp"

using namespace std;

int main()
{
	//Input parameters
	double RetreatRate1 = 0.1;            //Retreat Rate (m/yr) at the start of the model run
	double PlatformGradient = 1./100.;    		//Platform gradient (average if stepped platform)
	double Amp = 1;                     //Tidal amplitude (1/2 tidal range)
	double CliffHeight = 0.;             //Cliff height for shielding
	double BeachWidth = 0.;               //Beach width (currently assumes total shielding)
	int BeachType = 0;
	double BermHeight = 0.;
	double ElevInit = 0.;                 //Elevation of the platform/cliff junction
	int SteppedPlatformFlag = 0;          //Flag for a stepped platform (1 = steps, 0=no steps)
	double StepSize = 0.0;                //size of step (0=no steps)

	//Create Platform CRN object
	RockyCoastCRN RockyCoastCRNModel(RetreatRate1, BeachWidth, BeachType, BermHeight, PlatformGradient, CliffHeight, ElevInit, Amp, SteppedPlatformFlag, StepSize);

	//Run the model
	//First for no steps
	string OutFileName = "TypeA_Platform.xn";
	RockyCoastCRNModel.RunModel(OutFileName);
	
	
	cout << endl;
	
	return 0;
}


