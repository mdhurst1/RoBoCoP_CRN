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
#include "RockyCoastCRN.hpp"

using namespace std;

int main()
{
	//Input parameters
	double RetreatRate1 = 0.1;
	double RetreatRate2 = 0.1;
	double ChangeTime = 0;
	double PlatformGradient = 1./60.;
	double Amp = 2.4;
	double CliffHeight = 50.;
	double BeachWidth = 10.;
	double ElevInit = 0.;
	int RetreatType = 0;	
	int SteppedPlatformFlag = 0;
	double StepSize = 0;
	
	//Create Platform CRN object
	RockyCoastCRN RockyCoastCRNModel(RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, PlatformGradient, CliffHeight, ElevInit, Amp, SteppedPlatformFlag, StepSize);

  //Run the model
	RockyCoastCRNModel.RunModel(RetreatType);
	
	return 0;
}


