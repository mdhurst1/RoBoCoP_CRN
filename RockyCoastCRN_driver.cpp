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
#include "PlatformCRN.hpp"

using namespace std;

int main()
{
	double RetreatRate1 = 0.25;
	double RetreatRate2 = 0.025;
	double ChangeTime = 0;
	double PlatformGradient = 1./60.;
	double Amp = 2.4;
	double CliffHeight = 50.;
	double BeachWidth = 10.;
	double ElevInit = 0.;
	
	PlatformCRN PlatformCRNModel(RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, PlatformGradient, CliffHeight, ElevInit, Amp);
	int RetreatType = 2;
	PlatformCRNModel.RunModel(RetreatType);
	
	return 0;
}


