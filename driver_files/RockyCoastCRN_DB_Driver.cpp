#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include "../RockyCoastCRN.hpp"

using namespace std;

int main()
{
	//Input parameters
	double RetreatRate1 = 0.01;            //Retreat Rate (m/yr) at the start of the model run
	double RetreatRate2 = 0.01;            //Retreat Rate (m/yr) at the end of the model run
	int RetreatType = 0;	              //Scenario of retreat 0 = constant, 1 = step change, 2 = gradual change
	double ChangeTime = 0;                //Time to change retreat rates if a step change (years))
	
	double PlatformGradient = 1./50.;     // Platform gradient (average if stepped platform)
	double Amp = 2;                       // Tidal amplitude (1/2 tidal range)
	double CliffHeight = 3.;             // Cliff height for shielding
	double CliffGradient = 25./35.;       // slope of the coastal bluff
	double BeachWidth = 4.;               // Beach width 
	double BermHeight = 1.;               // Height of the beach berm
    double BeachSteepnessFactor = 0.6;    // Scaling factor related to grain size controlling beach steepness
	int BeachType = 0;                    // Constant beach width = 0
	double ElevInit = 2.;                 // Elevation of the platform/cliff junction
	double SeaLevelRise = 0.001;          // Rate of sea level rise
	
	int SteppedPlatformFlag = 0;          // Flag for a stepped platform (1 = steps, 0=no steps)
	double StepSize = 0.0;                // size of step (0=no steps)
	
	//Set up which nuclides to track, 10 is 10Be, 14 is 14C, 26 is 26Al, 36 is 36Cl
	vector<int> WhichNuclides;
	WhichNuclides.push_back(10);
	
	//Create Platform CRN object
	RockyCoastCRN RockyCoastCRNModel(RetreatRate1, RetreatRate2, RetreatType, ChangeTime, BeachWidth, BeachType, BermHeight, BeachSteepnessFactor, PlatformGradient, CliffHeight, CliffGradient, ElevInit, Amp, SeaLevelRise, SteppedPlatformFlag, StepSize, WhichNuclides);

    // Read in sea level data ? //comment out these two lines for constant RSL change
    string RSLFilename = "dunbar_RSL_bradley_model.data";
    RockyCoastCRNModel.InitialiseRSLData(RSLFilename);
    
    //Run the model
    //First for no steps
    string OutFileName = "Dunbar_1cm.dat";
	bool Flag = RockyCoastCRNModel.RunModel(OutFileName);
	if (Flag == true) cout << "RetreatRate modified during model run due to sea level rise being too rapid" << endl;
	cout << endl;
	
	return 0;
}
Wayne Stephenson

