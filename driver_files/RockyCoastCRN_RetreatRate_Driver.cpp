/*==============================================================

RockyCoastCRN_RetreatRate_Driver.hpp

A driver function to simulate the evolution of a shore platform and 
predict the resultant concentrations of 10Be in the platform

Developed by:
Martin D. Hurst

Copyright (C) 2016, Martin Hurst

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

==============================================================*/
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
	double RetreatRate1 = 0.07455212;            //Retreat Rate (m/yr) at the start of the model run
	double RetreatRate2 = 0.01156387;            //Retreat Rate (m/yr) at the end of the model run
        int RetreatType = 2;	             //Scenario of retreat 0 = constant, 1 = step change, 2 = gradual change

	double ChangeTime = 0;                //Time to change retreat rates if a step change
	double PlatformGradient = 1./50.;    //Platform gradient (average if stepped platform)
	double Amp = 4.;                     	//Tidal amplitude (1/2 tidal range)
	double CliffHeight = 35.;             //Cliff height for shielding
        double CliffGradient = 25./35.;       //slope of the coastal bluff
	double BeachWidth = 2.;               //Beach width (currently assumes total shielding)
	double BermHeight = 5.;               //Elevation of Beach Berm above junction (set beachwidth = 0 and bermheight = 0 for no beaches)
        double BeachSteepnessFactor = 0.6;    // Scaling factor related to grain size controlling beach steepness
	int BeachType = 0;                    //0 = constant beach width, 1 = sinusoidal beach width, 2 = thinning beach
	double ElevInit = 2.;                 //Elevation of the platform/cliff junction
        double SeaLevelRise = 0.001;           //Rate of sea level rise


	int SteppedPlatformFlag = 0;          //Flag for a stepped platform (1 = steps, 0=no steps)
	double StepSize = 0.0;                //size of step (0=no steps)

        //Set up which nuclides to track, 10 is 10Be, 14 is 14C, 26 is 26Al, 36 is 36Cl
	vector<int> WhichNuclides;
	WhichNuclides.push_back(10);
	
	
	//Create Platform CRN object
	RockyCoastCRN RockyCoastCRNModel (RetreatRate1, RetreatRate2, RetreatType, ChangeTime, BeachWidth, BeachType, BermHeight, BeachSteepnessFactor, PlatformGradient, CliffHeight, CliffGradient, ElevInit, Amp, SeaLevelRise, SteppedPlatformFlag, StepSize, WhichNuclides);

    // Read in sea level data ? //comment out these two lines for constant RSL change
    string RSLFilename = "bideford_RSL_bradley_model.data";
    RockyCoastCRNModel.InitialiseRSLData(RSLFilename);

  //Run the model for RetreatRate 5_1 + 5_2
  string OutFileName = "RetreatRate_1_test.xzn";
	RockyCoastCRNModel.RunModel(OutFileName);
	
	//RetreatRate 5_1 + M_2
	RetreatRate1 = 0.07455212;
	RetreatRate2 = 0.01920925;
	OutFileName = "RetreatRate_2_test.xzn";
	RockyCoastCRNModel.UpdateParameters(RetreatRate1,RetreatRate2,ChangeTime,BeachWidth);
	RockyCoastCRNModel.RunModel(OutFileName);
	
	//RetreatRate 5_1 + 95_2
	RetreatRate1 = 0.07455212;
	RetreatRate2 = 0.02963008;
	OutFileName = "RetreatRate_3_test.xzn";
	RockyCoastCRNModel.UpdateParameters(RetreatRate1,RetreatRate2,ChangeTime,BeachWidth);
	RockyCoastCRNModel.RunModel(OutFileName);
	
	//RetreatRate M_1 + 5_2
	RetreatRate1 = 0.08327285;
	RetreatRate2 = 0.01156387;
	OutFileName = "RetreatRate_4_test.xzn";
	RockyCoastCRNModel.UpdateParameters(RetreatRate1,RetreatRate2,ChangeTime,BeachWidth);
	RockyCoastCRNModel.RunModel(OutFileName);
	
	//RetreatRate M_1 + M_2
	RetreatRate1 = 0.08327285;
	RetreatRate2 = 0.01920925;
	OutFileName = "RetreatRate_5_test.xzn";
	RockyCoastCRNModel.UpdateParameters(RetreatRate1,RetreatRate2,ChangeTime,BeachWidth);
	RockyCoastCRNModel.RunModel(OutFileName);

        //RetreatRate M_1 + 95_2
	RetreatRate1 = 0.08327285;
	RetreatRate2 = 0.02963008;
	OutFileName = "RetreatRate_6_test.xzn";
        RockyCoastCRNModel.UpdateParameters(RetreatRate1,RetreatRate2,ChangeTime,BeachWidth);
	RockyCoastCRNModel.RunModel(OutFileName);

        //RetreatRate 95_1 + 5_2
	RetreatRate1 = 0.09386583;
	RetreatRate2 = 0.01156387;
	OutFileName = "RetreatRate_7_test.xzn";
	RockyCoastCRNModel.UpdateParameters(RetreatRate1,RetreatRate2,ChangeTime,BeachWidth);
	RockyCoastCRNModel.RunModel(OutFileName);

        //RetreatRate 95_1 + M_2
	RetreatRate1 = 0.09386583;
	RetreatRate2 = 0.01920925;
	OutFileName = "RetreatRate_8_test.xzn";
	RockyCoastCRNModel.UpdateParameters(RetreatRate1,RetreatRate2,ChangeTime,BeachWidth);
	RockyCoastCRNModel.RunModel(OutFileName);

        //RetreatRate 95_1 + 95_2
	RetreatRate1 = 0.09386583;
	RetreatRate2 = 0.02963008;
	OutFileName = "RetreatRate_9_test.xzn";
	RockyCoastCRNModel.UpdateParameters(RetreatRate1,RetreatRate2,ChangeTime,BeachWidth);
	RockyCoastCRNModel.RunModel(OutFileName);
	
	cout << endl;
	
	return 0;
}


