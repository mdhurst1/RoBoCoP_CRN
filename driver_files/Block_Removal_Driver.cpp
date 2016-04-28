/*==============================================================

Block_Removal_Driver.hpp

A driver function to simulate the evolution of a stepped shore platform and 
predict the resultant concentrations of 10Be

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
#include <omp.h>
#include <unistd.h>
#include "../RockyCoastCRN.hpp"

using namespace std;

int main()
{
	cout << "RockyCoastCRN: Running Block Removal Experiments..." << endl;
	
	//Input parameters
	double RetreatRate1 = 0.1;            //Retreat Rate (m/yr) at the start of the model run
	double RetreatRate2 = 0.1;            //Retreat Rate (m/yr) at the end of the model run
	double ChangeTime = 0;                //Time to change retreat rates if a step change
	double PlatformGradient = 1./100.;    //Platform gradient (average if stepped platform)
	double Amp = 1.;                     //Tidal amplitude (1/2 tidal range)
	double CliffHeight = 20.;             //Cliff height for shielding
	double BeachWidth = 0.;               //Beach width (currently assumes total shielding)
	int BeachType = 0;                    //Constant beach width
	double BermHeight = 0.;               //No Beach in these runs
	double ElevInit = 0.;                 //Elevation of the platform/cliff junction
	int RetreatType = 0;	                //Scenario of retreat 0 = constant, 1 = step change, 2 = gradual change
	int SteppedPlatformFlag = 0;          //Flag for a stepped platform (1 = steps, 0=no steps)
	double StepSize = 0.0;                //size of step (0=no steps)
	
	//Create Platform CRN object
	RockyCoastCRN RockyCoastCRNModel(RetreatRate1, RetreatRate2, RetreatType, ChangeTime, BeachWidth, BeachType, BermHeight, PlatformGradient, CliffHeight, ElevInit, Amp, SteppedPlatformFlag, StepSize);

   //Run the model
  //First for no steps
  cout << "\tNo steps" << endl;
  string OutFileName = "BlockRemoval_NoSteps.pdat";
	RockyCoastCRNModel.RunModel(OutFileName);
	
	//10cm Steps
	cout << "\t10 cm steps" << endl;
	SteppedPlatformFlag = 1;
	StepSize = 0.1;
	OutFileName = "BlockRemoval_Steps10cm.pdat";
	RockyCoastCRN RockyCoastCRNModel2(RetreatRate1, RetreatRate2, RetreatType, ChangeTime, BeachWidth, BermHeight, PlatformGradient, CliffHeight, ElevInit, Amp, SteppedPlatformFlag, StepSize);
	RockyCoastCRNModel2.RunModel(OutFileName);
	
	//20cm Steps
	cout << "\t20 cm steps" << endl;
	SteppedPlatformFlag = 1;
	StepSize = 0.2;
	OutFileName = "BlockRemoval_Steps20cm.pdat";
	RockyCoastCRN RockyCoastCRNModel3(RetreatRate1, RetreatRate2, RetreatType, ChangeTime, BeachWidth, BermHeight, PlatformGradient, CliffHeight, ElevInit, Amp, SteppedPlatformFlag, StepSize);
	RockyCoastCRNModel3.RunModel(OutFileName);
	
	//40cm Steps
	cout << "\t50 cm steps" << endl;
	SteppedPlatformFlag = 1;
	StepSize = 0.5;
	OutFileName = "BlockRemoval_Steps50cm.pdat";
	RockyCoastCRN RockyCoastCRNModel4(RetreatRate1, RetreatRate2, RetreatType, ChangeTime, BeachWidth, BermHeight, PlatformGradient, CliffHeight, ElevInit, Amp, SteppedPlatformFlag, StepSize);
	RockyCoastCRNModel4.RunModel(OutFileName);
		
//100cm Steps
	cout << "\t100 cm steps" << endl;
	SteppedPlatformFlag = 1;
	StepSize = 1.0;
	OutFileName = "BlockRemoval_Steps100cm.pdat";
	RockyCoastCRN RockyCoastCRNModel5(RetreatRate1, RetreatRate2, RetreatType, ChangeTime, BeachWidth, BermHeight, PlatformGradient, CliffHeight, ElevInit, Amp, SteppedPlatformFlag, StepSize);
	RockyCoastCRNModel5.RunModel(OutFileName);
	
	cout << "\tDone" << endl << endl;
	
	return 0;
}


