/*==============================================================

BeachWidth_Driver.hpp

A driver function to simulate the evolution of a beach fronted cliff
and shore platform and predict the resultant concentrations of 10Be
in the platform

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
	cout << "RockyCoastCRN: Running Variable Beach Width Experiments..." << endl;
	
	//Input parameters
	double RetreatRate1 = 0.1;            //Retreat Rate (m/yr) at the start of the model run
	double RetreatRate2 = 0.1;            //Retreat Rate (m/yr) at the end of the model run
	double ChangeTime = 0;                //Time to change retreat rates if a step change
	double PlatformGradient = 1./100.;    //Platform gradient (average if stepped platform)
	double Amp = 1.0;                     //Tidal amplitude (1/2 tidal range)
	double CliffHeight = 10.;             //Cliff height for shielding
	double BeachWidth = 50.;             	//Beach width landward of berm
	double BermHeight = 1.;               //Elevation of Beach Berm above junction (set beachwidth = 0 and bermheight = 0 for no beaches)
	int BeachType = 1;                    //0 = constant beach width, 1 = sinusoidal beach width, 2 = thinning beach
	double JunctionElevation = 0.;        //Elevation of the platform/cliff junction
	int RetreatType = 0;	                //Scenario of retreat 0 = constant, 1 = step change, 2 = gradual change
	int SteppedPlatformFlag = 0;          //Flag for a stepped platform (1 = steps, 0=no steps)
	double StepSize = 0.0;                //size of step (0=no steps)
	
	//Create Platform CRN object
	RockyCoastCRN RockyCoastCRNModel(RetreatRate1, RetreatRate2, RetreatType, ChangeTime, BeachWidth, BeachType, BermHeight, PlatformGradient, CliffHeight, JunctionElevation, Amp, SteppedPlatformFlag, StepSize);
	
	//Run the model
  //First for no beach
  cout << "\tBeach Width is 50 m variable" << endl;
  string OutFileName = "Variable_Bh1_Bw50.pdat";
	RockyCoastCRNModel.RunModel(OutFileName);
	
	//Create Platform CRN object
	BeachWidth = 50;
	BeachType = 0;
	RockyCoastCRN RockyCoastCRNModel2(RetreatRate1, RetreatRate2, RetreatType, ChangeTime, BeachWidth, BeachType, BermHeight, PlatformGradient, CliffHeight, JunctionElevation, Amp, SteppedPlatformFlag, StepSize);
	
	//Run the model
  //First for no beach
  cout << "\tBeach Width is 50 m" << endl;
  OutFileName = "BeachWidth_Bh1_Bw50.pdat";
	RockyCoastCRNModel2.RunModel(OutFileName);
	
	//Create Platform CRN object
	BeachWidth = 0;
	BermHeight = 0;
	RockyCoastCRN RockyCoastCRNModel3(RetreatRate1, RetreatRate2, RetreatType, ChangeTime, BeachWidth, BeachType, BermHeight, PlatformGradient, CliffHeight, JunctionElevation, Amp, SteppedPlatformFlag, StepSize);
	
	//Run the model
  //First for no beach
  cout << "\tBeach Width is 0 m" << endl;
  OutFileName = "BeachWidth_Bh0_Bw0.pdat";
	RockyCoastCRNModel3.RunModel(OutFileName);
	
	cout << "Done" << endl << endl;
}


