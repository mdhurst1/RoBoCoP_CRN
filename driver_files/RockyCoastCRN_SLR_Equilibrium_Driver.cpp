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
#include <omp.h>
#include <unistd.h>
#include "../RockyCoastCRN.hpp"

using namespace std;

int main()
{
	//Input parameters
	double RetreatRate1 = 0.325;          //Retreat Rate (m/yr) at the start of the model run
	double RetreatRate2 = 0.1;            //Retreat Rate (m/yr) at the end of the model run
	double ChangeTime = 0;                //Time to change retreat rates if a step change
	double PlatformGradient = 0.002;    	//Platform gradient (average if stepped platform)
	double Amp = 1.;                     	//Tidal amplitude (1/2 tidal range)
	double SLR = 0.0005;											//Rate of sea level rise (m/yr)
	double CliffHeight = 0.;             //Cliff height for shielding
	double BeachWidth = 0.;               //Beach width (currently assumes total shielding)
	double BermHeight = 0.;               //Elevation of Beach Berm above junction (set beachwidth = 0 and bermheight = 0 for no beaches)
	int BeachType = 0;                    //0 = constant beach width, 1 = sinusoidal beach width, 2 = thinning beach
	double ElevInit = 0.5;                //Elevation of the platform/cliff junction
	int RetreatType = 0;	                //Scenario of retreat 0 = constant, 1 = step change, 2 = gradual change
	int SteppedPlatformFlag = 0;          //Flag for a stepped platform (1 = steps, 0=no steps)
	double StepSize = 0.0;                //size of step (0=no steps)

	string OutputFileName = "Equilibrium_R0.32_SLR0.5.xzn";
	
	//Declare model
	RockyCoastCRN RockyCoastCRNModel(RetreatRate1, RetreatRate2, RetreatType, ChangeTime,	
											BeachWidth, BeachType, BermHeight, PlatformGradient, CliffHeight, 		
											ElevInit, SLR, Amp, SteppedPlatformFlag, StepSize);
											
	//Run the model
	RockyCoastCRNModel.RunModel(OutputFileName);
	
}
