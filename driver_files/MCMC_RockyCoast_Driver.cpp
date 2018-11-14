/*-----------------------------------------------------------------------------------------

	MCMC_RockyCoast_Driver.cpp
	
	Evolves a coastal profile in order to find the most likely history of coastal retreat
	in keeping with 10Be CRN measurements across a coastal platform

	Martin D. Hurst, British Geological Survey, January 2014

------------------------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "../MCMC_RockyCoast.hpp"
#include "../RockyCoastCRN.hpp"

using namespace std;

int main (int nNumberofArgs,char *argv[])
{
    //argv1 should be the input dataset
	//Test for correct input arguments
	if (nNumberofArgs!=5)
	{
		cout 	<< "FATAL ERROR: not enough inputs. The program needs:" << endl
					<< "\t1) the CRN input data filename" << endl
					<< "\t2) the paramfilename filename, and" << endl
					<< "\t3) the output chainfile " << endl
					<< "\t4) the number of iterations " << endl;
		exit(EXIT_SUCCESS);
	}

	// the name of the data
	char* DataFilename = argv[1];
	char* ParamFilename = argv[2];
    char* ChainFilename = argv[3];
    int NIterations = atoi(argv[4]);
  	
	// load an MCMC driver object
    MCMC_RockyCoast My_MCMC_Coast = MCMC_RockyCoast(DataFilename);
  
    //now run the metropolis algorithm along a chain
	My_MCMC_Coast.RunMetropolisChain(NIterations, ParamFilename, ChainFilename);
	
	cout << "Done!" << endl;
	
	return 0;
}

