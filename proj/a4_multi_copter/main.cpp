//#####################################################################
// Main
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#include <iostream>
#include "MultiCopterDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
	////0-Test simulation (LV1, 3.1-3.2)
	////1-Test controller (LV2, 3.3-3.5)
	int flag = 1;

	MultiCopterDriver<3> driver;
	driver.Initialize(flag);
	driver.Run();
}

#endif
