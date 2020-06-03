//#####################################################################
// Main
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#include <iostream>
#include <cstdio>
#include "MassSpringInteractiveDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
	std::string output_dir="output";
	int test=1;
	for(int i=0;i<argc;i++){
		if(strcmp(argv[i],"-o")==0){
			output_dir=std::string(argv[++i]);}
		else if(strcmp(argv[i],"-test")==0){
			test=atoi(argv[++i]);}}
	std::cout<<"[Mass spring simulation driver arguments]: -o ="<<output_dir<<", -test ="<<test<<std::endl;

	MassSpringInteractivDriver<3> driver;
	driver.scale=1;
	driver.test=test;

	driver.Initialize();
	driver.Run();	
}

#endif