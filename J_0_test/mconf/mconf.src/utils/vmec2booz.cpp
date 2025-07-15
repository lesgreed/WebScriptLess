#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "../include/CStconfig.h"

//*********************************************************************
//
int main(int argc, char* argv[])
{
  const char * Version = MConf::CStconfig::getVersion();
  std::string option = (argc>=2)?argv[1]:" ";

  if(argc<2||option=="--ver"||option=="-ver") {
    std::cout<<"Version "<<Version<<"  year.month"<<std::endl;
    return 1;
  }


  if(argc<2||option=="--help"||option=="-help"||option=="help"||option=="?"||option=="-?") {
    std::cerr<<"Version "<<Version<<std::endl;
    std::cerr<<"Convert vmec-wout file to boozer file. 2012, Turkin."<<std::endl;
    std::cerr<<"  Usage1: vmec2booz vmecFile  [boozFile [M N]]"<<std::endl;    
    std::cerr<<"   where vmecFile is the VMEC wout-file,"<<std::endl; 
    std::cerr<<"         boozFile is the output file with Boozer coordinates,"<<std::endl; 
    std::cerr<<"         M is the number of poloidal mods,"<<std::endl; 
    std::cerr<<"         N is the number of toroidal mods."<<std::endl; 
    std::cerr<<"   default values:"<<std::endl; 
    std::cerr<<"         M = 6*M_vmec"<<std::endl; 
    std::cerr<<"         N = 3*N_vmec"<<std::endl; 

    return 1;
  }

  int M = 0, N = 0;
  const char *fname   = argv[1];
  std::string fnameNew = "boozer.";
  fnameNew += fname;

  if(argc > 2) fnameNew = argv[2];
  if(argc > 3) M = atoi(argv[3]);
  if(argc > 4) N = atoi(argv[4]);

  MConf::CStconfig mConf;
  mConf.load (fname);
  if(!mConf.isOK()) {
    std::cerr<<fname<<": Loading error"<<std::endl; 
    return 2;
  }


  bool ok = mConf.writeVmec2Boozer(fnameNew.c_str(),0,M,N);
  if(!ok) std::cerr<<fnameNew<<": Converting error"<<std::endl;

  return ok?0:2;
}

