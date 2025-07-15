#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "../include/CStconfig.h"

using namespace std;

//*********************************************************************
//
int main(int argc, char* argv[])
{
  const char * Version = MConf::CStconfig::getVersion();
  std::string option = (argc>=2)?argv[1]:" ";

  if(argc<2||option=="--ver"||option=="-ver") {
    std::cout<<"MConf version "<<Version<<"  year.month"<<std::endl;
    return 1;
  }


  if(argc<2||option=="--help"||option=="-help"||option=="help"||option=="?"||option=="-?") {
    std::cerr<<"MConf version "<<Version<<std::endl;
    std::cerr<<"Calculate S matrix for magnetic configuration file."<<std::endl;
    std::cerr <<" see S matrix in P.I.Strand,W.A.Houlberg Physics of Plasmas 8,2782(2001)"<<std::endl;
    std::cerr<<"  https://doi.org/10.1063/1.1366618 " <<std::endl;
    std::cerr<<"  Usage: smatrix infile [outfile]"<<std::endl;    
    std::cerr<<"    where infile is the VMEC wout-file or Boozer file."<<std::endl; 
    std::cerr<<"         outfile is the output file with the S matrix components"<<std::endl; 

    return 1;
  }

  int M = 0, N = 0;
  const char *fname   = argv[1];
  std::string fnameNew = "Sij.";
  fnameNew += fname;

  if(argc > 2) fnameNew = argv[2];
  //  if(argc > 3) M = atoi(argv[3]);
  //  if(argc > 4) N = atoi(argv[4]);

  MConf::CStconfig magConf;
  magConf.load (fname);
  if(!magConf.isOK()) {
    std::cerr<<fname<<": Loading error"<<std::endl; 
    return 2;
  }
  
  double  truncSav = magConf.truncate(1e-8);
  magConf.setAccuracy(1e-6,1e-8);
  double s1st = magConf.get1st_s(); //get first point on s, which is stored in a bc-file
  double sLst = magConf.getLast_s();
  s1st = s1st<0.0002?0.0002:s1st;
  double x1st = sqrt(s1st);
  double xLst = sqrt(sLst);
  x1st = 0.01;
  xLst = 0.999;

  double FluxMax = magConf.Flux(1);
  double R0 = magConf.R0();
  double a = magConf.r(1.);

  int maxPnt = magConf.nSurfaces();
  if(maxPnt>50) maxPnt=50; 

  std::ofstream out(fnameNew.c_str(),std::ios::out);
  if(!out.is_open()) return 1;
  
  std::string confname = magConf.fname();
  out <<"# Sij for "<<confname<<std::endl;
  out <<"# see Sij in:"<<std::endl;
  out <<"#   Magnetic flux evolution in highly shaped plasma"<<std::endl; 
  out <<"#   P.I.Strand,W.A.Houlberg Physics of Plasmas 8,2782(2001)"<<std::endl;
  out <<"#   https://doi.org/10.1063/1.1366618"<<std::endl;
  out <<"# Sij(s) is for coordinates with flux surface label s=Flux/FluxMax, FluxMax="<<FluxMax<<std::endl;
  out <<"# Sij(r) is Sij when flux surface label is r=a*sqrt(s), a="<<a<<std::endl;
  out <<"# Sij(r) = Sij(s)*dr/ds=Sij(s)*r/(2*s)"<<std::endl;
  out <<"# Sij(Flux) is Sij when flux surface label is Flux (toroidal flux)"<<std::endl;
  out <<"# Sij(Flux) = Sij(s)*FluxMax, FluxMax="<<FluxMax<<std::endl;
  out <<"# delta_iota = mu0*Itor(s)/(Sij(s)*FluxMax)"<<std::endl;
  out <<"# rows="<<maxPnt<<std::endl;
  out <<"# s           iota          iota_CF     delta_iota     S11(s)       S12(s)       S22(s)       mu0*Itor     r            S11(r)       S12(r)       S22(r)"<<std::endl;
  out.precision(4);
  int w=13;
  out.setf(std::ios::scientific);

  for(int i=0; i<maxPnt; i++) {
    double r = x1st+i*(xLst-x1st)/(maxPnt-1);
    double s = r*r;
    r = magConf.r(s);

    double SF11 = magConf.S11(s);  //S-matrix in (Flux,th,ph)-coordinates
    double SF12 = magConf.S12(s);
    double SF21 = magConf.S21(s);
    double SF22 = magConf.S22(s);

    const double pi  = 3.14159265358979323846;
    const double mu0 = 4e-7*pi;
    double mu0It = mu0*magConf.It(s);    
    double It = magConf.It(s);     
    double Ip = magConf.Ip(s);  

    double iCF  = -SF12/SF11;  // current free iota
    double dIota = mu0*It/SF11;
    double iota  = magConf.iota(s);
    double iota1 = dIota+iCF;  

    out<<std::setw(w-1)<<s<<std::setw(w)  //s
      <<iota<<std::setw(w)   // iota 
      <<iCF<<std::setw(w)    // iota_CF;
      <<dIota<<std::setw(w)  // delta iota
      <<SF11/FluxMax<<std::setw(w)    // s11(s)  is the S-matrix in (s,th,ph)-coordinates
      <<SF12/FluxMax<<std::setw(w)    // s12
      <<SF22/FluxMax<<std::setw(w)    // s22
      <<mu0It<< std::setw(w)          // mu0It
      <<r<< std::setw(w)              // a*sqrt(s)
      <<SF11/FluxMax*r/(2*s)<< std::setw(w)    // sr11(r) sr11 is the S-matrix in (r,th,ph)-coordinates r=a*sqrt(s)
      <<SF12/FluxMax*r/(2*s)<< std::setw(w)    // sr12
      <<SF22/FluxMax*r/(2*s)                   // sr22
      <<std::endl;
  }
    return 0;
}



