/*
  Calculate bootstrap current geometric factor for several magnetic configurations

based on code:

  MConf::C3dMesh mc("w7x-sc1.bc");
  std::cout <<"#r/a         Fbs"<<std::endl;  
  const int maxP = 33;
  for(int i=0; i<maxP; i++) {
    double x = 0.05+(0.975-0.05)*i/(maxP-1);    // x=sqrt(s)
    double fbs = mc.FbsDkes(x*x);
    std::cout<<x<<" "<<fbs<<std::endl;
  }

 */

//*********************************************************************

#include "../include/C3dMesh.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#ifdef WIN32
  #include <io.h>
  #include "XGetopt.h"
// #endif
#else
// #if (defined( __unix__ ) || (defined(linux)) || defined(__sun) )
  #define _T(x) x 
  #include <dirent.h>
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <unistd.h>
#endif

const double degree = 3.14159265358979323846/180;

//*********************************************************************
class Property  {
public:
  double xmin;
  double xmax;
  double dphi;  // in degree;
  double xlambdamax;
  double Btruncation;
  int lambdaLevels;
  int maxP;
  int turns;
  bool avoidResonances;
  bool accuracyTest;
  bool useRationalIota;
  int rationalIotaDenominator;

  Property()
  {
    maxP = 51;
    xmin = .05;
    xmax = .95;
    turns= 60;
    dphi = 1;  // in degree;
    avoidResonances=true;
    accuracyTest=true;
    useRationalIota=false;
    rationalIotaDenominator=150;
    xlambdamax=15;
    lambdaLevels=257;
    Btruncation = 1e-5;
  }
};

//*********************************************************************
//*********************************************************************
//*********************************************************************
class Device  {
  int rdef;
  Property prop;
  std::string option;
  std::string option3;
  std::string option4;
public:
  std::string FbsName;
  std::string name; // config. name trancated 
  std::string nameOriginal; 
  std::vector<double> x;
  std::vector<double> Fbs;
  std::vector<double> FbsN;
  std::vector<double> f_trap;
  std::vector<double> eps_eff;
  std::vector<MConf::Vector3d> Fbs_all;
  std::vector<MConf::Vector3d> FbsN_all;

  Device(Property &parm, std::string opt, std::string opt3,std::string opt4) {
    prop = parm;
    option = opt.substr(0,1);
    option3= opt3.substr(0,1);
    option4= opt4.substr(0,1);
    rdef = 1;
    if(option=="d") rdef = 1;
    else if(option=="v") rdef = 2;
    else if(option=="g") rdef = 3;
    else if(option=="a") rdef = 4;
  }
  
  int points() { return int(x.size());  } 

  bool load(std::string &path,std::string &fname) {
    x.clear();
    f_trap.clear();
    Fbs.clear();
    FbsN.clear();
    Fbs_all.clear();
    FbsN_all.clear();
    MConf::C3dMesh mc;
    std::string pathname = path+fname;
    if(false==mc.load(pathname.c_str())) {
      std::cerr << "Error loading "<<pathname<<std::endl;
      return false;
    }
    mc.truncate(prop.Btruncation); // truncate Bmn spectrum at level Bmn/B00 = 1e-6
    name = fname;
    nameOriginal = fname;

    std::cerr<<"Processing file: "<<pathname<<std::endl;
    clock_t tStart = clock();  //timing
    //std::cerr<<"   Calculating bootstrap current geometric factor F_bs..."<<std::endl;
    mc.FbsInvalidate();
    mc.FbsSetXeffParam   (prop.xmin, prop.xmax, prop.maxP);
    mc.FbsSetTracingParam(prop.turns, prop.accuracyTest, prop.dphi*degree);
    mc.FbsSetIotaParam   (prop.avoidResonances, prop.useRationalIota, prop.rationalIotaDenominator); 
    mc.FbsSetMagnMomentParam (prop.lambdaLevels, prop.xlambdamax);
    mc.FbsCreate();
    std::cerr<<"   F_bs done, time/s="<<double(clock() - tStart)/CLOCKS_PER_SEC<<std::endl;

    int len=11;
    if(option3=="b") {
      len=22;
    }
    else if(option3=="n") {
      len=11;
    }
    else {
      len=11;
    }

    if(option=="a") len *= 3;
    if(option4=="p") len += 11;
    if(len>32) 
    {
      name.resize(len-1,' ');  // adjust name length for printing
      //name.resize(len,' ');  // adjust name length for printing
    }

    for(int i=0; i<prop.maxP; i++) {
      double x1 = prop.xmin+(prop.xmax-prop.xmin)*i/double(prop.maxP-1); 
      x.push_back(x1);

      double s = x1*x1;
      f_trap.push_back(mc.ftrapped(s));

      double ev  = mc.r(s)/mc.R0();
      double ed  = mc.rdkes(s)/mc.R0();
      double eg  = mc.rgraz(s)/mc.R0();
      double eav = mc.r(1)/mc.R0();    // (a/R)_vmec
      double ead = mc.rdkes(1)/mc.R0();
      double eag = mc.rgraz(1)/mc.R0();
      double iota = fabs(mc.iota(s));
      double ftraptok = 1.469*sqrt(ev)*(1-0.325*ev-0.034*ev*ev);
      double fbstok = ftraptok*(1+sqrt(1-eav*eav))/(2*iota*ev);

      double fbstokNormv = 1.46/(iota*sqrt(ev)); // ftrap/(iota*e)
      double fbstokNormd = 1.46/(iota*sqrt(ed));
      double fbstokNormg = 1.46/(iota*sqrt(eg));

      ////double fbstokNormv = (1+sqrt(1-eav*eav))/2*fabs( 1.46/(iota*sqrt(ev)) );
      ////double fbstokNormd = (1+sqrt(1-ead*ead))/2*fabs( 1.46/(iota*sqrt(ed)) );
      ////double fbstokNormg = (1+sqrt(1-eag*eag))/2*fabs( 1.46/(iota*sqrt(eg)) );

      double lambdab = mc.GrazLambdab(s);
      double fbsVmec = mc.FbsVmec(s);
      double fbsDkes = mc.FbsDkes(s);
      
      double lambdabN = lambdab/fbstokNormg;
      double fbsVmecN = fbsVmec/fbstokNormv;
      double fbsDkesN = fbsDkes/fbstokNormd;

      MConf::Vector3d e, eN;
      switch(rdef) {
        case 1:  //r_dkes
          FbsName = "Fbs_dkes   ";
          if(option3=="b") {
            FbsName+= "FbsN_dkes  ";
          }
          else if(option3=="n") {
            FbsName = "FbsN_dkes  ";
          }
          Fbs.push_back(fbsDkes);
          FbsN.push_back(fbsDkesN);
          break;
        case 2:  //r_vmec
          FbsName = "Fbs_vmec   ";
          if(option3=="b") {
            FbsName+= "FbsN_vmec  ";
          }
          else if(option3=="n") {
            FbsName = "FbsN_vmec  ";
          }
          Fbs.push_back(fbsVmec);
          FbsN.push_back(fbsVmecN);
          break;
        case 3:  //lambda_b
          FbsName = "lambda_b   ";
          FbsName+= "lambda_bN  ";
          Fbs.push_back(lambdab);
          FbsN.push_back(lambdabN);
          break;
        case 4:  // all
          FbsName = "Fbs_dkes   Fbs_vmec   lambda_b   ";
          if(option3=="b") {
            FbsName+= "FbsN_dkes  FbsN_vmec  lambda_bN  ";
          }
          else if(option3=="n") {
            FbsName = "FbsN_dkes  FbsN_vmec  lambda_bN  ";
          }
          e.set (fbsDkes, fbsVmec, lambdab);
          eN.set(fbsDkesN,fbsVmecN,lambdabN);
          Fbs_all.push_back(e);
          FbsN_all.push_back(eN);
          break;
        default:
          FbsName = "Fbs_dkes   ";
          if(option3=="b") {
            FbsName+= "FbsN_dkes  ";
          }
          else if(option3=="n") {
            FbsName = "FbsN_dkes  ";
          }
          Fbs.push_back(fbsDkes);
          FbsN.push_back(fbsDkesN);
          break;
      }
      if(rdef!=4&&option3=="b") {  
        FbsName = " " + FbsName;
        FbsName += " ";
      }
      if(rdef!=4) {  
        FbsName += " ";
      }
      if(option4=="p") {
        if(option3=="b")
          FbsName += "f_trap    ";
        else
          FbsName += "f_trap     ";
      }           
    }
    return true;
  } 
};

static bool createList(const std::string &fname, std::vector<std::string> &list);
static void createFileList(std::string dir, std::string mask, std::vector<std::string> &list);
static std::string getPath(std::string fileName);
static const char *fl_filename_name(const char *name);

template <class T> void swap (T &v1, T &v2) {T w = v1;  v1 = v2;  v2 = w;};

//*********************************************************************
//  
//
int main(int argc, char *argv[])
{
  Property prop;
  Property propFastest;
  Property propFast;
  Property propStandard;
  Property propDetailed;

  //case FAST: see constructor
  propFastest  = propFast;
  propStandard = propFast;
  propDetailed = propFast;
  { //case FASTEST:   for testing
    Property &prp = propFastest;
    prp.maxP = 11;
    prp.turns= 10;
    prp.Btruncation = 1e-1;
  }  
  { //case STANDARD:
    Property &prp = propStandard;
    prp.maxP = 101;
    prp.turns= 100;
    prp.lambdaLevels=257;
  }
  { //case DETAILED:
    Property &prp = propDetailed;
    prp.maxP = 201;
    prp.turns= 200;
    prp.dphi = 0.5; // in degree;
    prp.avoidResonances=false;
    prp.accuracyTest=false;
    prp.lambdaLevels=513;
  }
  prop = propStandard; 

  std::string option = (argc>=2)?argv[1]:" ";

  const char * Version = MConf::CStconfig::getVersion();

  if(argc<2||option=="--ver"||option=="-ver") {
    std::cout<<"Version "<<Version<<"  year.month"<<std::endl;
    return 1;
  }


  if(argc<2||option=="--help"||option=="-help"||option=="help"||option=="?"||option=="-?") {
    std::cerr<<"Version "<<Version<<std::endl;
    std::cerr<<"Calculate bootstrap current geometric factor. 2011, Turkin."<<std::endl;
    std::cerr<<"The program is based on papers by"<<std::endl;
    
std::cerr<<"P.Helander et al,Phys.Plasmas 18,092505(2011)."<<std::endl;
std::cerr<<"N.Nakajima et al,Nuclear Fusion,29,605(1989)."<<std::endl;
std::cerr<<"V.Nemov,V.Kalyuzhnyj,S.Kasilov et al,PPCF,46,179(2004)."<<std::endl<<std::endl;
    
    std::cerr<<"Usage1:  Fbs [option1] [option2] [option3] [-p] [-m] filename"<<std::endl;
    std::cerr<<"   where filename is the directory or path, and the file name,"<<std::endl; 
    std::cerr<<"   which can include wildcard characters, for example,"<<std::endl; 
    std::cerr<<"   an asterisk or a question mark. The files must be in"<<std::endl; 
    std::cerr<<"   format of VMEC-wout-file, bc-file or new_boz-file."<<std::endl<<std::endl;
    std::cerr<<"Usage2:  Fbs [option1] [option2] [option3] [-p] [-m] @filelist"<<std::endl;
    std::cerr<<"   where filelist is the name of a file that contains list of files"<<std::endl;
    std::cerr<<"   for which the calculations will be made "<<std::endl<<std::endl;;

    std::cerr<<"Options:"<<std::endl;
    std::cerr<<"  -m      produce minimal output if this option is supplied"<<std::endl;
//    std::cerr<<"  -e      output eps_eff(vmec)"<<std::endl;
    std::cerr<<"  -p      output fraction of trapped particles"<<std::endl;
//    std::cerr<<"  -t   output average gradB drift velocity of trapped particles"<<std::endl;
    std::cerr<<"option1:"<<std::endl;
    std::cerr<<"  -d   default option,"<<std::endl;
    std::cerr<<"       use DKES definition of the minor radius r=sqrt(F(s)/(pi*B_00(s)))"<<std::endl;
    std::cerr<<"  -v   use VMEC definition of the minor radius r=a*sqrt(s)"<<std::endl;
    std::cerr<<"  -g   calculate lambda_b with Graz definition of the minor radius dr=dV(s)/A(s)=ds/<|grad(s)|>"<<std::endl;
    std::cerr<<"  -a   output Fbs for all definitions of the minor radius"<<std::endl;
    std::cerr<<" where F is the toroidal flux,"<<std::endl; 
    std::cerr<<"       s is the toroidal flux normalized to the value"<<std::endl; 
    std::cerr<<"         at the last closed magnetic surface,"<<std::endl;
    std::cerr<<"       a is the plasma radius,"<<std::endl;
    std::cerr<<"       A is the area of the flux surface."<<std::endl; 
    std::cerr<<"option2 controls the calculation accuracy:"<<std::endl;
    std::cerr<<"  -s      default option,"<<std::endl;
    std::cerr<<"          use standard set of parameters,"<<std::endl;
    std::cerr<<"  -f      use set of parameters for fast calculation"<<std::endl;
    std::cerr<<"  -c      use set of parameters for careful calculation"<<std::endl;
    std::cerr<<"option3 controls the output:"<<std::endl;
    std::cerr<<"  -n      output values normalized to the tokamak values"<<std::endl;
    std::cerr<<"  -b      output both values: non-normalized and normalized to the tokamak values"<<std::endl<<std::endl;
    std::cerr<<"Note: Only the first letter of the option must be supplied."<<std::endl;
    std::cerr<<"   Option letters may be combined, e.g., '-vcn' is equivalent to '-v -c -n'."<<std::endl;
    std::cerr<<"   Option letters are case sensitive."<<std::endl;

    std::cerr<<"Examples:"<<std::endl;
    std::cerr<<"        Fbs  -ver"<<std::endl;
    std::cerr<<"        Fbs  -help"<<std::endl;
    std::cerr<<"        Fbs  w7x-sc1.bc"<<std::endl;
    std::cerr<<"        Fbs -d -c w7x-sc1.bc"<<std::endl;
    std::cerr<<"        Fbs -a w7x-sc1.bc"<<std::endl;
    std::cerr<<"        Fbs  path"<<std::endl;
    std::cerr<<"        Fbs  \"path/wout*\""<<std::endl;
    std::cerr<<"        Fbs  \"path/*.bc\""<<std::endl;
    std::cerr<<"        Fbs  \"path/pattern\""<<std::endl;
    std::cerr<<"        Fbs  \"path/pattern\" > outputfile"<<std::endl;
    std::cerr<<"        Fbs  -vns  @w7xlist > w7xFbs.txt"<<std::endl<<std::endl;
    return 1;
  }
   
  std::string option1="d"; 
  std::string option2="s"; 
  std::string option3=" ";  
  std::string option4=" ";  
  std::string option5=" ";  
  std::string option_m=" "; 
  int c;
  while ((c = getopt(argc, argv, _T("dvgasfc/bmnp"))) != -1) {
    switch (c) {
      case _T('d'):
      case _T('v'):
      case _T('g'):
      case _T('a'):
        option1 = c;
        break;
      case _T('s'):
      case _T('f'):
      case _T('c'):
      case _T('/'):
        option2 = c;
        break;
      case _T('n'):
      case _T('b'):
        option3 = c;
        break;
      case _T('m'):
        option_m = c;
        break;
      case _T('p'):
        option4 = c;
        break;
      case _T('t'):
        option5 = c;
        break;
      case _T('?'):
#ifdef WIN32        
        std::cerr<<_T("ERROR: invalid option ")<<argv[optind-1]<<std::endl;
#endif        
        break;
      default:
        std::cerr<<_T("WARNING: no handler for option ")<<c<<std::endl;
        break;
    }
  }

  std::string pathname = "";
  
  if (optind < argc) {  // check for non-option args here
    pathname = argv[optind];
  }
  else {
    std::cerr<<"No files"<<std::endl;
    return 1;
  }
    
  const char *name = fl_filename_name(pathname.c_str());
  std::string mask = name[0]!=0?name:"*";
  std::string path = getPath(pathname);

  std::vector<std::string> filelist;
  
  if(pathname[0]=='@') {
    path="";
    createList(pathname.c_str()+1, filelist);
  }
  else {
    createFileList(path,mask,filelist);
  }
   
  if(filelist.empty()) {
    std::cerr<<"No files"<<std::endl;
    return 1;
  }
  
#ifdef _DEBUG
  std::cerr<<"option1 : "<<option1<<std::endl;
  std::cerr<<"option2 : "<<option2<<std::endl;
  std::cerr<<"option3 : "<<option3<<std::endl;
  std::cerr<<"option_m: "<<option_m<<std::endl;
  std::cerr<<path<<std::endl;
  std::cerr<<mask<<std::endl;
  std::cerr<<"******************************"<<std::endl;
  // for (unsigned int i=0;i<filelist.size();i++) {
  //   std::cerr<<filelist[i]<<std::endl;
  // }
  // std::cerr<<"-------------------------------------"<<std::endl;
  for (unsigned int i=0;i<filelist.size();i++) {
    FILE *fp = fopen((path+filelist[i]).c_str(),"r");
    std::cerr<<(fp?"opened: ":"error opening: ")<<path+filelist[i]<<std::endl;
    if(fp) fclose(fp);
  }
  return 0;
#endif

  if(option2=="/") prop = propFastest;
  if(option2=="f") prop = propFast;
  if(option2=="s") prop = propStandard;
  if(option2=="c") prop = propDetailed;

  if(option_m!="m") {
    std::cerr<<std::endl;
    std::cerr<<std::setw(8)<<prop.maxP         <<"  - number of r_eff-grid points"<< std::endl;
    std::cerr<<std::setw(8)<<prop.xmin         <<"  - sqrt(s)_min"<< std::endl;
    std::cerr<<std::setw(8)<<prop.xmax         <<"  - sqrt(s)_max"<< std::endl;
    std::cerr<<std::setw(8)<<prop.turns        <<"  - max. number of toroidal turns"<< std::endl;
    std::cerr<<std::setw(8)<<prop.dphi         <<"  - field line tracing step in degree"<< std::endl;  // in degree;
    std::cerr<<std::setw(8)<<prop.Btruncation  <<"  - Tolerance B_mn/B_00"<< std::endl;
    std::cerr<<std::setw(8)<<prop.lambdaLevels <<"  - number of integration points over magn.moment lambda"<< std::endl;
    std::cerr<<std::setw(8)<<prop.xlambdamax   <<"  - y_max, lambda_max=tanh(y_max)"<< std::endl;
    
    std::cerr<<std::setw(8)<<(prop.accuracyTest?"yes":"no")    <<"  - use accuracy test"<< std::endl;
    std::cerr<<std::setw(8)<<(prop.avoidResonances?"yes":"no") <<"  - try to avoid resonances"<< std::endl;
    if(prop.useRationalIota)
      std::cerr<<std::setw(8)<<(prop.useRationalIota?"yes":"no") <<"  - use rational approximation to iota"<< std::endl;
    if(prop.useRationalIota)
      std::cerr<<std::setw(8)<<prop.rationalIotaDenominator      <<"  - rational iota denominator close to this value"<< std::endl;
  }

  Device wrk(prop,option1,option3,option4);
  std::vector<Device> dev;

  for (unsigned int i=0;i<filelist.size();i++) {
    if(wrk.load(path,filelist[i])) 
      dev.push_back(wrk);
  }

  if(dev.empty()) {
    std::cerr<<"No valid files"<<std::endl;
    return 1;
  }
  else if(dev.size()==1) dev[0].name = dev[0].nameOriginal;

  if(option_m!="m") {
    std::cout <<"% Bootstrap current geometric factor."<<std::endl;
    std::cout <<"% Fbs_vmec is calculated using VMEC definition of the minor radius r=a*sqrt(s),"<<std::endl;
    std::cout <<"% Fbs_dkes is calculated using DKES definition of the minor radius r=sqrt(F(s)/(pi*B_00(s))),"<<std::endl;
    std::cout <<"% lambda_b is calculated using Graz definition of BS factor and the minor radius dr=ds/<|grad(s)|>,"<<std::endl;
    std::cout <<"% FbsN_dkes, FbsN_vmec, lambda_bN are the above values normalized to the tokamak value 1.46sqrt(e)/(e*iotabar), where e=r_eff/R"<<std::endl;
    std::cout <<"% F_bs = f_t*G_bs*dr/dV, where G_bs is defined in N.Nakajima et al,NF, vol.29, 605(1989), http://dx.doi.org/10.1088/0029-5515/29/4/006" <<std::endl;
    std::cout <<"% lambda_b is defined in V.Nemov et al, PPCF, vol.46, 179(2004), http://dx.doi.org/10.1088/0741-3335/46/1/011" <<std::endl;
    std::cout <<"% s is the normalized toroidal flux."<<std::endl;
  }

  // print header
  if(option_m!="m") {
    std::cout<<"% "<<dev.size()<<" - number of devices (columns)"<<std::endl;
    std::cout<<"% "<<dev[0].points()<<" - number of radial points (rows)"<<std::endl;
    std::cout<<"%xaxis: r/a"<<std::endl;
    std::cout<<"%yaxis: F_bs"<<std::endl;
    std::cout<<"%names:   ";
    for(unsigned int k=0; k<dev.size(); k++) {   // print device names
      std::cout<<" "<<dev[k].name;
    }
    std::cout<<std::endl;
  }
  
  std::cout<<"% version "<<Version<<std::endl;
  std::cout<<"% sqrt(s)  ";
  for(unsigned int k=0; k<dev.size(); k++)  {    // print eps names
    std::cout<<dev[k].FbsName;
  }
  std::cout<<std::endl;

  // print data
#ifdef _MSC_VER
  std::cout.precision(3);
#else
  std::cout.precision(4);
#endif
  std::cout.setf(std::ios::scientific);
  for(int i=0; i<dev[0].points(); i++) {
    std::cout<<dev[0].x[i];
    for(unsigned int k=0; k<dev.size(); k++) {
      if(option1=="a") {
        if(option3=="b") {
          std::cout<<" "<<dev[k].Fbs_all[i]<<" "<<dev[k].FbsN_all[i];
        }
        else if(option3=="n") {
          std::cout<<" "<<dev[k].FbsN_all[i];
        }
        else {
          std::cout<<" "<<dev[k].Fbs_all[i];
        }
      }
      else {
        if(option3=="b") {
          std::cout<<" "<<dev[k].Fbs[i]<<" "<<dev[k].FbsN[i];
        }
        else if(option3=="n") {
          std::cout<<" "<<dev[k].FbsN[i];
        }
        else {
          std::cout<<" "<<dev[k].Fbs[i];
        }
      }
      if(option4=="p") {
        std::cout<<" "<<dev[k].f_trap[i];
      }
    }
    std::cout<<std::endl;
  }

  return 0;
}

//*********************************************************************
//*********************************************************************
//*********************************************************************
// "$Id: filename_match.cxx 5190 2006-06-09 16:16:34Z mike $"
//
// Pattern matching routines for the Fast Light Tool Kit (FLTK).
//
// Copyright 1998-2005 by Bill Spitzak and others.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA.
//
// Please report all bugs and problems on the following page:
//
//     http://www.fltk.org/str.php
//

/* Adapted from Rich Salz. */
#include <ctype.h>

static int fl_filename_match(const char *s, const char *p) 
{
  int matched(0);

  for (;;) {
    switch(*p++) {

    case '?' :	// match any single character
      if (!*s++) return 0;
      break;

    case '*' :	// match 0-n of any characters
      if (!*p) return 1; // do trailing * quickly
      while (!fl_filename_match(s, p)) if (!*s++) return 0;
      return 1;

    case '[': {	// match one character in set of form [abc-d] or [^a-b]
      if (!*s) return 0;
      int reverse = (*p=='^' || *p=='!'); if (reverse) p++;
      matched = 0;
      char last = 0;
      while (*p) {
        if (*p=='-' && last) {
          if (*s <= *++p && *s >= last ) matched = 1;
          last = 0;
        } else {
          if (*s == *p) matched = 1;
        }
        last = *p++;
        if (*p==']') break;
      }
      if (matched == reverse) return 0;
      s++; p++;}
    break;

    case '{' : // {pattern1|pattern2|pattern3}
    NEXTCASE:
    if (fl_filename_match(s,p)) return 1;
    for (matched = 0;;) {
      switch (*p++) {
      case '\\': if (*p) p++; break;
      case '{': matched++; break;
      case '}': if (!matched--) return 0; break;
      case '|': case ',': if (matched==0) goto NEXTCASE;
      case 0: return 0;
      }
    }
    case '|':	// skip rest of |pattern|pattern} when called recursively
    case ',':
      for (matched = 0; *p && matched >= 0;) {
        switch (*p++) {
        case '\\': if (*p) p++; break;
        case '{': matched++; break;
        case '}': matched--; break;
        }
      }
      break;
    case '}':
      break;

    case 0:	// end of pattern
      return (!*s)?1:0;

    case '\\':	// quote next character
      if (*p) p++;
    default:
      if (tolower(*s) != tolower(*(p-1))) return 0;
      s++;
      break;
    }
  }
}

// end "$Id: filename_match.cxx 5190 2006-06-09 16:16:34Z mike $"


// returns pointer to the filename, or pointer to null if name ends with /
static const char *fl_filename_name(const char *name) 
{
  const char *p,*q;
  if (!name) return (0);
  q = name;
  if (q[0] && q[1]==':') q += 2; // skip leading drive letter
  for (p = q; *p; p++) if (*p == '/' || *p == '\\') q = p+1;
  return q;
}

// get path from the full filename
static std::string getPath(std::string fileName)
{
  const int szfn = 2048;
  char *fn = new char[szfn];
  //if( fl_filename_absolute(fn,szfn,fileName.c_str()) )  fileName = fn;
  const char * c_path = fileName.c_str();
  size_t i = fl_filename_name(c_path) - c_path;
  std::string path = i==0?"./":fileName.substr(0,i);
  delete fn;
  return path;
}

#if (defined( __unix__ ) || (defined(linux)) || defined(__sun) )
static mode_t getMode(const std::string& path)
{
  struct stat _stat;
  if ( 0 == lstat(path.c_str(), &_stat) )  {
    return _stat.st_mode;
  }
  return 0;
}
    
static bool isDirectory(const std::string& path)
{
  return S_ISDIR(getMode(path));
}
#endif

//*******************************************************************************
// function creates 'list' of files belonging to directory 'dir' and satisfying to the 'mask'
static void createFileList(std::string dir, std::string mask, std::vector<std::string> &list)
{
  list.clear();

#ifdef WIN32
  if     (dir[dir.size()-1 ]=='/')   dir += "*";
  else if(dir[dir.size()-1 ]=='\\')  dir += "*";
  struct _finddata_t file;
  intptr_t hFile =_findfirst(dir.c_str(), &file );
  if(hFile==-1) return;
  do {
    if( (file.attrib&_A_SUBDIR)==0 ) {
      std::string name = file.name;
      if(!(name=="."||name=="..") ) {
        if(fl_filename_match(name.c_str(), mask.c_str())) 
          list.push_back(name);
      }
    }
  } while(_findnext(hFile,&file)==0);
  _findclose( hFile );
#elif (defined( __unix__ ) || (defined(linux)) || defined(__sun) )
  DIR *dp;
  if((dp=opendir(dir.c_str()))==0) {
    std::cerr<<"Error("<<errno<<") opening "<<dir<<std::endl;
    return;
  }
  struct dirent *dirp;
  while( (dirp=readdir(dp))!=0 ) {
    std::string name = dirp->d_name;
    std::string pathname = dir + name;
    if(isDirectory(pathname)) continue;
    if( fl_filename_match(name.c_str(), mask.c_str()) ) list.push_back(name);
  }
  closedir(dp);
#endif
  sort (list.begin(), list.end());
}

//*******************************************************************************
// function creates 'list' of files reading the names from the file 'fname'
static bool createList(const std::string &fname, std::vector<std::string> &list)
{
  list.clear();
  std::string str;
  std::ifstream in(fname.c_str());
  if(!in.is_open()) return false;
  std::getline(in,str);   
  while(in) {
    list.push_back(str);
    std::getline(in,str);
  }
  return true;
}


#if 0

class mcout {
public:
  std::ostream *myout;
  mcout() {myout = &std::cout;};
  void Redirect(std::ostream &out) {myout = &out;};
};

mcout Mcout;
#define mout (*Mcout.myout)

std::ostream *savedcout;
savedcout = &std::cout;
*savedcout<<"jgjgjgjgjgjgg"<<std::endl;  

mout << "57572573257537532"<<std::endl;

return 1;

#endif

