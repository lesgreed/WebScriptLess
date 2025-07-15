
/*
  calculate effective helical ripple for several magnetic configurations

based on code:

  MConf::C3dMesh mc("w7x-sc1.bc");
  std::cout <<"#r/a         eps"<<std::endl;  
  const int maxP = 33;
  for(int i=0; i<maxP; i++) {
    double x = 0.05+(0.975-0.05)*i/(maxP-1);    // x=sqrt(s)
    double eps = mc.epsEffDkes(x*x);
    std::cout<<x<<" "<<eps<<std::endl;
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


//*********************************************************************
class Device  {
  int rdef;
  int maxP;
  double trunc;
  std::string option;
  std::string option2;
  std::string option3;
  std::string meshType;
public:
  std::string epsName;
  std::string name; // config. name truncated 
  std::string nameOriginal; 
  std::vector<double> x;
  std::vector<double> f_trap;
  std::vector<double> eps_eff;
  std::vector<MConf::Vector3d> eps_eff_all;
  std::vector<MConf::Vector3d> Gamma_factor;
  Device(int maxPoints, double truncation, std::string opt, std::string opt2, std::string opt3, std::string meshTyp) {
    option = opt.substr(0,1);
    option2= opt2.substr(0,1);
    option3= opt3.substr(0,1);
    meshType= meshTyp;
    rdef = 1;
    maxP = maxPoints;
    trunc = truncation;
    if(option=="d") rdef = 1;
    else if(option=="v") rdef = 2;
    else if(option=="g") rdef = 3;
    else if(option=="a") rdef = 4;
  }
  size_t points() {
    return x.size();
  } 
  bool load(std::string &path,std::string &fname) {
    x.clear();
    f_trap.clear();
    eps_eff.clear();
    eps_eff_all.clear();
    Gamma_factor.clear();
    MConf::C3dMesh mc;
    std::string pathname = path+fname;
    if(false==mc.load(pathname.c_str())) {
      std::cerr << "Error loading "<<pathname<<std::endl;
      return false;
    }
    std::cerr<<"Processing file: "<<pathname<<std::endl;
    clock_t tStart = clock();  //timing
    mc.truncate(trunc); // truncate Bmn spectrum at level Bmn/B00 = 1e-6
    name = fname;
    nameOriginal = fname;

    int len=11;
    if(option=="a") len *= 3;
    if(option2=="t") len += 33;
    if(option3=="p") len += 11;
    if(len>32) 
    {
      name.resize(len-1,' ');  // adjust name length for printing
      //name.resize(len,' ');  // adjust name length for printing
    }
 
    int i1=0,i2=maxP-1;

    if(meshType=="original") {
      maxP = mc.nSurfacesInitial(); 
      i1 = mc.getIdx1st_s();
      i2 = mc.getIdxLst_s();
    }

    for(int i=i1; i<=i2; i++) {
      double s,x1;
      double xmin=0.05, xmax=0.975;
      if(meshType=="s") { 
        x1 = xmin+(xmax-xmin)*i/(maxP-1);  // x=s
        s = x1;
        x.push_back(s);
      } else if(meshType=="sqrt(s)") {
        x1 = xmin+(xmax-xmin)*i/(maxP-1);  // x=sqrt(s)
        s = x1*x1;
        x.push_back(x1);
      } else if(meshType=="original") {
        s=mc.get_s(i);
        x.push_back(s);
      }
      f_trap.push_back(mc.ftrapped(s));
      MConf::Vector3d e;
      e.set(mc.Gw(s),mc.Gs(s),mc.Gv(s));
      Gamma_factor.push_back(e);     
      switch(rdef) {
        case 1:
          epsName = "eps(dkes)  ";
          eps_eff.push_back(mc.epsEffDkes(s));
          break;
        case 2:               
          epsName = "eps(vmec)  ";
          eps_eff.push_back(mc.epsEffVmec(s));
          break;
        case 3:
          epsName = "eps(graz)  ";
          eps_eff.push_back(mc.epsEffGraz(s));
          break;
        case 4:
          epsName = "eps(dkes)  eps(vmec)  eps(graz)  ";
          e.set(mc.epsEffDkes(s),mc.epsEffVmec(s),mc.epsEffGraz(s));
          eps_eff_all.push_back(e);
          break;
        default:
          epsName = "eps(dkes)  ";
          eps_eff.push_back(mc.epsEffDkes(s));
          break;
       }
    }
    if(option3=="p") {
      epsName += "f_trap     ";
    }
    if(option2=="t") {
      epsName+= "Gw         Gs         Gv         ";
    }
    std::cerr<<" done, time/s="<<double(clock() - tStart)/CLOCKS_PER_SEC<<std::endl;
    return true;
  } 
};

static bool createList(const std::string &fname, std::vector<std::string> &list);
static void createFileList(std::string dir, std::string mask, std::vector<std::string> &list);
static std::string getPath(std::string fileName);
static const char *fl_filename_name(const char *name);
template <class T> void swap (T &v1, T &v2) {T w = v1;  v1 = v2;  v2 = w;};


//*********************************************************************
//  effective helical ripple
//
int main(int argc, char *argv[])
{
  const char * Version = MConf::CStconfig::getVersion();

  std::string option = (argc>=2)?argv[1]:" ";

  if(argc<2||option=="--ver"||option=="-ver") {
    std::cout<<"Version "<<Version<<"  year.month"<<std::endl;
    return 1;
  }

  if(argc<2||option=="--help"||option=="-help"||option=="help"||option=="?"||option=="-?") {
    std::cout<<"Version "<<Version<<std::endl;
    std::cerr<<"Calculate effective helical ripple and gradB drift velocity. 2010, Turkin."<<std::endl;
    
    std::cerr<<"The program is based on papers:"<<std::endl;
    std::cerr<<"V.V.Nemov,S.V.Kasilov,W.Kernbichler,M.F.Heyn,"<<std::endl;
    std::cerr<<"Physics of Plasmas, vol.6,4622(1999);"<<std::endl<<std::endl;
    
    std::cerr<<"The gradB drift velocity of trapped particles in stellarators"<<std::endl;
    std::cerr<<"Physics of Plasmas, vol.12,112507(2005)"<<std::endl;
    std::cerr<<"http://link.aip.org/link/doi/10.1063/1.2131848"<<std::endl<<std::endl; 
    
    std::cerr<<"Usage:  eps [-m] [-p] [-t] [option1] [option2] filename"<<std::endl;
    std::cerr<<"   where filename is the directory or path, and the file name,"<<std::endl; 
    std::cerr<<"   which can include wildcard characters, for example,"<<std::endl; 
    std::cerr<<"   an asterisk or a question mark. The files must be in"<<std::endl; 
    std::cerr<<"   format of VMEC-wout-file, bc-file or new_boz-file."<<std::endl<<std::endl; 
    std::cerr<<"Usage2:  eps [option1] [-p] [-t] [-m] @filelist"<<std::endl;
    std::cerr<<"   where filelist is the name of a file that contains list of files"<<std::endl;
    std::cerr<<"   for which the calculations will be made."<<std::endl<<std::endl;;
    std::cerr<<"Options:"<<std::endl;
    std::cerr<<"  -m      produce minimal output if this option is supplied"<<std::endl;
    std::cerr<<"  -p      output fraction of trapped particles"<<std::endl;
    std::cerr<<"  -t      output average gradB drift velocity of trapped particles"<<std::endl;
    std::cerr<<"option1:"<<std::endl;
    std::cerr<<"  -d   default option,"<<std::endl;
    std::cerr<<"          use DKES definition of the minor radius r=sqrt(F(s)/(pi*B_00(s)))"<<std::endl;
    std::cerr<<"  -v   use VMEC definition of the minor radius r=a*sqrt(s)"<<std::endl;
    std::cerr<<"  -g   use Graz definition of the minor radius dr=dV(s)/A(s)=ds/<|grad(s)|>"<<std::endl;
    std::cerr<<"  -a   output eps_eff for all definitions of the minor radius"<<std::endl;
    std::cerr<<" where F is the toroidal flux,"<<std::endl; 
    std::cerr<<"       s is the toroidal flux normalized to the value"<<std::endl; 
    std::cerr<<"         at the last closed magnetic surface,"<<std::endl;
    std::cerr<<"       a is the plasma radius,"<<std::endl;
    std::cerr<<"       A is the area of the flux surface."<<std::endl<<std::endl;
    std::cerr<<"option2:"<<std::endl;
    std::cerr<<"  -Gtype  type of the the radial grid,"<<std::endl;
    std::cerr<<"          where type is one from: original, s, sqrt(s);"<<std::endl;
    std::cerr<<"            original instruct to output eps_eff using s-grid from the file,"<<std::endl;
    std::cerr<<"            s instruct to output eps_eff using equidistant s-grid,"<<std::endl;
    std::cerr<<"            sqrt(s) instruct to output eps_eff using equidistant sqrt(s)-grid ,"<<std::endl;
    std::cerr<<"  -nxx    number of radial grid points, where xx is the number;"<<std::endl;
    std::cerr<<"          this option is ignored if -Goriginal is used"<<std::endl;
    std::cerr<<"Note: Only the first letter of the option must be supplied."<<std::endl;
    std::cerr<<"   Option letters may be combined, e.g., '-vtm' is equivalent to '-v -t -m'."<<std::endl;
    std::cerr<<"   Option letters are case sensitive."<<std::endl;
    std::cerr<<"Examples:"<<std::endl;
    std::cerr<<"     eps  -ver"<<std::endl;
    std::cerr<<"     eps  -help"<<std::endl;
    std::cerr<<"     eps  w7x-sc1.bc"<<std::endl;
    std::cerr<<"     eps -m -Goriginal w7x-sc1.bc"<<std::endl;
    std::cerr<<"     eps -m -Gs -n20 w7x-sc1.bc"<<std::endl;
    std::cerr<<"     eps -a w7x-sc1.bc"<<std::endl;
    std::cerr<<"     eps  path"<<std::endl;
    std::cerr<<"     eps  \"path/wout*\""<<std::endl;
    std::cerr<<"     eps  \"path/*.bc\""<<std::endl;
    std::cerr<<"     eps  \"path/pattern\""<<std::endl;
    std::cerr<<"     eps  \"path/pattern\" > outputfile"<<std::endl;
    std::cerr<<"     eps  -v @w7xlist > w7xeps.txt"<<std::endl<<std::endl;
    return 1;
  }
 
  std::string temp=""; 
  std::string option1="d"; 
  std::string option2=" "; 
  std::string option3=" "; 
  std::string option_m=" "; 
  std::string option_G=" "; 
  std::string meshType="sqrt(s)"; 
  int c;
  int nPoints=0;
  while ((c = getopt(argc, argv, _T("dvgamptG:n:"))) != -1) {
    switch (c) {
      case _T('d'):
      case _T('v'):
      case _T('g'):
      case _T('a'):
        option1 = c;
        break;
      case _T('t'):
        option2 = c;
        break;
      case _T('p'):
        option3 = c;
        break;
      case _T('m'):
        option_m = c;
        break;
      case _T('G'):
        option_G = c;
        meshType = optarg;
        break;      
      case _T('n'):
        nPoints = atoi(optarg);
        break;      
      case _T('?'):
#ifdef WIN32        
       std::cerr<<_T("ERROR: invalid option ")<<argv[optind-1]<<std::endl;
#endif        
        break;
      default: 
        temp = c;
        std::cerr<<_T("WARNING: no handler for option ")<<temp<<std::endl;
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

//#define _DEBUG  1

#ifdef _DEBUG
  std::cerr<<"option1 : "<<option1<<std::endl;
  std::cerr<<"option2 : "<<option2<<std::endl;
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

  int maxPoints = nPoints==0?26:nPoints;
  if( !(meshType=="original"||meshType=="s"||meshType=="sqrt(s)") )meshType="sqrt(s)";
  if(filelist.size()>1 && meshType=="original") meshType="s";

  Device wrk(maxPoints,1e-6,option1,option2,option3,meshType);
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

  // print header
  if(option_m!="m") {
    std::cout<<"% "<<dev.size()<<"  - number of devices (columns)"<<std::endl;
    std::cout<<"% "<<dev[0].points()<<" - number of radial points (rows)"<<std::endl;
    if(meshType=="sqrt(s)") 
      std::cout<<"%xaxis: r/a"<<std::endl;
    else
      std::cout<<"%xaxis: s"<<std::endl;
    std::cout<<"%yaxis: eps_eff"<<std::endl;
    std::cout<<"%names:   ";
    for(unsigned int k=0; k<dev.size(); k++) {   // print device names
      std::cout<<" "<<dev[k].name;
    }
    std::cout<<std::endl;
  }

  std::cout<<"% version "<<Version<<std::endl;
  if(meshType=="sqrt(s)") 
    std::cout<<"% sqrt(s)  ";
  else
    std::cout<<"% s        ";
  for(unsigned int k=0; k<dev.size(); k++)  {    // print eps_eff names
    std::cout<<dev[k].epsName;
  }
  std::cout<<std::endl;

  // print data
#ifdef _MSC_VER
  std::cout.precision(3);
#else
  std::cout.precision(4);
#endif
  std::cout.setf(std::ios::scientific);
  for(unsigned int i=0; i<dev[0].points(); i++) {
    std::cout<<dev[0].x[i];
    for(unsigned int k=0; k<dev.size(); k++) {
      if(option1=="a") {   
        std::cout<<" "<<dev[k].eps_eff_all[i];
      }
      else { 
        std::cout<<" "<<dev[k].eps_eff[i];
      }
      if(option3=="p") {
          std::cout<<" "<<dev[k].f_trap[i];
      }
      if(option2=="t") {
        std::cout<<" "<<dev[k].Gamma_factor[i];
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
