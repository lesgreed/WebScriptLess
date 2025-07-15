%
% 
% setenv LD_PRELOAD  /usr/lib/libstdc++.so.6
%
% This is working example of how to use the mconf_matlab.dll
  clear all
  if strcmp(computer,'GLNX86')
    loadlibrary('mconf_matlab.so','mconf_matlab.h');
    mconf = 'mconf_matlab';
  end
  if strcmp(computer,'PCWIN')
    loadlibrary('mconf_matlab.dll','mconf_matlab.h');
    mconf = 'mconf_matlab';
  end  
  if strcmp(computer,'PCWIN64')
    loadlibrary('mconf_matlab64.dll','mconf_matlab64.h');
    mconf = 'mconf_matlab64';
  end  
  if strcmp(computer,'GLNXA64')
    loadlibrary('mconf_matlab64.so','mconf_matlab64.h');
    mconf = 'mconf_matlab64';
  end  
  
  if ~libisloaded(mconf)
    error('Could not find mconf_matlab')
  end
  libfunctions(mconf,'-full');

         
  fname='w7x-sc1beta=0.02.bc';
  
  
  %fname='wout_w7x.txt';
  % load the magnetic configuration file
  % @return -- if the function succeeds, the return value is 
  % the address of C3dMesh object;  zero otherwise.
  MC = calllib(mconf,'MCload',fname); 
  % test the MC before next calls, it must be non-zero
  if MC == 0
    error('mconf_matlab: Could not load magnetic configuration')
  end
 
 
 
  pi = 3.1415926535897932384626433832795
  
  degree = pi/180
  
  a = 1*degree 
  
  mag = [0.5,0,0];  % s, theta, cyl_angle
  r0 = [0,0,0];
  r  = r0;
  l = 0;
  
  [mag,r0]=calllib(mconf,'MCmix2xyz',MC,mag,r0); 
  
   N= 300;
   dtheta=20*degree/N;
   
  for i=1:N
    mag(2) =  mag(2) + dtheta;
    [mag,r]=calllib(mconf,'MCmix2xyz',MC,mag,r); 
    w = r-r0;
    dl = norm(w);
    l = l + dl;
    r0 = r;
  end

  l
  
  calllib(mconf,'MCfree',MC);
  unloadlibrary(mconf)

