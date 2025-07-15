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

         
  fname='w7x_ref_9.bc';
  %fname='wout_w7x.txt';
  % load the magnetic configuration file
  % @return -- if the function succeeds, the return value is 
  % the address of C3dMesh object;  zero otherwise.
  MC = calllib(mconf,'MCload',fname); 
  % test the MC before next calls, it must be non-zero
  if MC == 0
    error('mconf_matlab: Could not load magnetic configuration')
  end

  % set accuracies and create the mesh
  epsTrunc = 1e-6;
  epsA = 1e-4; 
  calllib(mconf,'MCtruncate',MC,epsTrunc);  % truncate spectrum
  calllib(mconf,'MCsetAccuracy',MC,epsA);  % set accuracy of coordinate transformation im meters
  dr=0.02;
  dz=0.02; % 2cm
  dfi=3.14/180; % 1degree
  tic;
  disp(['3d-mesh is creating.....' ]);
  calllib(mconf,'MCcreateMeshUsingSymmetry',MC,dr,dz,dfi);  
  toc;

%***********************************************************
%***********************************************************

  tic;
  Npol = 181;
  Ntor = 721;
  x = nan(Ntor,Npol); % phi,theta
  y = nan(Ntor,Npol);
  z = nan(Ntor,Npol);
  Test = nan(Ntor,Npol);

  mix=[0,0,0];
  r=[0,0,0];
  dGdx=[0,0,0];
  dGdy=[0,0,0];
  dGdz=[0,0,0];
  gradS=[0,0,0];

  sLabel = 1;
  for j=1:Npol
    for i=1:Ntor
      phi = (i-1)*2*pi/(Ntor-1);
      theta = (j-1)*2*pi/(Npol-1);
      mix = [sLabel,theta,phi];
      [mix,r]=calllib(mconf,'MCmix2xyz',MC,mix,r); 
      
 %get grads and partial derivatives in cartesian coordinates
      [s,r,gradS,dGdx,dGdy,dGdz]=calllib(mconf,'M3DgetdGradsxyz',MC,r,gradS,dGdx,dGdy,dGdz);
      
      x(i,j) = r(1);
      y(i,j) = r(2);
      z(i,j) = r(3);
      Test1 = sqrt(dot(gradS,gradS));
      Test(i,j) = Test1;
    end
  end
  toc;
  figure;
  surf(x,y,z,Test);
  colormap parula
  axis equal;
  shading interp;  
%  camlight left;
  

  
  calllib(mconf,'MCfree',MC);
  unloadlibrary(mconf)

