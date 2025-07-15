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
  epsTrunc = 1e-9;
  epsA = 1e-6; 
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

  phi = 0;
  sLabel = 1;
% find  min & max of R and Z
  mix=[0,0,0];
  r=[0,0,0];
  for i=1:360
    th = i*6.28318531/360;  
    mix = [1,th,phi];
    [mix,r]=calllib(mconf,'MCmix2xyz',MC,mix,r); 
    R1(i) = sqrt(r(1)^2+r(2)^2);
	Z1(i) = r(3);
  end

  Rmin = min(R1);
  Rmax = max(R1);
  Zmin = min(Z1);
  Zmax = max(Z1);
  

  tic;
  NR = floor((Rmax-Rmin)/0.01)+	1
  NZ = floor((Zmax-Zmin)/0.01)+ 1

  %R = nan(NR); 
  %Z = nan(NZ);
  jper = nan(NZ,NR);
  jpar = nan(NZ,NR);
  js   = nan(NZ,NR);

  B=[0,0,0];
  dBdx=[0,0,0];
  dBdy=[0,0,0];
  dBdz=[0,0,0];
  gradS=[0,0,0];
  curlB = [0,0,0];
  b = [0,0,0];
  n = [0,0,0];
  
  for j=1:NR
      R(j) = Rmin+(j-1)*0.01;
  end
  
  for i=1:NZ
    Z(i) = Zmin+(i-1)*0.01;
  end

  for j=1:NR
    for i=1:NZ
      r = [R(j)*cos(phi),R(j)*sin(phi),Z(i)];
      [s,r,B,dBdx,dBdy,dBdz,gradS]=calllib(mconf,'MCgetdB_Gradsxyz',MC,r,B,dBdx,dBdy,dBdz,gradS);
      
	  curlx = dBdy(3)-dBdz(2);  % curl(B)_x = dBz/dy - dBy/dz
      curly = dBdz(1)-dBdx(3);  % curl(B)_y = dBx/dz - dBz/dx
      curlz = dBdx(2)-dBdy(1);  % curl(B)_z
      curlB = [curlx,curly,curlz];

	  b = B/sqrt(dot(B,B));
      n = gradS/sqrt(dot(gradS,gradS));

      j_perToS = dot(curlB,n);
      j_par    = dot(curlB,b);
      j_per    = curlB - j_perToS*n -  j_par*b;
      js(i,j)    =  j_perToS;  % must be zero
      jpar(i,j)  =  j_par;
      jper(i,j)  =  sqrt(dot(j_per,j_per));
    end
  end
  toc;
  
  size(R)
  size(Z)
  size(jpar)
  
  figure;
  contourf(R,Z,jpar);
  
  figure;
  contourf(R,Z,jper);

% % x = linspace(-2*pi,2*pi);
% % y = linspace(0,4*pi);
% % [X,Y] = meshgrid(x,y);
% % A = sin(X)+cos(Y);

% % figure
% % contour(X,Y,A)
  
  
  calllib(mconf,'MCfree',MC);
  unloadlibrary(mconf)

