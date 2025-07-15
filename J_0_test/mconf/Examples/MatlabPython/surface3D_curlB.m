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

  tic;
  Npol = 181;
  Ntor = 721;
  x = nan(Ntor,Npol); % phi,theta
  y = nan(Ntor,Npol);
  z = nan(Ntor,Npol);
  jper = nan(Ntor,Npol);
  jpar = nan(Ntor,Npol);
  js   = nan(Ntor,Npol);

  mix=[0,0,0];
  r=[0,0,0];
  B=[0,0,0];
  dBdx=[0,0,0];
  dBdy=[0,0,0];
  dBdz=[0,0,0];
  gradS=[0,0,0];
  curlB = [0,0,0];
  b = [0,0,0];
  n = [0,0,0];

  sLabel = 0.25;
  for j=1:Npol
    for i=1:Ntor
      phi = (i-1)*2*pi/(Ntor-1);
      theta = (j-1)*2*pi/(Npol-1);
      mix = [sLabel,theta,phi];
      [mix,r]=calllib(mconf,'MCmix2xyz',MC,mix,r); 
      [s,r,B,dBdx,dBdy,dBdz,gradS]=calllib(mconf,'MCgetdB_Gradsxyz',MC,r,B,dBdx,dBdy,dBdz,gradS);
      
      x(i,j) = r(1);
      y(i,j) = r(2);
      z(i,j) = r(3);
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
  figure;
  surf(x,y,z,jpar);
  colormap parula
  axis equal;
  shading interp;  
%  camlight left;
  
  figure;
  surf(x,y,z,jper);
  colormap parula
  axis equal;
  shading interp;  
%  camlight left;

  
  calllib(mconf,'MCfree',MC);
  unloadlibrary(mconf)

