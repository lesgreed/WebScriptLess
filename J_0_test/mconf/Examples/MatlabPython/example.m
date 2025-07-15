%
% Under Linux run the following command before starting Matlab
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

  % port AEO41
  r0=[-3.92678,-6.29125,1.38104];   % 1st point
  r1=[-3.14302,-4.97023,0.52986];   % 2nd point

  % trace plasma along the ray through the port AEL41  
  r0=[-2.39133,-2.32718,-0.12071];  % 1st point 
  r1=[-3.37847,-4.27681, 0.17038];  % 2nd point
  rd=r1-r0;                         % ray direction
  rd=rd/norm(rd);                   % normalize
  % find entry point of the ray into plasma
  entry=[0,0,0]; % ray entry
  exit =[0,0,0]; % ray exit
  [retcode,r0,rd,entry,exit]=calllib(mconf,'MCgetRayIntersectionPoints',MC,r0,rd,entry,exit);
  %  test retcode, it must be non-zero
  if retcode == 0
    unloadlibrary(mconf)
    error('mconf_matlab: ray does not hit plasma')
  end

  w = entry-r0;
  l0 = norm(w); 
  r = entry;    % r is the cartesian coordinates of a point lying on the ray.
  dl = 0.0001;    % 0.1mm step along the ray (only for test)
  dr = rd*dl; 
  disp(['  x          y            z        l         x        |B|' ])
  B = [0,0,0];                    
  gradB = [0,0,0];                    
  count=0;
  tic;
  for i=1:10000 
    % get  s, dV/ds, V, and vectors  B, gradB
    [s,r,    B]=calllib(mconf,'MCgetBxyz',    MC,r,B); 
    [s,r,gradB]=calllib(mconf,'MCgetGradBxyz',MC,r,gradB); 
  % V          =calllib(mconf,'MCVolume',     MC,s);  % V  is the  volume inside the surface s 
    Vp         =calllib(mconf,'MCVprime',     MC,s);  % Vp is the  dV/ds 
    x = sqrt(s); % x is the normalized plasma radius x=reff/a             
                 % s is the normalized toroidal flux 0<=s<=1
    dVdx = 2*x*Vp;
    x1(i)=x;
    l1(i)=l0+(i-1)*dl;
    B1(i)=norm(B);
    gB(i)=norm(gradB);
    if x>1, break, end
  %  disp([' ' num2str(r) '   ' num2str(l0+(i-1)*dl) '   ' num2str(x) '   ' num2str(norm(B))])    
    % get temperature and density
    %    T = temperature(x)    
    %    N = density(x)
    % do something using  B, gradB, s, dV/ds, V     
    r = r + dr; % advance along the ray
    count = count+1;
  end 
  disp(['tracing done, # of iteration '  num2str(count)  ]);
  toc;

  plot(l1,x1);
  pause;
  plot(l1,B1);
  pause;
  plot(l1,gB);
  
  calllib(mconf,'MCfree',MC);
  unloadlibrary(mconf)

