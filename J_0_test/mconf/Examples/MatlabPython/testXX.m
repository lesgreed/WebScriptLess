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

  % set accuracies and create the mesh
  epsTrunc = 1e-6;
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
  
  tic;
  % create copy of MConf object for multithreaded computation.
  %   No data are copied! Only the references inside objects are copied. 
  % This operation is very fast.
  MC2=calllib(mconf,'MCcopy',MC);  % create copy of MConf object for Multithreaded Computation
  % Now we have two instances of object mconf.
  % They can be used in parallel in multithreaded application, 
  % but I don't know how to do this in Matlab script.   
  % Vp1=calllib(mconf,'MCVprime',MC, s1);
  % Vp2=calllib(mconf,'MCVprime',MC2,s2);
  disp(['object copied' ]);
  toc;

  

  r0=[7.0,  0.0,-0.01];     % 1st point 
  r1=[6.93,0.07, 0.005];     % 2nd point
  rd=r1-r0;                         % ray direction
  rd=rd/norm(rd);                   % normalize
  % find entry point of the ray into plasma
  entry=[0,0,0]; % ray entry
  [retcode,r0,rd,entry]=calllib(mconf,'M3DgetRayEntryPoint',MC,r0,rd,entry);
 
 %  test retcode, it must be non-zero
  if retcode == 0
    unloadlibrary(mconf)
    error('mconf_matlab: ray does not hit plasma')
  end

  w = entry-r0;
  l0 = norm(w); 
  r = entry;     % r is the cartesian coordinates of a point lying on the ray.
  dl = 0.001;    % 0.1mm step along the ray
  dr = rd*dl; 
  disp(['  x          y            z        l         x        |B|' ])
  
  B = [0,0,0];                    
  gradB = [0,0,0];                    
  gradS = [0,0,0];                    
  count=0;
  tic;
  for i=1:10000 
   % get  s, dV/ds, V, and vectors  B, gradB

    [s,r,    B]=calllib(mconf,'MCgetBxyz', MC,r,B); 
    [s3,r,   B]=calllib(mconf,'M3DgetBxyz',MC,r,B);
   % [s,r,B,gradB,gradS]=calllib(mconf,'M3DgetBandGradBxyz',MC,r,B,gradB,gradS);
    
    x1(i)=sqrt(s);
    x3(i)=sqrt(s3);
    dx(i)=sqrt(s)-sqrt(s3);
    l1(i)=l0+(i-1)*dl;
    if s>1, break, end
    r = r + dr; % advance along the ray
    count = count+1;
  end 
  disp(['tracing done, # of iteration '  num2str(count)  ]);
  toc;

  plot(l1,x1);
  pause;
  plot(l1,x3);
  pause;
  plot(l1,dx);
  
  calllib(mconf,'MCfree',MC2);
  calllib(mconf,'MCfree',MC);
  unloadlibrary(mconf)

