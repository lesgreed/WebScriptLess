
%  Linux users: run the following command 
% 
% setenv LD_PRELOAD  /usr/lib/libstdc++.so.6
%
% before starting matlab

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
       
fname='w7x-sc1(reduced).bc';
% load the magnetic configuration file
% @return -- if the function succeeds, the return value is 
% the address of C3dMesh object;  zero otherwise.
MC = calllib(mconf,'MCload',fname); 

% please test the MC before next calls, it must be non-zero

B00 = calllib(mconf,'MCgetB00',MC);

% set accuracies and create the mesh
epsTrunc = 1e-7;
epsA = 1e-4; 
calllib(mconf,'MCtruncate',MC,epsTrunc);  % truncate spectrum
calllib(mconf,'MCsetAccuracy',MC,epsA);  % set accuracy of coordinate transformation im meters
dr=0.02;
dz=0.02; % 2cm
dfi=3.14/180; % 1degree
tic;
calllib(mconf,'MCcreateMeshUsingSymmetry',MC,dr,dz,dfi);  
toc;


% get B-field
xyz=[6,0,0];
B = [0,0,0];                    
pB = libpointer('doublePtr',B);
s=calllib(mconf,'M3DgetBxyz',MC,xyz,pB);
s; % print s -- norm. tor.flux
B = pB.Value; % one have to assign value, pB is a real pointer!!
B;

% next is also possible call
[s,xyz,B]=calllib(mconf,'M3DgetBxyz',MC,xyz,B); 
B;

% find entry point of the ray into plasma
r0=[7,0,0];  % initial point
rd=[-1,0,0]; % ray direction
entry=[0,0,0]; % ray entry
pentry = libpointer('doublePtr',entry);
calllib(mconf,'M3DgetRayEntryPoint',MC,r0,rd,pentry);
pentry.Value;          % print result
entry = pentry.Value;  % assign result
% next is also possible call
[retcode,r0,rd,entry]=calllib(mconf,'M3DgetRayEntryPoint',MC,r0,rd,entry);

% get s,B,grads  in point Entry
grads=1:3
[s,entry,B,grads]=calllib(mconf,'M3DgetBGradsxyz',MC,entry,B,grads); grads

% Advance entry
entry = entry + 0.1*rd; % 10cm inside
[s,entry,B,grads]=calllib(mconf,'M3DgetBGradsxyz',MC,entry,B,grads); grads

% get r_effective
r_eff=calllib(mconf,'MCreff',MC,s);

unloadlibrary(mconf)

