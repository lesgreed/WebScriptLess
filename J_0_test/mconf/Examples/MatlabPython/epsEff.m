

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
 
  fname = 'w7x_ref_9.bc'; 
% load the magnetic configuration file
% @return -- if the function succeeds, the return value is 
% the address of mconf object;  zero otherwise.
  MC = calllib(mconf,'MCload',fname); 

% test the MC before next calls, it must be non-zero
  if MC == 0
    error('W7-X library: Could not load magnetic configuration')
    return
  end
 
  B00 = calllib(mconf,'MCgetB00',MC);

  epsTrunc = 1e-7;
  calllib(mconf,'MCtruncate',MC,epsTrunc);  % truncate spectrum


  N= 100;
  ds=.95/N;
   
  for i=1:N
    s = i*ds;
    x(i) = sqrt(s);
    eps(i)=calllib(mconf,'MCepsEff',MC,s); 
  end

  plot(x,eps);
  unloadlibrary(mconf)
