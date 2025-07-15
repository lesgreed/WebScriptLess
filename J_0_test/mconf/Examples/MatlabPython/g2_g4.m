

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
 
  fname = 'boozer.jose_tok1.data'; 
%  fname='w7x-sc1(reduced).bc';
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
   
  calllib(mconf,'MCFbsSetIotaParam',      MC, 0 , 1, 150);
  calllib(mconf,'MCFbsSetMagnMomentParam',MC, 257, 15.0);
  calllib(mconf,'MCFbsSetSlabelParam',    MC, 0.0025, 0.9, 101);
  calllib(mconf,'MCFbsSetTracingParam',   MC, 100, 1, 0.01); % 0.01745->1degree  
  for i=1:N
    s1 = i*ds;
    s(i) = s1;
    g2(i)=calllib(mconf,'MCFbsg2',MC,s1); 
    g4(i)=calllib(mconf,'MCFbsg4',MC,s1,0.1); 
    ft(i)=calllib(mconf,'MCftrapped',MC,s1); 
    fbs(i)=calllib(mconf,'MCFbs',MC,s1); 
    B2(i)=calllib(mconf,'MCFbsB2',MC,s1); 
  end

  plot(s,B2);
  plot(s,g2);
  pause;
  plot(s,g4);
  % pause;
  % plot(s,ft);
  % pause;
  % plot(s,fbs);

  unloadlibrary(mconf)
