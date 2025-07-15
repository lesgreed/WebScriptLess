# This example needs Python 3 
# This is working example of how to use the mconf.dll under Python

import sys
import platform
import math
from ctypes import *

os = platform.system()
is_64bits = sys.maxsize > 2**32

if is_64bits:
  if os=='Windows':
    mconf = cdll.LoadLibrary("mconf_matlab64.dll")
  elif os=='Linux':
    mconf = cdll.LoadLibrary("mconf_matlab64.so")  
  elif os=='Darwin':
    mconf = cdll.LoadLibrary("mconf_matlab64.dylib",".")
else:
  if os=='Windows':
    mconf = cdll.LoadLibrary("mconf_matlab.dll")
  elif os=='Linux':
    mconf = cdll.LoadLibrary("mconf_matlab.so")  

if is_64bits:
  mconf.MCload.restype = c_longlong   # for 64-bit architecture
else:
  mconf.MCload.restype = c_long
  
mconf.MCgetBxyz.restype = c_double
mconf.MCgetB00.restype = c_double
mconf.MCgetB0.restype  = c_double
mconf.MCVprime.restype = c_double
mconf.MCVolume.restype = c_double
mconf.MCtorFlux2polFlux.restype = c_double
mconf.MCIp.restype = c_double
mconf.MCreff.restype = c_double


fname=c_char_p(b"w7x-sc1(reduced).bc")

# load the magnetic configuration file
# @return -- if the function succeeds, the return value is 
# the address of C++ object;  zero otherwise.
mc = mconf.MCload(fname) # mc is like self in python
if mc == 0:
  print ('mconf: Could not load magnetic configuration')
print(fname); 

B0 = mconf.MCgetB0(mc,c_double(0.));
print("Bax="+str(B0)); 

# set magnetic field 2.52T at toroidal angle 0
mconf.MCsetB0 (mc,c_double(2.52), c_double(0.));
B0 = mconf.MCgetB0(mc,c_double(0.));
print("Bax="+str(B0)); 

Ip = mconf.MCIp(mc);
print("Ip="+str(Ip)); 

mconf.MCfree(mc);

# load configuration with beta
fname=c_char_p(b"w7x-sc1beta=0.02.bc")

mc = mconf.MCload(fname) 
if mc == 0:
  print ('mconf: Could not load magnetic configuration')
print(fname); 

mconf.MCsetIpLCMS(mc,c_double(Ip))
B0 = mconf.MCgetB0(mc,c_double(0.));
print("Bax="+str(B0)); 

r = mconf.MCreff(mc,c_double(1.))
print("a="+str(r))

fname=c_char_p(b"w7x-sc1beta=0.02,s=1.69.bc")
s_new = 1.69
mconf.MCsetLCMS(mc,c_double(s_new))
##mconf.MCfree(mc);

##mc = mconf.MCload(fname) 
if mc == 0:
  print ('mconf: Could not load magnetic configuration')
print(fname); 

s=1/s_new
r = mconf.MCreff(mc,c_double(s))
print("r="+str(r)+" at s_new="+str(s))

fname=c_char_p(b"w7x-sc1beta=0.02,s=1.69__x.bc")
mconf.MCwrite(mc,fname)


 
  