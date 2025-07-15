# This example needs Python 3 
# This is working example of how to use the mconf.dll under Python

import sys
import platform
import math
import numpy as np
from ctypes import *

os = platform.system()
is_64bits = sys.maxsize > 2**32

if is_64bits:
  mc_handle = c_longlong
  if os=='Windows':
    mconf = cdll.LoadLibrary("mconf_matlab64.dll")
  elif os=='Linux':
    mconf = cdll.LoadLibrary("mconf_matlab64.so")  
  elif os=='Darwin':
    mconf = cdll.LoadLibrary("mconf_matlab64.dylib",".")
else:
  mc_handle = c_long
  if os=='Windows':
    mconf = cdll.LoadLibrary("mconf_matlab.dll")
  elif os=='Linux':
    mconf = cdll.LoadLibrary("mconf_matlab.so")  

# https://docs.python.org/3/library/ctypes.html#module-ctypes
    
mconf.MCgetBxyz.restype = c_double
mconf.MCgetB00.restype = c_double
mconf.MCVprime.restype = c_double
mconf.MCVolume.restype = c_double
mconf.MCtorFlux2polFlux.restype = c_double
mconf.M3Dxyz2s.restype = c_double
mconf.MCxyz2s.restype = c_double
mconf.M3Dxyz2s.argtypes = [mc_handle,c_double*3]

if is_64bits:
  mconf.MCload.restype = c_longlong   # for 64-bit architecture
else:
  mconf.MCload.restype = c_long

fname=c_char_p(b"w7x-sc1beta=0.02.bc")

# load the magnetic configuration file
# @return -- if the function succeeds, the return value is 
# the address of C++ object;  zero otherwise.
#mc = mconf.MCload(fname) # mc is like self in python

mc = mconf.MCload(c_char_p(b"w7x_ref_1"))
if mc == 0:
  print ('mconf: Could not load magnetic configuration')

B00 = mconf.MCgetB00(mc);
print(B00) 

# set magnetic field 2.5T at toroidal angle 0
##mconf.MCsetB0 (mc,c_double(2.5), c_double(0.));
B00 = mconf.MCgetB00(mc);
print(B00) 

mconf.MCtruncate(mc,c_double(1e-7))
mconf.MCsetAccuracy(mc,c_double(1e-4))
mconf.MCcreateMeshUsingSymmetry(mc,c_double(0.02),c_double(0.02),c_double(np.pi/180.))

m = (c_double*3)()
m[0]=6.2
m[1]=0.
m[2]=0.

s = c_double
s = mconf.M3Dxyz2s(mc,m);
print(s) 
s = mconf.MCxyz2s(mc,m);
print(s) 