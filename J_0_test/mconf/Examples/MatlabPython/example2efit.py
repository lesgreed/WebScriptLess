
# This is working example of how to use the mconf.dll

#  Numpy support for interfacing with ctypes.

import sys
import platform
import math
import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

os = platform.system()
is_64bits = sys.maxsize > 2**32

# load the library, using numpy mechanisms
if is_64bits:   # 64-bit architecture
  mc_handle = c_longlong
  if os=='Windows':
    mconf = npct.load_library("mconf_matlab64.dll",".")
  elif os=='Linux':
    mconf = npct.load_library("mconf_matlab64.so",".")  
  elif os=='Darwin':
    mconf = npct.load_library("mconf_matlab64.dylib",".")
else:
  mc_handle = c_long
  if os=='Windows':
    mconf = npct.load_library("mconf_matlab.dll",".")
  elif os=='Linux':
    mconf = npct.load_library("mconf_matlab.so",".")  

vec3 = npct.ndpointer(dtype=np.float64, ndim=1, flags='CONTIGUOUS')

# setup the return typs and argument types
mconf.MCload.restype = mc_handle  
mconf.MCload.argtypes = [c_char_p] 

mconf.MCgetRayIntersectionPoints.restype = c_int
mconf.MCgetRayIntersectionPoints.argtypes = [mc_handle,vec3,vec3,vec3,vec3]

mconf.MCgetB00.restype = c_double
mconf.MCgetB00.argtypes = [mc_handle]

mconf.MCsetB0.restype = c_double
mconf.MCsetB0.argtypes = [mc_handle,c_double,c_double]

mconf.MCgetBxyz.restype = c_double
mconf.MCgetBxyz.argtypes = [mc_handle,vec3,vec3]

mconf.MCVprime.restype = c_double
mconf.MCVprime.argtypes = [mc_handle,c_double]

mconf.MCVolume.restype = c_double
mconf.MCVolume.argtypes = [mc_handle,c_double]

mconf.MCtorFlux2polFlux.restype = c_double
mconf.MCtorFlux2polFlux.argtypes = [mc_handle,c_double]

# EQDSK file from EFIT
fname=c_char_p("g033068_02750.txt")
# load the magnetic configuration file
# @return -- if the function succeeds, the return value is 
# the address of C++ object;  zero otherwise.
mc = mconf.MCload(fname) # mc is like self in python
if mc == 0:
  print 'mconf: Could not load magnetic configuration'

B00 = mconf.MCgetB00(mc);
print B00 

# set magnetic field 2.5T at toroidal angle 0
##mconf.MCsetB0 (mc,c_double(2.5), c_double(0.));
B00 = mconf.MCgetB00(mc);
print B00 

# trace plasma along the ray through the port AEL41(W7X only)
# using Cartesian coordinates
r0 = np.array([-2.39133,-2.32718,-0.12071],dtype=np.float64)
r1 = np.array([-3.37847,-4.27681, 0.17038],dtype=np.float64)
rd = r1 - r0
rd = rd/np.sqrt(np.dot(rd,rd)) # normalize
#  !! Reverse direction for testing EFIT file
print "direction",rd
print "!! Reverse direction for testing EFIT file ", fname.value
rd = - rd
print "direction",rd

# find entry point of the ray into plasma
entry= np.empty_like(r0) # ray entry
exit = np.empty_like(r0) # ray exit
retcode = mconf.MCgetRayIntersectionPoints(mc,r0,rd,entry,exit)
#  test retcode, it must be non-zero
if retcode == 0:
  print 'mconf: ray does not hit plasma'

print entry, exit, retcode

# ########################################################

maxPnt = 2000
dl = 0.01;     #  1cm
r0 = entry
          
# trace plasma along the ray
r = np.empty_like(r0)
rB = np.empty_like(r0)
for i in xrange(0,maxPnt):    
  lng  = i*dl          # length from r0 to r
  r = r0 + lng*rd      # move along the ray
  s = mconf.MCgetBxyz(mc,r,rB)
  if s>1:  break     # break if not inside plasma
  V =  mconf.MCVolume(mc,c_double(s))  # V  is the  volume inside the surface s 
  Vp = mconf.MCVprime(mc,c_double(s)) # Vp is the  dV/ds 
  #x = sqrt(s)       # x is the normalized plasma radius x=reff/a  
  #n = ne(x)         # density
  #t = Te(x)         # temperature
  sPol = mconf.MCtorFlux2polFlux(mc,c_double(s)) # sPol is the normalized poloidal flux, where s is the normalized toroidal flux.
  print lng, s, sPol, rB 

  
mconf.MCmix2xyz.argtypes = [mc_handle,vec3,vec3]
 
# ****************************************************
# plot flux surface s=0.5 at cyl. angle 2degree 
phi = 2*6.28318531/360    # 2 degree 
s   = 0.5         

m = np.empty_like(r)
maxPnt = 200
dth =6.28318531/(maxPnt-1)
for i in range(0,maxPnt):    
  th = i*dth  
  m[0] = s
  m[1] = th
  m[2] = phi
  
  mconf.MCmix2xyz(mc,m,r)
  R = math.sqrt(r[0]**2+r[1]**2)
  Z = r[2]
  print R, Z  # cyl. coordinates R,Z
