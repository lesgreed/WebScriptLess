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

mconf.MCgetBxyz.restype = c_double
mconf.MCgetB00.restype = c_double
mconf.MCVprime.restype = c_double
mconf.MCVolume.restype = c_double
mconf.MCtorFlux2polFlux.restype = c_double
mconf.MCreff.restype = c_double
mconf.MCxyz2s.restype = c_double

if is_64bits:
  mconf.MCload.restype = c_longlong   # for 64-bit architecture
else:
  mconf.MCload.restype = c_long

fname=c_char_p(b"w7x-sc1beta=0.02.bc")

# load the magnetic configuration file
# @return -- if the function succeeds, the return value is 
# the address of C++ object;  zero otherwise.
mc = mconf.MCload(fname) # mc is like self in python
if mc == 0:
  print ('mconf: Could not load magnetic configuration')

B00 = mconf.MCgetB00(mc);
print(B00) 

# set magnetic field 2.5T at toroidal angle 0
mconf.MCsetB0 (mc,c_double(2.5), c_double(0.));
B00 = mconf.MCgetB00(mc);
print(B00) 

# trace plasma along the ray through the port AEL41
r0 = (c_double*3)(-2.39133,-2.32718,-0.12071)  # 1st point 
r1 = (c_double*3)(-3.37847,-4.27681, 0.17038)  # 2nd point
rd = (c_double*3)()

#rd = r1 - r0                          # ray direction
rd[0] = r1[0] - r0[0]
rd[1] = r1[1] - r0[1]
rd[2] = r1[2] - r0[2]
#rd = rd/rd.abs()     # normalize
rdabs = math.sqrt(rd[0]**2+rd[1]**2+rd[2]**2)
rd[0] = rd[0]/rdabs  # normalize
rd[1] = rd[1]/rdabs  # normalize
rd[2] = rd[2]/rdabs  # normalize

print (rd[0], rd[1], rd[2])

# find entry point of the ray into plasma
entry=(c_double*3)() # ray entry
exit =(c_double*3)() # ray exit
retcode = mconf.MCgetRayIntersectionPoints(mc,r0,rd,entry,exit)
#  test retcode, it must be non-zero
if retcode == 0:
  print ('mconf: ray does not hit plasma')

print( entry[0], entry[1], entry[2])

# ****************************************************

maxPnt = 2000
dl = 0.01;     #  1cm
r0 = entry
          
# trace plasma along the ray
r = (c_double*3)()
B = (c_double*3)()
for i in range(0,maxPnt):    
  lng  = i*dl          # length from r0 to r
  r[0] = r0[0] + lng*rd[0]      # move along the ray
  r[1] = r0[1] + lng*rd[1]      # move along the ray
  r[2] = r0[2] + lng*rd[2]      # move along the ray
  
######################################  
  ##s = mconf.MCxyz2s(mc,r) 
  s = mconf.MCgetBxyz(mc,r,B)
  r = mconf.MCreff(mc,s) 
 ## Er = ErFunc(r)
######################################  

  if s>1:  break     # break if not inside plasma
  V = mconf.MCVolume(mc,c_double(s))  # V  is the  volume inside the surface s 
  Vp = mconf.MCVprime(mc,c_double(s)) # Vp is the  dV/ds 
  #x = sqrt(s)        # x is the normalized plasma radius x=reff/a  
  #n = ne(x)         # density
  #t = Te(x)         # temperature
  sPol = mconf.MCtorFlux2polFlux(mc,c_double(s)) # sPol is the normalized poloidal flux, where s is the normalized toroidal flux.
  print( lng, s,r,sPol, B[0], B[1], B[2]  )


# ****************************************************
# plot flux surface s=0.5 at cyl. angle 2degree 
phi = 2*6.28318531/360    # 2 degree 
s   = 0.5         

r = (c_double*3)()
m = (c_double*3)()
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
  print( R, Z ) # cyl. coordinates R,Z
 
  