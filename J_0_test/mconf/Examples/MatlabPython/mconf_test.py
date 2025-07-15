# This example needs Python 3 
# This is working example of how to use the mconf.dll under Python

import sys
import os
import platform
import math
from ctypes import *

from osa import Client

os_type = platform.system()
is_64bits = sys.maxsize > 2**32

# it doesn't seem to work without the abspath (even normpath doesn't)
fname_lib = os.path.abspath(os.path.join('mconf', "mconf_matlab"))

if is_64bits: fname_lib += "64"
  
if os_type=='Windows': fname_lib += ".dll"
elif os_type=='Linux': fname_lib += ".so"
elif os_type=='Darwin':  fname_lib += ".dylib"

print("Loading library file : " + fname_lib)
mconf = cdll.LoadLibrary(fname_lib) 
#mconf = CDLL(fname_lib) 

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
mconf.MCepsEff = c_double

fname_vac   = "w7x-sc1(reduced).bc"
print("loading vacuum equilibrium : " + fname_vac)
fname_ctype = c_char_p(fname_vac.encode("utf-8"))

# load the magnetic configuration file
# @return -- if the function succeeds, the return value is 
# the address of C++ object;  zero otherwise.
mc = mconf.MCload(fname_ctype) # mc is like self in python
assert mc != 0, 'mconf: Could not load magnetic configuration'

unscaled_b0 = mconf.MCgetB0(mc,c_double(0.))
print("vacuum Bax="+str(unscaled_b0));

# set magnetic field 2.52T at toroidal angle 0
mconf.MCsetB0 (mc, c_double(2.52), c_double(0.))
#mconf.MCsetB00 (mc, c_double(2.52)) # does not seem to do the same thing!
scaled_b0 = mconf.MCgetB0(mc, c_double(0.))
print("scaled vacuum Bax="+str(scaled_b0)) 
assert scaled_b0 == 2.52, "Something is wrong with the B-Field scaling!" 

vac_Ip = mconf.MCIp(mc)
print("yields vac_Ip="+str(vac_Ip)) 

mconf.MCfree(mc)


print("\n#####################\n     beta file   \n#####################\n")

fname_beta = "w7x-sc1beta=0.02.bc"
print("loading beta case : " + fname_beta)
# load configuration with beta
fname_ctype=c_char_p(fname_beta.encode("utf-8"))

mc = mconf.MCload(fname_ctype) 
assert mc != 0, 'mconf: Could not load magnetic configuration'

unscaled_b0 = mconf.MCgetB0(mc, c_double(0.))
print("old Bax="+str(unscaled_b0)); 
print("old Ip = {:f}".format(mconf.MCIp(mc)))
mconf.MCsetIpLCMS(mc, c_double(vac_Ip))
scaled_b0 = mconf.MCgetB0(mc, c_double(0.))
print("new Bax="+str(scaled_b0))
print("new Ip = {:f}".format(mconf.MCIp(mc)))
assert mconf.MCIp(mc) == vac_Ip, "The poloidal current has not been updated correctly! We expected vac_Ip="+str(vac_Ip)

r = mconf.MCreff(mc, c_double(1.))
print("a="+str(r))

print("\n######################\n scale s & write file \n######################\n")

fname_out = "w7x-sc1beta=0.02,s=1.69.bc"
print("saving to : " + fname_beta)

fname_ctype=c_char_p(fname_out.encode("utf-8"))
s_new = 1.69
# this also writes a file
mconf.MCsetLCMS(mc, c_double(s_new), fname_ctype)
mconf.MCfree(mc);

print("Checking stored file")
mc = mconf.MCload(fname_ctype) 
assert mc != 0, 'mconf: Could not load magnetic configuration'


s=1/s_new
r = mconf.MCreff(mc,c_double(s))
print("r="+str(r)+" at s_new="+str(s))
print("vac_Ip = {:f}".format(mconf.MCIp(mc)))
assert mconf.MCIp(mc) == vac_Ip, "The updated file should have the vacuum field's Ip of " + \
                                 "{:f}.But it contains {:f}".format(vac_Ip, mconf.MCIp(mc))




###########################################################
################## NETCDF test
# load VMEC
vmec = Client("http://esb.ipp-hgw.mpg.de:8280/services/vmec_v8?wsdl")
vac_eq = "w7x_ref_27" 

# get the NETCDF file 
print("Getting netcdf file for vacuum equilibrium : " + vac_eq, end=" ... ")
wout_netcdf = vmec.service.getVmecOutputNetcdf(vac_eq)
print("Done!")

if not(os.path.exists("netcdf")) : os.makedirs("netcdf")
fname_vac_eq = os.path.join('netcdf', vac_eq + ".nc")
print(" Writing it to file " + fname_vac_eq, end="...")
with open(fname_vac_eq, "wb") as f : f.write(wout_netcdf)
print("Done!")

# test whether we can load the file 
fname_ctype=c_char_p(fname_vac_eq.encode("utf-8"))
mc = mconf.MCload(fname_ctype) 
assert mc != 0, 'mconf: Could not load magnetic configuration'


print("Passed!")

eps_eff = mconf.MCepsEff(mc, c_double(0.25))
print("eps_eff="+str(eps_eff))



 
  
