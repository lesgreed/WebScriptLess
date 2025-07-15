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

# load the magnetic configuration file
# @return -- if the function succeeds, the return value is 
# the address of C++ object;  zero otherwise.
mc = mconf.MCload(fname_ctype) 
assert mc != 0, 'mconf: Could not load magnetic configuration'


print("Passed!")

eps_eff = mconf.MCepsEff(mc, c_double(0.25))
print("eps_eff="+str(eps_eff))

mconf.MCfree(mc)
