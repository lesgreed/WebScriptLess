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
fname_lib = os.path.abspath(os.path.join("mconf_matlab64_netcdf"))

#if is_64bits: fname_lib += "64"
  
if os_type=='Windows': fname_lib += ".dll"
elif os_type=='Linux': fname_lib += ".so"
elif os_type=='Darwin':  fname_lib += ".dylib"

print("Loading library file : " + fname_lib)

mconf = cdll.LoadLibrary(fname_lib) 
#mconf = CDLL(fname_lib) 

mconf.MCload.argtypes =[c_char_p]
if is_64bits:
  mconf.MCload.restype = c_longlong   # for 64-bit architecture
else:
  mconf.MCload.restype = c_long
  
mconf.MCepsEff.argtypes = [c_longlong, c_double]
mconf.MCepsEff.restype = c_double
mconf.MCfree.argtypes = [c_longlong]

###########################################################
################## NETCDF test
# load VMEC
vmec = Client("http://esb.ipp-hgw.mpg.de:8280/services/vmec_v8?wsdl")
w7x_eq = "w7x_ref_27" 

# get the NETCDF file 
print("Getting netcdf file for vacuum equilibrium : " + w7x_eq, end=" ... ")
wout_netcdf = vmec.service.getVmecOutputNetcdf(w7x_eq)
print("Done!")

if not(os.path.exists("netcdf")) : os.makedirs("netcdf")
fname_w7x_eq = os.path.join('netcdf', w7x_eq + ".nc")
print("Writing it to file " + fname_w7x_eq, end="...")
with open(fname_w7x_eq, "wb") as f : f.write(wout_netcdf)
print("Done!")

fname_ctype=c_char_p(fname_w7x_eq.encode("utf-8"))

# load the magnetic configuration file
# @return -- if the function succeeds, the return value is 
# the address of C++ object;  zero otherwise.
mc = mconf.MCload(fname_ctype) 
assert mc != 0, 'mconf: Could not load magnetic configuration'
print("The file in netCDF format is loaded!")

eps_eff = mconf.MCepsEff(mc, c_double(0.25))
print("eps_eff="+str(eps_eff))

mconf.MCfree(mc)
