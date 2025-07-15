#!/usr/bin/env python
import platform
import os
import ctypes as ct
import numpy as np

class Mconf_equilibrium:
    """
    Class provides access to VMEC (and other) equilibriums via Juri Turkin's mconf c-based libary.
    
    Takes an optional dic of parameters:
    mconf_config = {'B0': 2.52, #B field at magentic axis at toroidal angle B0_angle
              'B0_angle': 0.0,
             'extraLCMS': 1.2,   #Flux surfaces extrapolation parameter (s_max)
              'accuracy': 1e-10, #accuracy of magnetic to cartesian coordinat transformation
            'truncation': 1e-10, #trancation of mn harmonics
             'angleStep': 2.0,   #if these parameters are present the 3D grid will be generated to speed up
                  'step': 0.015, #the calculations.
         'grid_accuracy': 0.0005,  #accuracy when generating grid
       'grid_truncation': 2e-6}  #truncation when generating grid 
    """

    def __init__(self, equilibrium_name, mconf_config=None):
        
        if platform.system()=='Windows':
            libname='mconf_matlab64.dll'
        elif platform.system()=='Linux':
            libname = os.path.join(os.path.dirname(__file__),"unix/mconf_matlab64.so")
        elif platform.system()=='Darwin':
            libname= os.path.join(os.path.dirname(__file__),"unix/mconf_matlab64.dylib")
    
        self.mconf = self.import_mconf(libname=libname) #import mconf (VMEC equilibrium)
    
        # setting up equilibrium
        self.equi_data = self.mconf.MCload(equilibrium_name.encode('utf-8')) #read VMEC equilibrium
        if mconf_config is not None:
            if 'B0' in mconf_config.keys():
                self.mconf.MCsetB0(self.equi_data,ct.c_double(mconf_config['B0']), ct.c_double(mconf_config['B0_angle']))#renormilize B field to B0 T and toroidal angle 
            
            if 'extraLCMS' in mconf_config.keys():
                self.mconf.MCsetsmax(self.equi_data,ct.c_double(mconf_config['extraLCMS']))
            
            if 'step' in mconf_config.keys() and 'angleStep' in mconf_config.keys():
                if 'grid_truncation' in mconf_config.keys():
                    self.mconf.MCtruncate(self.equi_data,ct.c_double(mconf_config['grid_truncation']))
                elif 'truncation' in mconf_config.keys():
                    self.mconf.MCtruncate(self.equi_data,ct.c_double(mconf_config['truncation']))
               
                if 'grid_accuracy' in mconf_config.keys():
                    self.mconf.MCsetAccuracy(self.equi_data,ct.c_double(mconf_config['grid_accuracy']))
                elif 'accuracy' in mconf_config.keys():
                    self.mconf.MCsetAccuracy(self.equi_data,ct.c_double(mconf_config['accuracy'])) #trancate fourier harmonics
                ecode = self.mconf.MCcreateMeshUsingSymmetry( self.equi_data, ct.c_double(mconf_config['step']), ct.c_double(mconf_config['step']), ct.c_double(mconf_config['angleStep']/180*np.pi))
            
            if 'truncation' in mconf_config.keys():
                self.mconf.MCtruncate(self.equi_data,ct.c_double(mconf_config['truncation'])) #trancate fourier harmonics
            
            if 'accuracy' in mconf_config.keys():
                self.mconf.MCsetAccuracy(self.equi_data,ct.c_double(mconf_config['accuracy'])) #trancate fourier harmonics
            
           
        
        #self.mconf.MCtruncate(self.equi_data,ct.c_double(1.e-15)) #trancate fourier harmonics
        #self.mconf.MCsetAccuracy(self.equi_data,ct.c_double(1.e-13)) #trancate fourier harmonics
        
        self.get_s_B_T         = np.vectorize(self.get_s_and_B,signature='(),(),()->(),(n)')
        self.get_M3D_s_B_T     = np.vectorize(self.M3D_get_B2,signature='(),(),()->(),(n)')
        self.get_iota          = np.vectorize(self.MCiota)
        self.get_iota_prime    = np.vectorize(self.MCiotaPrime)
        self.get_toroidal_flux = np.vectorize(self.MCToroidalFlux)
        self.get_reff          = np.vectorize(self.MCreff)
        self.get_Ip            = np.vectorize(self.MCIp)
        self.get_It            = np.vectorize(self.MCIt)
        self.get_Volume        = np.vectorize(self.MCVolume)
        self.get_norm_pol_flux = np.vectorize(self.MCtorFlux2polFlux)
        self.mag2xyz           = np.vectorize(self.MCmag2xyz,signature='(),(),()->(),(),()')
        self.xyz2mag           = np.vectorize(self.MCxyz2mag,signature='(),(),()->(),(),()')
        self.get_s_and_B       = np.vectorize(self.get_s_and_B,signature='(),(),()->(),(n)')

        #retcode = self.mconf.MCgetRayIntersectionPoints(self.equi_data,r0,rd,entry,exit1)
        #retcode = self.mconf.M3DgetRayEntryPoint(self.equi_data,r0,rd,entry)
        self.B0 = None
        self.get_B = self.get_B_vmec
        self.grad_B_grad_s = self.grad_B_grad_s_vmec
    
    def import_mconf(self,libname='mconf_matlab64.so',path='.'):
        mconf = np.ctypeslib.load_library(libname,path)
        vec3  = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='CONTIGUOUS')
        
        # setup the return typs and argument types
        mc_handle = ct.c_longlong
        mconf.MCload.restype = mc_handle
        mconf.MCload.argtypes = [ct.c_char_p]
        
        mconf.MCcreateMeshUsingSymmetry.restype = ct.c_int
        mconf.MCcreateMeshUsingSymmetry.argtypes = [mc_handle,ct.c_double,ct.c_double,ct.c_double]
        
        mconf.MCgetRayIntersectionPoints.restype = ct.c_int
        mconf.MCgetRayIntersectionPoints.argtypes = [mc_handle,vec3,vec3,vec3,vec3]
        mconf.M3DgetRayEntryPoint.restype = ct.c_int
        mconf.M3DgetRayEntryPoint.argtypes = [mc_handle,vec3,vec3,vec3]
        
        mconf.MCgetB00.restype = ct.c_double
        mconf.MCgetB00.argtypes = [mc_handle]
        
        mconf.MCsetB0.restype = None
        mconf.MCsetB0.argtypes = [mc_handle,ct.c_double,ct.c_double]
        
        mconf.MCtruncate.restype = None
        mconf.MCtruncate.argtypes = [mc_handle,ct.c_double]
    
        mconf.MCsetAccuracy.restype = None
        mconf.MCsetAccuracy.argtypes = [mc_handle,ct.c_double]
    
        mconf.MCsetsmax.restype = None
        mconf.MCsetsmax.argtypes = [mc_handle,ct.c_double]
        
        mconf.MCgetBxyz.restype = ct.c_double
        mconf.MCgetBxyz.argtypes = [mc_handle,vec3,vec3]
    
        mconf.M3DgetBxyz.restype = ct.c_double
        mconf.M3DgetBxyz.argtypes = [mc_handle,vec3,vec3]
    
        mconf.M3DgetdB_Gradsxyz.restype = ct.c_double
        mconf.M3DgetdB_Gradsxyz.argtypes = [mc_handle,vec3,vec3,vec3,vec3,vec3,vec3]

        mconf.MCgetdB_Gradsxyz.restype  = ct.c_double
        mconf.MCgetdB_Gradsxyz.argtypes  = [mc_handle,vec3,vec3,vec3,vec3,vec3,vec3]

        mconf.MCgetBandGradientsxyz.restype = None
        mconf.MCgetBandGradientsxyz.argtypes = [mc_handle,vec3,vec3,vec3,vec3,vec3,vec3]

        mconf.MCmag2xyz.restype = None
        mconf.MCmag2xyz.argtypes = [mc_handle,vec3,vec3]

        #mconf.MCxyz2mag.restype = None
        #mconf.MCxyz2mag.argtypes = [mc_handle,vec3,vec3]

        mconf.MCtorFlux2polFlux.restype = ct.c_double
        mconf.MCtorFlux2polFlux.argtypes = [mc_handle,ct.c_double]

        mconf.MCiota.restype = ct.c_double
        mconf.MCiota.argtypes = [mc_handle,ct.c_double]

        mconf.MCiotaPrime.restype = ct.c_double
        mconf.MCiotaPrime.argtypes = [mc_handle,ct.c_double]

        mconf.MCFlux.restype = ct.c_double
        mconf.MCFlux.argtypes = [mc_handle,ct.c_double]

        mconf.MCreff.restype = ct.c_double
        mconf.MCreff.argtypes = [mc_handle,ct.c_double]

        mconf.MCIp.restype = ct.c_double
        mconf.MCIp.argtypes = [mc_handle,ct.c_double]
        
        mconf.MCIt.restype = ct.c_double
        mconf.MCIt.argtypes = [mc_handle,ct.c_double]

        mconf.MCVolume.restype = ct.c_double
        mconf.MCVolume.argtypes = [mc_handle,ct.c_double]
    
        return mconf

    def set_B_replace(self,B_replace,dBxds=None,dByds=None,dBzds=None):
        self.B_replace = B_replace
        self.get_B = self.get_B_replace
        
        if dBxds is not None:
            self.dBxds = dBxds
            self.dByds = dByds
            self.dBzds = dBzds
            self.grad_B_grad_s = self.grad_B_grad_s_replace
    
    def getRayIntersectionPoints(self,origin,direction):
        entry = np.zeros(3)
        exit  = np.zeros(3)
        code  = self.mconf.MCgetRayIntersectionPoints(self.equi_data,origin,direction,entry,exit)
        #if not code:
        #    print 'Code:', code, '; no intercection'
        return entry,code#,exit __ uncomment if you want exit also
    
    def get_B_vmec(self,X):
        B = np.zeros(3)
        s  = self.mconf.MCgetBxyz(self.equi_data,np.array(X),B)
        return s,B
    
    def get_B_replace(self,X):
        B = np.zeros(3)
        s  = self.mconf.MCgetBxyz(self.equi_data,np.array(X),B)
        return s,self.B_replace(s,X)

    def grad_B_grad_s_vmec(self,X):
        B = np.zeros(3)
        dBdx = np.zeros(3)
        dBdy = np.zeros(3)
        dBdz = np.zeros(3)
        grad_s = np.zeros(3)
        
        s  = self.mconf.MCgetdB_Gradsxyz(self.equi_data,np.array(X),B,dBdx,dBdy,dBdz,grad_s)
        return s,B,np.vstack((dBdx,dBdy,dBdz)).T,grad_s
    
    def grad_B_grad_s_replace(self,X):
        B = np.zeros(3)
        dBdx = np.zeros(3)
        dBdy = np.zeros(3)
        dBdz = np.zeros(3)
        grad_s = np.zeros(3)
        
        s  = self.mconf.MCgetdB_Gradsxyz(self.equi_data,np.array(X),B,dBdx,dBdy,dBdz,grad_s)
        
        #dBdx = np.zeros(3)
        #dBdy = np.zeros(3)
        #dBdz = np.zeros(3)
        #return     s, B, np.array((self.dBxds(s),self.dByds(s),self.dBzds(s)))*grad_s[0],\
        #                 np.array((self.dBxds(s),self.dByds(s),self.dBzds(s)))*grad_s[1],\
        #                 np.array((self.dBxds(s),self.dByds(s),self.dBzds(s)))*grad_s[2],grad_s
        return     s, B, self.dBxds(s) * grad_s, self.dByds(s) * grad_s, self.dBzds(s) * grad_s, grad_s 
    
    def M3D_grad_B_grad_s(self,X):
        B = np.zeros(3)
        dBdx = np.zeros(3)
        dBdy = np.zeros(3)
        dBdz = np.zeros(3)
        grad_s = np.zeros(3)
        
        s  = self.mconf.M3DgetdB_Gradsxyz(self.equi_data,np.array(X),B,dBdx,dBdy,dBdz,grad_s)
        return s,B,dBdx,dBdy,dBdz,grad_s

    def get_grads_s_theta_phi(self,X):
        B = np.zeros(3)
        gradB = np.zeros(3)
        gradS = np.zeros(3)
        gradTh = np.zeros(3)
        gradPh = np.zeros(3)
        
        self.mconf.MCgetBandGradientsxyz(self.equi_data, np.array(X), B, gradB, gradS, gradTh, gradPh)
        return B, gradB, gradS, gradTh, gradPh
   
    def get_s_and_B(self,x,y,z):
        B = np.zeros(3)
        s  = self.mconf.MCgetBxyz(self.equi_data,np.array([x,y,z]),B)
        return s,B

    def M3D_get_B(self,X):
        B = np.zeros(3)
        s = self.mconf.M3DgetBxyz(self.equi_data,np.array(X),B)
        #s  = self.mconf.MCgetBxyz(self.equi_data,np.array(X),B)
        return s,B
    
    def M3D_get_B2(self,x,y,z):
        B = np.zeros(3)
        s = self.mconf.M3DgetBxyz(self.equi_data,np.array([x,y,z]),B)
        return s,B
    
    def get_s_B(self,x,y,z):
        sT,BT = self.get_s_B_T(x,y,z)
        return sT.copy(), BT.copy()
    
    def s(self,x,y,z):
        sT,BT = self.get_s_B_T(x,y,z)
        return sT.copy()
    
    def MCmag2xyz(self,x,y,z):
        xyz = np.zeros(3)
        self.mconf.MCmag2xyz(self.equi_data,np.array([x,y,z]),xyz)
        return xyz[0],xyz[1],xyz[2]

    def MCxyz2mag(self,x,y,z):
        mag = np.zeros(3)
        self.mconf.MCxyz2mag(self.equi_data,np.array([x,y,z]),mag)
        return mag[0],mag[1],mag[2]

    def MCtorFlux2polFlux(self,s):
        sPol = self.mconf.MCtorFlux2polFlux(self.equi_data,s)
        return sPol

    def MCiota(self,s):
        iota = self.mconf.MCiota(self.equi_data,s)
        return iota

    def MCiotaPrime(self,s):
        iotaP = self.mconf.MCiotaPrime(self.equi_data,s)
        return iotaP

    def MCToroidalFlux(self,s):
        return self.mconf.MCFlux(self.equi_data,s)

    def MCreff(self,s):
        return self.mconf.MCreff(self.equi_data,s)

    def MCIp(self,s):
        return self.mconf.MCIp(self.equi_data,s)

    def MCIt(self,s):
        return self.mconf.MCIt(self.equi_data,s)

    def MCVolume(self,s):
        return self.mconf.MCVolume(self.equi_data,s)
