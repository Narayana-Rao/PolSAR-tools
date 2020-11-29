
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtWidgets import *
from qgis.PyQt.QtGui import *

from qgis.PyQt import *
from qgis.core import *
import requests
import numpy as np
import multiprocessing

# Initialize Qt resources from file resources.py
from .resources import *
# Import the code for the dialog
from .SAR_Tools_dialog import MRSLabDialog
import os.path
from osgeo import gdal
import time
import os.path


##############################################################################################

class iS_Omega(QtCore.QObject):
    '''DOP CP '''
    def __init__(self,iFolder,C2,ws,tau):
        QtCore.QObject.__init__(self)

        self.iFolder = iFolder
        self.C2 = C2
        self.ws=ws
        self.tau=tau
        self.killed = False
        # self.mainObj = MRSLab()
        
    def conv2d(self,a, f):
        filt = np.zeros(a.shape)
        wspad = int(f.shape[0]/2)
        s = f.shape + tuple(np.subtract(a.shape, f.shape) + 1)
        strd = np.lib.stride_tricks.as_strided
        subM = strd(a, shape = s, strides = a.strides * 2)
        filt_data = np.einsum('ij,ijkl->kl', f, subM)
        filt[wspad:wspad+filt_data.shape[0],wspad:wspad+filt_data.shape[1]] = filt_data
        return filt
    
    def run(self):
        finish_cond = 0
        try:
            def iS_Omega_fn(C2_stack,ws):

                if self.tau==0:                    
                    chi_in = -45.0
                else:
                    chi_in = 45.0

                kernel = np.ones((ws,ws),np.float32)/(ws*ws)
                c11_T1 = C2_stack[:,:,0]
                c12_T1 = C2_stack[:,:,1]
                c21_T1 = C2_stack[:,:,2]
                c22_T1 = C2_stack[:,:,3]
            
                c11_T1r = self.conv2d(np.real(c11_T1),kernel)
                c11_T1i = self.conv2d(np.imag(c11_T1),kernel)
                c11s = c11_T1r+1j*c11_T1i

                c12_T1r = self.conv2d(np.real(c12_T1),kernel)
                c12_T1i = self.conv2d(np.imag(c12_T1),kernel)
                c12s = c12_T1r+1j*c12_T1i
                self.pBar.emit(25)

                c21_T1r = self.conv2d(np.real(c21_T1),kernel)
                c21_T1i = self.conv2d(np.imag(c21_T1),kernel)
                c21s = c21_T1r+1j*c21_T1i


                c22_T1r = self.conv2d(np.real(c22_T1),kernel)
                c22_T1i = self.conv2d(np.imag(c22_T1),kernel)
                c22s = c22_T1r+1j*c22_T1i

                self.pBar.emit(50)

                # c2_det = (c11s*c22s-c12s*c21s)
                # c2_trace = c11s+c22s
                # t2_span = t11s*t22s
                # m1 = np.real(np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))

                # Stokes Parameter
                #s0 = c11s + c22s;
                #s1 = c11s - c22s;
                #s2 = (c12s + c21s);

                #if (chi_in >= 0):
                #    s3 = (1j*(c12s - c21s)); # The sign is according to RC or LC sign !!
                #if (chi_in < 0):
                #    s3 = -(1j*(c12s - c21s)); # The sign is according to RC or LC sign !!
                    
                # Stokes Parameter
                s0 = np.float32(np.real(c11s + c22s))
                s1 = np.float32(np.real(c11s - c22s))
                s2 = np.float32(np.real(c12s + c21s))

                if (chi_in >= 0):
                    s3 = np.float32(np.real(1j*(c12s - c21s))) # The sign is according to RC or LC sign !!
                if (chi_in < 0):
                    s3 = np.float32(np.real(-(1j*(c12s - c21s)))) # The sign is according to RC or LC sign !!
                
                ## Stokes child parameters
                SC = ((s0)-(s3))/2;
                OC = ((s0)+(s3))/2;
                #old_err_state = np.seterr(divide='raise')
                #ignored_states = np.seterr(**old_err_state)
                CPR = np.divide(SC,OC)  ##SC/OC
                
                ##scattered fields    
                dop= np.sqrt(np.power(s1,2) + np.power(s2,2) + np.power(s3,2))/(s0)
                Psi = 0.5*((180/np.pi)*np.arctan2(s2,s1))
                DOCP = (-s3)/(dop*s0);
                Chi = 0.5*((180/np.pi)*np.arcsin(DOCP))
                ##---------------------------------
                psi_in = self.tau
                
                ##---------------------------------
                # Calculating Omega from S-Omega decomposition        
                x1 = np.cos(2*chi_in*np.pi/180)*np.cos(2*psi_in*np.pi/180)*np.cos(2*Chi*np.pi/180)*np.cos(2*Psi*np.pi/180)
                x2 = np.cos(2*chi_in*np.pi/180)*np.sin(2*psi_in*np.pi/180)*np.cos(2*Chi*np.pi/180)*np.sin(2*Psi*np.pi/180)
                x3 = np.abs(np.sin(2*chi_in*np.pi/180)*np.sin(2*Chi*np.pi/180))
                Prec =  dop*(1 + x1 + x2 + x3)
                Prec1 = (1 - dop) + dop*(1 + x1 + x2 + x3)
                omega = (Prec/Prec1)
                
                
                ## Improved S-Omega (i-SOmega powers
        
                if (CPR > 1.0):
                    surface_new = omega*s0 - omega*(1 - omega)*SC
                    double_bounce_new = omega*(1 - omega)*SC   ##depolarized of OC x polarized of SC
            
                elif (CPR < 1.0):
                    surface_new = omega*(1 - omega)*OC    ##depolarized of SC x polarized of OC
                    double_bounce_new = omega*s0 - omega*(1 - omega)*OC
            
                elif (CPR == 1.0):
                    surface_new = omega*OC
                    double_bounce_new = omega*SC
    
                else:
                    surface_new = np.nan
                    double_bounce_new = np.nan
                    diffused_new = np.nan
        
                diffused_new = (1 - omega)*s0; ##diffused scattering
                
                self.pBar.emit(90)                        
                
                self.progress.emit('->> Write files to disk...')
                """Write files to disk"""
                infile = self.iFolder+'/C11.bin'
                
                ofileps = self.iFolder+'/Ps_iSOmega.bin'
                write_bin(ofileps,surface_new,infile)
                
                ofilepd = self.iFolder+'/Pd_iSOmega.bin'
                write_bin(ofilepd,double_bounce_new,infile)
                
                ofilepv = self.iFolder+'/Pv_iSOmega.bin'
                write_bin(ofilepv,diffused_new,infile)
                
                self.pBar.emit(100)
                self.progress.emit('->> Finished i-SOmega power calculation!!')


            

            def write_bin(file,wdata,refData):
                
                ds = gdal.Open(refData)
                [cols, rows] = wdata.shape
            
                driver = gdal.GetDriverByName("ENVI")
                outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
                outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
                outdata.SetProjection(ds.GetProjection())##sets same projection as input
                
                outdata.SetDescription(file)
                outdata.GetRasterBand(1).WriteArray(wdata)
                # outdata.GetRasterBand(1).SetNoDataValue(np.NaN)##if you want these values transparent
                outdata.FlushCache() ##saves to disk!!    
        
            # self.dop_fp(self.T3)
            iS_Omega_fn(self.C2,self.ws)
            
            finish_cond = 1
            
        
        except Exception as e:
            # forward the exception upstream
            self.error.emit(e, traceback.format_exc())
            
        self.finished.emit(finish_cond)
        
    # def kill(self):
    #     self.killed = True
        

        
    """***************************************"""
    finished = QtCore.pyqtSignal(object)
    error = QtCore.pyqtSignal(Exception, str)
    progress = QtCore.pyqtSignal(str)     
    pBar =  QtCore.pyqtSignal(int)             
