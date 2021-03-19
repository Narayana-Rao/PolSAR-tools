
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtWidgets import *
from qgis.PyQt.QtGui import *

from qgis.PyQt import *
from qgis.core import *
import requests
import numpy as np
import multiprocessing

# Initialize Qt resources from file resources.py
from ...resources import *
# Import the code for the dialog
from ...SAR_Tools_dialog import MRSLabDialog
import os.path
from osgeo import gdal
import time
import os.path


##############################################################################################

class MF4CF(QtCore.QObject):
    '''MF4CF '''
    def __init__(self,iFolder,T3,ws):
        QtCore.QObject.__init__(self)

        self.iFolder = iFolder
        self.T3 = T3
        self.ws=ws
        self.killed = False
        
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
            def MF4CF_fn(T3_stack,ws):
        
                t11_T1 = T3_stack[:,:,0]
                t12_T1 = T3_stack[:,:,1]
                t13_T1 = T3_stack[:,:,2]
                t21_T1 = T3_stack[:,:,3]
                t22_T1 = T3_stack[:,:,4]
                t23_T1 = T3_stack[:,:,5]
                t31_T1 = T3_stack[:,:,6]
                t32_T1 = T3_stack[:,:,7]
                t33_T1 = T3_stack[:,:,8]
                
                
                kernel = np.ones((ws,ws),np.float32)/(ws*ws)

            
                t11_T1r = self.conv2d(np.real(t11_T1),kernel)
                t11_T1i = self.conv2d(np.imag(t11_T1),kernel)
                t11s = t11_T1r+1j*t11_T1i

                t12_T1r = self.conv2d(np.real(t12_T1),kernel)
                t12_T1i = self.conv2d(np.imag(t12_T1),kernel)
                t12s = t12_T1r+1j*t12_T1i

                t13_T1r = self.conv2d(np.real(t13_T1),kernel)
                t13_T1i = self.conv2d(np.imag(t13_T1),kernel)
                t13s = t13_T1r+1j*t13_T1i
                self.pBar.emit(25)

                t21_T1r = self.conv2d(np.real(t21_T1),kernel)
                t21_T1i = self.conv2d(np.imag(t21_T1),kernel)
                t21s = t21_T1r+1j*t21_T1i

                t22_T1r = self.conv2d(np.real(t22_T1),kernel)
                t22_T1i = self.conv2d(np.imag(t22_T1),kernel)
                t22s = t22_T1r+1j*t22_T1i

                t23_T1r = self.conv2d(np.real(t23_T1),kernel)
                t23_T1i = self.conv2d(np.imag(t23_T1),kernel)
                t23s = t23_T1r+1j*t23_T1i
                
                self.pBar.emit(35)
                
                t31_T1r = self.conv2d(np.real(t31_T1),kernel)
                t31_T1i = self.conv2d(np.imag(t31_T1),kernel)
                t31s = t31_T1r+1j*t31_T1i

                t32_T1r = self.conv2d(np.real(t32_T1),kernel)
                t32_T1i = self.conv2d(np.imag(t32_T1),kernel)
                t32s = t32_T1r+1j*t32_T1i

                t33_T1r = self.conv2d(np.real(t33_T1),kernel)
                t33_T1i = self.conv2d(np.imag(t33_T1),kernel)
                t33s = t33_T1r+1j*t33_T1i
                
                self.pBar.emit(52)
                
                det_T3 = t11s*(t22s*t33s-t23s*t32s)-t12s*(t21s*t33s-t23s*t31s)+t13s*(t21s*t32s-t22s*t31s)
                trace_T3 = t11s + t22s + t33s
                m1 = np.real(np.sqrt(1-(27*(det_T3/(trace_T3**3)))))

                k11_f = (t11s + t22s + t33s)/2
                k44_f = (-t11s + t22s + t33s)/2
                k14_f = np.imag(t23s)

                trace_T3 = t11s + t22s + t33s
                self.pBar.emit(75)
                

                s0_f = trace_T3
                dop_f = m1

                val1 = (4*dop_f*k11_f*k44_f)/(k44_f**2 - (1 + 4*dop_f**2)*k11_f**2)
                val2 = np.abs(k14_f)/(k11_f)
                
                self.pBar.emit(78)

                theta_f = np.real(np.arctan(val1)) # separation for surface and dbl
                tau_f = np.real(np.arctan(val2)) # separation for helix
                # thet = np.rad2deg(thet)
                theta_FP = np.rad2deg(theta_f)
                tau_FP = np.rad2deg(tau_f)
                self.pBar.emit(81)
                

                pc_f = dop_f*s0_f*(np.sin(2*tau_f))
                pv_f = (1-dop_f)*s0_f
                res_pow = s0_f - (pc_f + pv_f)
                ps_f = (res_pow/2)*(1+np.sin((2*theta_f)))
                pd_f = (res_pow/2)*(1-np.sin((2*theta_f)))
                
                
                
                self.pBar.emit(90)
                self.progress.emit('>>> Writing files to disk...')
                """Write files to disk"""
                if os.path.exists(self.iFolder+'/T11.bin'):
                    infile = self.iFolder+'/T11.bin'
                elif os.path.exists(self.iFolder+'/C11.bin'):
                    infile = self.iFolder+'/C11.bin'

                ofilegrvi0 = self.iFolder+'/Tau_FP_4c.bin'
                write_bin(ofilegrvi0,tau_FP,infile)

                ofilegrvi = self.iFolder+'/Theta_FP_4c.bin'
                write_bin(ofilegrvi,theta_FP,infile)
                
                ofilegrvi1 = self.iFolder+'/Pd_FP_4c.bin'
                write_bin(ofilegrvi1,pd_f,infile)
                
                ofilegrvi2 = self.iFolder+'/Ps_FP_4c.bin'
                write_bin(ofilegrvi2,ps_f,infile)
                
                ofilegrvi3 = self.iFolder+'/Pv_FP_4c.bin'
                write_bin(ofilegrvi3,pv_f,infile)

                ofilegrvi4 = self.iFolder+'/Pc_FP_4c.bin'
                write_bin(ofilegrvi4,pc_f,infile)     
                
                self.pBar.emit(100)
                self.progress.emit('>>> Finished MF4CF calculation!!')
                

            
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
            MF4CF_fn(self.T3,self.ws)
            
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
