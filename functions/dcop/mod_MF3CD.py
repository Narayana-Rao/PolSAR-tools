
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

class MF3CD(QtCore.QObject):
    '''MF3CD '''
    def __init__(self,iFolder,T2,ws):
        QtCore.QObject.__init__(self)

        self.iFolder = iFolder
        self.T2 = T2
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
            def MF3CD_fn(T3_stack,ws):
        
                t11_T1 = T3_stack[:,:,0]
                t12_T1 = T3_stack[:,:,1]
                t21_T1 = T3_stack[:,:,2]
                t22_T1 = T3_stack[:,:,3]
                
                
                kernel = np.ones((ws,ws),np.float32)/(ws*ws)

            
                t11_T1r = self.conv2d(np.real(t11_T1),kernel)
                t11_T1i = self.conv2d(np.imag(t11_T1),kernel)
                t11s = t11_T1r+1j*t11_T1i

                t12_T1r = self.conv2d(np.real(t12_T1),kernel)
                t12_T1i = self.conv2d(np.imag(t12_T1),kernel)
                t12s = t12_T1r+1j*t12_T1i

                t21_T1r = self.conv2d(np.real(t21_T1),kernel)
                t21_T1i = self.conv2d(np.imag(t21_T1),kernel)
                t21s = t21_T1r+1j*t21_T1i

                t22_T1r = self.conv2d(np.real(t22_T1),kernel)
                t22_T1i = self.conv2d(np.imag(t22_T1),kernel)
                t22s = t22_T1r+1j*t22_T1i
                
                self.pBar.emit(52)
                
                det_T2 = t11s*t22s-t12s*t21s
                trace_T2 = t11s + t22s

                m1 = np.real(np.sqrt(1-(4*(det_T2/(trace_T2**2)))))
                h = (t11s - t22s)
                g = t22s
                span = t11s + t22s
                self.pBar.emit(75)
                
                val = (m1*span*h)/(t11s*g+m1**2*span**2);
                thet = np.real(np.arctan(val))
                # thet = np.rad2deg(thet)
                theta_DP = np.rad2deg(thet)
                self.pBar.emit(81)
                
                Ps_DP = (((m1*(span)*(1+np.sin(2*thet))/2)))
                Pd_DP = (((m1*(span)*(1-np.sin(2*thet))/2)))
                Pv_DP = (span*(1-m1))
                
                
                
                self.pBar.emit(90)
                self.progress.emit('>>> Writing files to disk...')
                """Write files to disk"""
                if os.path.exists(self.iFolder+'/T11.bin'):
                    infile = self.iFolder+'/T11.bin'
                elif os.path.exists(self.iFolder+'/C11.bin'):
                    infile = self.iFolder+'/C11.bin'

                
                ofilegrvi = self.iFolder+'/Theta_DP.bin'
                write_bin(ofilegrvi,theta_DP,infile)
                
                ofilegrvi1 = self.iFolder+'/Pd_DP.bin'
                write_bin(ofilegrvi1,Pd_DP,infile)
                
                ofilegrvi2 = self.iFolder+'/Ps_DP.bin'
                write_bin(ofilegrvi2,Ps_DP,infile)
                
                ofilegrvi3 = self.iFolder+'/Pv_DP.bin'
                write_bin(ofilegrvi3,Pv_DP,infile)     
                
                self.pBar.emit(100)
                self.progress.emit('>>> Finished MF3CD calculation!!')
                

            
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
            MF3CD_fn(self.T2,self.ws)
            
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
