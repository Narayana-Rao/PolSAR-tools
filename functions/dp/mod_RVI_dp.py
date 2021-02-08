
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtWidgets import *
from qgis.PyQt.QtGui import *

from qgis.PyQt import *
from qgis.core import *
import requests
import numpy as np
import multiprocessing

from ...resources import *
# Import the code for the dialog
from ...SAR_Tools_dialog import MRSLabDialog
import os.path
from osgeo import gdal
import time
import os.path


class RVIdp(QtCore.QObject):
    '''DpRVI '''
    def __init__(self,iFolder,C2,ws):
        QtCore.QObject.__init__(self)

        self.iFolder = iFolder
        
        self.C2 = C2
        self.ws=ws
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
    
    def eig22(self,c2):
        c11 = c2[:,:,0].flatten()
        c12 = c2[:,:,1].flatten()
        c21 = c2[:,:,2].flatten()
        c22 = c2[:,:,3].flatten()
        trace = -(c11+c22)
        det = c11*c22-c12*c21
        # const= 1
        sqdiscr = np.sqrt(trace*trace - 4*det);
        lambda1 = -(trace + sqdiscr)*0.5;
        lambda2 = -(trace - sqdiscr)*0.5;
        
        return lambda1,lambda2

    def run(self):
        finish_cond = 0
        try:
            def RVIdp_fn(C2_stack,ws):
                
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

                c2_det = (c11s*c22s-c12s*c21s)
                c2_trace = c11s+c22s
                rvi = 4*c22s/c2_trace
            
                self.pBar.emit(90)
                self.progress.emit('->> Write files to disk...')
                """Write files to disk"""
                
                infile = self.iFolder+'/C11.bin'
                
                ofilervidp = self.iFolder+'/RVI_dp.bin'
                write_bin(ofilervidp,rvi,infile)
                
                self.pBar.emit(100)
                self.progress.emit('->> Finished RVI calculation!!')        
                
                
            
            def read_bin(file):  
                ds = gdal.Open(file)
                band = ds.GetRasterBand(1)
                arr = band.ReadAsArray()
                return arr
            
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
        
            RVIdp_fn(self.C2,self.ws)
            
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
