
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


class PRVI(QtCore.QObject):
    '''PRVI '''
    def __init__(self,iFolder,T3,ws):
        QtCore.QObject.__init__(self)

        self.iFolder = iFolder
        
        self.T3 = T3
        self.ws=ws
        self.killed = False
        # self.mainObj = MRSLab()
    def run(self):
        finish_cond = 0
        try:
            def PRVI_fn(T3,ws):
        
                DOP = np.zeros([np.size(T3,0),np.size(T3,1)])
                prvi = np.zeros([np.size(T3,0),np.size(T3,1)])
                "Special Unitary Matrix"
                D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
                
                for i in range(np.size(T3,0)):

                    self.pBar.emit((i/np.size(T3,0))*100)
                    for j in range(np.size(T3,1)):
                        det = np.abs(np.linalg.det(T3[i,j,:].reshape((3, 3))))
                        trace = np.abs(np.trace(T3[i,j,:].reshape((3, 3))))
                        if trace==0:
                            DOP[i][j]=0
                        # elif((27*det/trace**3)>1.0):
                        #     DOP[i][j]=0
                        else:
                            DOP[i][j] = np.sqrt(1-((27*det)/trace**3))
                        
                        tempT3 =  np.reshape(T3[i,j,:],(3,3))
                        C3 = np.matmul(np.matmul((D.T),tempT3),D).flatten()
                        
                        prvi[i][j] =(1- DOP[i][j])*C3[4]*0.5 # (1-dop)*vh
            
                """Write files to disk"""
                
                if os.path.exists(self.iFolder+'/T11.bin'):
                    infile = self.iFolder+'/T11.bin'
                elif os.path.exists(self.iFolder+'/C11.bin'):
                    infile = self.iFolder+'/C11.bin'
          
                ofileprvi = self.iFolder+'/PRVI_FP.bin'
                write_bin(ofileprvi,prvi,infile)
                self.pBar.emit(100)
                self.progress.emit('->> Finished PRVI calculation!!')
                


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
        
            PRVI_fn(self.T3,self.ws)
            
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
