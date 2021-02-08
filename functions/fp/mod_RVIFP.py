

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

class RVI_FP(QtCore.QObject):
    '''RVI FP '''
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
            def RVIFP_fn(T3_stack,ws):
        
                t11_T1 = T3_stack[:,:,0]
                t12_T1 = T3_stack[:,:,1]
                t13_T1 = T3_stack[:,:,2]
                t21_T1 = T3_stack[:,:,3]
                t22_T1 = T3_stack[:,:,4]
                t23_T1 = T3_stack[:,:,5]
                t31_T1 = T3_stack[:,:,6]
                t32_T1 = T3_stack[:,:,7]
                t33_T1 = T3_stack[:,:,8]
                
                nrows  = np.shape(T3_stack)[1]
                ncols = np.shape(T3_stack)[0]
                # nrows  = 100
                # ncols = 100
               
                temp_rvi = np.zeros((ncols,nrows))
                
                # %% for window processing
                wsi=wsj=ws
                
                inci=int(np.fix(wsi/2)) # Up & down movement margin from the central row
                incj=int(np.fix(wsj/2)) # Left & right movement from the central column
                # % Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999
                
                starti=int(np.fix(wsi/2)) # Starting row for window processing
                startj=int(np.fix(wsj/2)) # Starting column for window processing
                
                stopi= int(nrows-inci)-1 # Stop row for window processing
                stopj= int(ncols-incj)-1 # Stop column for window processing
                        
                for ii in np.arange(startj,stopj+1):
        
                    # self.progress.emit(str(ii)+'/'+str(nrows))
                    self.pBar.emit(int((ii/ncols)*90))
                    for jj in np.arange(starti,stopi+1):
                
                        t11s = np.nanmean(t11_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        t12s = np.nanmean(t12_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        t13s = np.nanmean(t13_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        
                        t21s = np.nanmean(t21_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        t22s = np.nanmean(t22_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        t23s = np.nanmean(t23_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        
                        t31s = np.nanmean(t31_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        t32s = np.nanmean(t32_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        t33s = np.nanmean(t33_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                
                        T_T1 = np.array([[t11s, t12s, t13s], [t21s, t22s, t23s], [t31s, t32s, t33s]])
                
                        # %% RVI
                        if np.isnan(np.real(T_T1)).any() or np.isinf(np.real(T_T1)).any() or np.isneginf(np.real(T_T1)).any():
                            T_T1 = np.array([[0,0],[0,0]])
                            temp_rvi[ii,jj] = 0
                            # self.progress.emit(str('invalid Value encountered!!'))
                            continue
                            
                        e_v = -np.sort(-np.linalg.eigvals(T_T1)); # sorting in descending order
                        e_v1 = e_v[0]; e_v2 = e_v[1]; e_v3 = e_v[2];

                        # self.progress.emit(str('Eigen val Done'))
                        
                        p1 = e_v1/(e_v1 + e_v2 + e_v3);
                        p2 = e_v2/(e_v1 + e_v2 + e_v3);
                        p3 = e_v3/(e_v1 + e_v2 + e_v3);
                        
                        p1=0 if p1<0 else p1
                        p2=0 if p2<0 else p2
                        p3=0 if p3<0 else p3
                            
                        p1=1 if p1>1 else p1
                        p2=1 if p2>1 else p2
                        p3=1 if p3>1 else p3
                        
                        
                        
                        temp_rvi[ii,jj] = np.real((4*p3)/(p1 + p2 + p3));
                
                # %% RVI scaled (0 - 1)   
                rvi = temp_rvi;   
                idx = np.argwhere(rvi>1)
           
                rvi[idx] = (3/4)*rvi[idx];
                rvi[~idx] = rvi[~idx];
                rvi[rvi==0] = np.NaN
                
                self.progress.emit('->> Write files to disk...')
                """Write files to disk"""
                if os.path.exists(self.iFolder+'/T11.bin'):
                    infile = self.iFolder+'/T11.bin'
                elif os.path.exists(self.iFolder+'/C11.bin'):
                    infile = self.iFolder+'/C11.bin'

                ofilervi = self.iFolder+'/RVI_FP.bin'

                write_bin(ofilervi,rvi,infile)
                self.pBar.emit(100)
                self.progress.emit('->> Finished RVI calculation!!')
            
            
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
        
            RVIFP_fn(self.T3,self.ws)
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

