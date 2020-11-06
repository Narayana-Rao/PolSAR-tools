
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtWidgets import *
from qgis.PyQt.QtGui import *

from qgis.PyQt import *
from qgis.core import *
import requests
import numpy as np
import multiprocessing

from .resources import *
# Import the code for the dialog
from .SAR_Tools_dialog import MRSLabDialog
import os.path
from osgeo import gdal
import time
import os.path


class DpRVI(QtCore.QObject):
    '''DpRVI '''
    def __init__(self,iFolder,C2,ws):
        QtCore.QObject.__init__(self)

        self.iFolder = iFolder
        
        self.C2 = C2
        self.ws=ws
        self.killed = False
        # self.mainObj = MRSLab()
    def run(self):
        finish_cond = 0
        try:
            def DpRVI_fn(C2_stack,ws):
                """Note:
                ncols=no.of rows 
                nrows=no.of cols"""
                nrows  = np.shape(C2_stack)[1]
                ncols = np.shape(C2_stack)[0]
                
                C11 = C2_stack[:,:,0]
                C12 = C2_stack[:,:,1]
                C21 = C2_stack[:,:,2]
                C22 = C2_stack[:,:,3]
                
                dprvi = np.zeros((ncols,nrows));
                
                dop_b = np.zeros((ncols,nrows));
                fp22 = np.zeros((ncols,nrows));
                rvi = np.zeros((ncols,nrows));
                
                wsi=wsj=ws
                
                inci=int(np.fix(wsi/2)) # Up & down movement margin from the central row
                incj=int(np.fix(wsj/2)) # Left & right movement from the central column
            
                
                starti=int(np.fix(wsi/2)) # Starting row for window processing
                startj=int(np.fix(wsj/2)) # Starting column for window processing
                
                stopi= int(nrows-inci)-1 # Stop row for window processing
                stopj= int(ncols-incj)-1 # Stop column for window processing
            
                for ii in np.arange(startj,stopj+1):
                    
                    self.pBar.emit(int((ii/ncols)*100))
                    for jj in np.arange(starti,stopi+1):
                        
                        C11c = np.nanmean(C11[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        C12c = np.nanmean(C12[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        C21c = np.nanmean(C21[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        C22c = np.nanmean(C22[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
             
                        C0 = np.array([[C11c,C12c], [C21c, C22c]]);
                        
                        # %% GD_VI -- VV-VH/VV-HH
                        if np.isnan(np.real(C0)).any() or np.isinf(np.real(C0)).any() or np.isneginf(np.real(C0)).any():
                            C0 = np.array([[0,0],[0,0]])
                            
                        e_v = -np.sort(-np.linalg.eigvals(C0)); # sorting in descending order
                        e_v1 = e_v[0]; e_v2 = e_v[1];
                        x_1 = e_v1/(e_v1 + e_v2);
                        fp22[ii,jj] = x_1;
                        
                         #%% dop Barakat-DPRVI
                        dop_b[ii,jj] = np.real(np.sqrt(1 - (4*np.linalg.det(C0)/(np.trace(C0)**2))));
                        span_c = C11c + C22c;
                        
                        dprvi[ii,jj] = 1 - dop_b[ii,jj]*(fp22[ii,jj]); #% c22c for S1/c11c for others
                        
                        rvi[ii, jj] = 4*C11c/span_c; #%RVI dual-pol for comparison for S1
                        
                
                """Write files to disk"""
                
                infile = self.iFolder+'/C11.bin'
                
                ofiledprvi = self.iFolder+'/DpRVI.bin'
                write_bin(ofiledprvi,dprvi,infile)
                
                ofilervidp = self.iFolder+'/RVI_dp.bin'
                write_bin(ofilervidp,rvi,infile)
                
                self.pBar.emit(100)
                self.progress.emit('>>> Finished DpRVI calculation!!')
                
            
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
        
            DpRVI_fn(self.C2,self.ws)
            
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
