
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

class NM3CC(QtCore.QObject):
    '''NM3CC '''
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
            def NM3CC_fn(C2_stack,ws):
        
                C11 = C2_stack[:,:,0]
                C12 = C2_stack[:,:,1]
                C21 = C2_stack[:,:,2]
                C22 = C2_stack[:,:,3]
                
                nrows  = np.shape(C2_stack)[1]
                ncols = np.shape(C2_stack)[0]
                # nrows  = 100
                # ncols = 100
                
                chi_in = -45 # change this to input from user
               
                theta_CP = np.zeros((ncols,nrows))
                Pd_CP = np.zeros((ncols,nrows))
                Pv_CP = np.zeros((ncols,nrows))
                Ps_CP = np.zeros((ncols,nrows))
                
                # D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
                # %% for window processing
                wsi=wsj=ws
                
                inci=int(np.fix(wsi/2)) # Up & down movement margin from the central row
                incj=int(np.fix(wsj/2)) # Left & right movement from the central column
                # % Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999
                
                starti=int(np.fix(wsi/2)) # Starting row for window processing
                startj=int(np.fix(wsj/2)) # Starting column for window processing
                
                stopi= int(nrows-inci)-1 # Stop row for window processing
                stopj= int(ncols-incj)-1 # Stop column for window processing
                        # %% Elementary targets
                             
                for ii in np.arange(startj,stopj+1):
        
                    # self.progress.emit(str(ii)+'/'+str(nrows))
                    self.pBar.emit(int((ii/ncols)*100))
                    for jj in np.arange(starti,stopi+1):
                
                        C11c = np.nanmean(C11[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        C12c = np.nanmean(C12[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        C21c = np.nanmean(C21[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        C22c = np.nanmean(C22[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                
                        C0 = np.array([[C11c,C12c], [C21c, C22c]]);
                        
                        m1 = np.real(np.sqrt(1-(4*(np.linalg.det(C0)/(np.trace(C0)**2)))))
                        
                        # Stokes Parameter
                        s0 = C11c + C22c;
                        s1 = C11c - C22c;
                        s2 = (C12c + C21c);
                        
                        if (chi_in >= 0):
                            s3 = (1j*(C12c - Bin_C21)); # The sign is according to RC or LC sign !!
        
                        if (chi_in < 0):
                            s3 = -(1j*(Bin_C12 - Bin_C21)); # The sign is according to RC or LC sign !!
                    
                        SC = ((s0)-(s3))/2;
                        OC = ((s0)+(s3))/2;
                        
                        h = (OC-SC)
                        
                        span = C11c + C22c
                        
                        val = ((m1*s0*h))/((SC*OC + (m1**2)*(s0**2)))
                        thet = np.real(np.arctan(val))
                        # thet = np.rad2deg(thet)
                        theta_CP[ii,jj] = np.rad2deg(thet)
                        Ps_CP[ii,jj] = (((m1*(span)*(1+np.sin(2*thet))/2)))
                        Pd_CP[ii,jj] = (((m1*(span)*(1-np.sin(2*thet))/2)))
                        Pv_CP[ii,jj] = (span*(1-m1))
                        
                
                
                """Write files to disk"""
                # ofilervi = self.iFolder+'/RVI.bin'
                infile = self.iFolder+'/C11.bin'
                # write_bin(ofilervi,rvi,infile)
                ofilegrvi = self.iFolder+'/Theta_CP.bin'
                write_bin(ofilegrvi,theta_CP,infile)
                ofilegrvi1 = self.iFolder+'/Pd_CP.bin'
                write_bin(ofilegrvi1,Pd_CP,infile)
                ofilegrvi2 = self.iFolder+'/Ps_CP.bin'
                write_bin(ofilegrvi2,Ps_CP,infile)
                ofilegrvi3 = self.iFolder+'/Pv_CP.bin'
                write_bin(ofilegrvi3,Pv_CP,infile)     
                self.pBar.emit(100)
                self.progress.emit('>>> Finished NM3CC calculation!!')
                # self.pBar.emit(0)
                # self.iface.addRasterLayer(self.inFolder+'\RVI.bin')
                # self.iface.addRasterLayer(self.inFolder+'\GRVI.bin')
                # return rvi,vi 
            
            
            
            
            def read_bin(file):
            
                # data, geodata=load_data(file_name, gdal_driver='GTiff')
                ds = gdal.Open(file)
                band = ds.GetRasterBand(1)
                arr = band.ReadAsArray()
                # [cols, rows] = arr.shape
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
        
            # self.dop_fp(self.T3)
            NM3CC_fn(self.C2,self.ws)
            
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
