
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

class CpRVI(QtCore.QObject):
    '''CpRVI '''
    def __init__(self,iFolder,C2,ws,tau):
        QtCore.QObject.__init__(self)

        self.iFolder = iFolder
        self.C2 = C2
        self.ws=ws
        self.tau = tau
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
            def CpRVI_fn(C2_stack,ws):
                
                nrows  = np.shape(C2_stack)[1]
                ncols = np.shape(C2_stack)[0]
                
                C11 = C2_stack[:,:,0]
                C12 = C2_stack[:,:,1]
                C21 = C2_stack[:,:,2]
                C22 = C2_stack[:,:,3]
                
                fp22 = np.zeros((ncols,nrows))
                l_lambda = np.zeros((ncols,nrows))

                if self.tau==0:                    
                    chi_in = -45.0
                else:
                    chi_in = 45.0
                self.pBar.emit(90)

                
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
                        
                        C11c = np.nanmean(C11[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        C12c = np.nanmean(C12[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        C21c = np.nanmean(C21[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
                        C22c = np.nanmean(C22[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
             
                        C0 = np.array([[C11c,C12c], [C21c, C22c]]);
                        
                        # %% GD_VI -- VV-VH/VV-HH
                        if np.isnan(np.real(C0)).any() or np.isinf(np.real(C0)).any() or np.isneginf(np.real(C0)).any():
                            C0 = np.array([[0,0],[0,0]])
                
                        # Stokes Parameter
                        s0 = C11c + C22c;
                        s1 = C11c - C22c;
                        s2 = (C12c + C21c);
        
                        if (chi_in >= 0):
                            s3 = (1j*(C12c - C21c)); # The sign is according to RC or LC sign !!
                        if (chi_in < 0):
                            s3 = -(1j*(C12c - C21c)); # The sign is according to RC or LC sign !!
                        
                        k11 = s0
                        k12 = 0
                        k13 = (1/2)*s2
                        k14 = 0
                        k21 = k12 
                        k22 = 0
                        k23 = 0
                        k24 = s1
                        k31 = k13
                        k32 = k23
                        k33 = 0 
                        k34 = 0;
                        k41 = k14; 
                        k42 = k24; 
                        k43 = k34; 
                        k44 = (1/2)*s3;
                
                        K_T = np.array([[k11,k12,k13,k14], [k21,k22,k23,k24], 
                            [k31, k32, k33, k34], [k41,k42,k43,k44]])       
        
                        # Stokes vector child products
                        SC = ((s0)-(s3))/2;
                        OC = ((s0)+(s3))/2;
                        
                        min_sc_oc = min(SC,OC);
                        max_sc_oc = max(SC,OC);                
                        
                        K_depol = np.array([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]);
                        
                        # GD_DEPOL
                        
                        num1 = np.matmul(K_T.T,K_depol);
                        num = np.trace(num1);
                        den1 = np.sqrt(abs(np.trace(np.matmul(K_T.T,K_T))));
                        den2 = np.sqrt(abs(np.trace(np.matmul(K_depol.T,K_depol))));
                        den = den1*den2;
                        
                        temp_aa = np.real(2*np.arccos(num/den)*180/np.pi);
                        GD_t1_depol = np.real(temp_aa/180);
                        
                                                                       
                        l_lambda[ii,jj] = (3/2)*GD_t1_depol;
                        
                        #GD_VI -- RH-RV/LH-LV
                        
                        fp22[ii,jj] = (min_sc_oc/max_sc_oc);

                
                vi_c = (1 - l_lambda)*np.power(fp22, 2*l_lambda);
                
                self.pBar.emit(90)                        
                
                self.progress.emit('>>> Write files to disk...')
                """Write files to disk"""
                
                infile = self.iFolder+'/C11.bin'
                
                ofilecprvi = self.iFolder+'/CpRVI.bin'
                write_bin(ofilecprvi,vi_c,infile)

                self.pBar.emit(100)
                self.progress.emit('->> Finished CpRVI calculation!!')
            

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
        

            CpRVI_fn(self.C2, self.ws)
            
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
