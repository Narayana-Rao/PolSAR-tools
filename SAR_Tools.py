# -*- coding: utf-8 -*-
"""
/***************************************************************************

                                 A QGIS plugin
This plugin generates derived SAR parameters from input polarimetric matrix (C3, T3, C2, T2).
                              -------------------
        begin                : 2020-02-03
        git sha              : $Format:%H$
        copyright            : (C) 2020 by MRSLab
        email                : bnarayanarao@iitb.ac.in
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtWidgets import *
from qgis.PyQt.QtGui import *

from qgis.PyQt import *
from qgis.core import *
import requests
import numpy as np
import multiprocessing
import webbrowser


from .resources import *
# Import the code for the dialog
from .SAR_Tools_dialog import MRSLabDialog
import os.path
import sys
# import os
# import os
from osgeo import gdal
import time
from PyQt5 import QtWidgets
# import QtCore
#################


from .functions.dp.mod_DpRVI import DpRVI
from .functions.dp.mod_PRVI_dp import PRVI_dp
from .functions.dp.mod_dop_dp import dop_dp
from .functions.dp.mod_RVI_dp import RVIdp

from .functions.fp.mod_NM3CF import NM3CF
from .functions.fp.mod_MF4CF import MF4CF
from .functions.fp.mod_GRVI import GRVI
from .functions.fp.mod_PRVI import PRVI
from .functions.fp.mod_dop_fp import dop_FP
from .functions.fp.mod_RVIFP import RVI_FP

from .functions.dcop.mod_MF3CD import MF3CD

from .functions.cp.mod_dop_cp import dop_cp
from .functions.cp.mod_CpRVI import CpRVI
from .functions.cp.mod_iS_Omega import iS_Omega
from .functions.cp.mod_NM3CC import NM3CC
#############################

# Create a lock for multiprocess
p_lock = multiprocessing.Lock()

############################################################################################################################################
############################################################################################################################################

class MRSLab(object):
    """QGIS Plugin Implementation."""
    sig_abort_workers = pyqtSignal()
    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface


        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'MRSLab_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)
            QCoreApplication.installTranslator(self.translator)

        self.dlg = MRSLabDialog()

        # Declare instance attributes
        self.actions = []
        icon_path = ':/plugins/SAR_Tools/icon.png'
        self.menu = self.tr(u'&PolSAR tools')

        # Check if plugin was started the first time in current QGIS session
        # Must be set in initGui() to survive plugin reloads
        self.first_start = None
        ##################################################################
        # USER VARIABLES
        
        self.inFolder=''
        # self.ws = 5
        
        self.toolbar = self.iface.addToolBar(u'PolSAR Tools')
        self.toolbar.setObjectName(u'PolSAR Tools')
        
        
        # self.dlg.fp_browse.setEnabled(True)
        # self.dlg.inFolder_fp.setEnabled(False)
        self.dlg.inFolder_fp.setEnabled(False)
        self.dlg.fp_browse.setEnabled(False)
        # self.dlg.fp_cb_C3.setEnabled(False)
        # self.dlg.fp_cb_T3.setEnabled(False)
       
        self.dlg.inFolder_cp.setEnabled(False)
        self.dlg.cp_browse.setEnabled(False)
        # self.dlg.cp_cb_C2.setEnabled(False)
        # self.dlg.cp_cb_T2.setEnabled(False)
        
        self.dlg.inFolder_dp.setEnabled(False)
        self.dlg.dp_browse.setEnabled(False)
        # self.dlg.dp_cb_C2.setEnabled(False)
        # self.dlg.dp_cb_T2.setEnabled(False)
                     
        
        self.dlg.fp_ws.setEnabled(True)
        self.dlg.cp_ws.setEnabled(True)
        self.dlg.dp_ws.setEnabled(True)
        self.dlg.fp_parm.setEnabled(True)
        self.dlg.cp_parm.setEnabled(True)
        self.dlg.dp_parm.setEnabled(True)
        # self.dlg.lineEdit.clear()
        
        self.dlg.pb_process.setEnabled(False)
        # self.dlg.cp_cb_tau.setCurrentText('Tau')
        # Set active tab background colour  
        self.dlg.tabWidget.setStyleSheet(
                """
            QTabBar::tab:selected {
                background: rgb(0, 175, 255)
            }
            """
            )
        
        
        
    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('MRSLab', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            # Adds plugin icon to Plugins toolbar
            self.iface.addToolBarIcon(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)
        self.Startup()
        

        """ Global variables"""
        global  logger
        logger = self.dlg.terminal
        global pol_mode
        pol_mode = self.dlg.tabWidget.currentIndex()




        # self.dlg.cb_mat_type.currentIndexChanged.connect(self.Cob_mode)

        # logger.append(str(self.dlg.tabWidget.currentIndex()))
        # self.dlg.fp_parm.currentIndexChanged.connect(self.Cob_mode)
        """CP"""
        # self.dlg.fp_cb_C3.setChecked(False)
        # self.dlg.fp_cb_T3.setChecked(True)
        # self.dlg.fp_cb_T3.stateChanged.connect(self.fpt3_state_changed)
        # self.dlg.fp_cb_C3.stateChanged.connect(self.fpc3_state_changed)

        self.dlg.fp_browse.clicked.connect(self.openRaster)

        self.dlg.fp_parm.currentIndexChanged.connect(self.Cob_parm)   
        self.ws = int(self.dlg.fp_ws.value())

        self.dlg.fp_ws.valueChanged.connect(self.ws_update)
        
        """CP"""
        # self.dlg.cp_cb_C2.setChecked(False)
        # self.dlg.cp_cb_T2.setChecked(True)
        # self.dlg.cp_cb_T2.stateChanged.connect(self.cpt2_state_changed)
        # self.dlg.cp_cb_C2.stateChanged.connect(self.cpc2_state_changed)
        self.dlg.cp_browse.clicked.connect(self.openRaster)
        self.dlg.cp_ws.valueChanged.connect(self.ws_update)
        self.dlg.cp_parm.currentIndexChanged.connect(self.Cob_parm)
        
        self.psi_val=0
        self.chi_val=45
        self.dlg.cp_sb_psi.valueChanged.connect(self.psi_update)
        self.dlg.cp_sb_chi.valueChanged.connect(self.chi_update)
        
        """DP"""
        # self.dlg.dp_cb_C2.setChecked(False)
        # self.dlg.dp_cb_T2.setChecked(True)
        # self.dlg.dp_cb_T2.stateChanged.connect(self.dpt2_state_changed)
        # self.dlg.dp_cb_C2.stateChanged.connect(self.dpc2_state_changed)
        self.dlg.dp_browse.clicked.connect(self.openRaster)
        self.dlg.dp_ws.valueChanged.connect(self.ws_update)
        self.dlg.dp_parm.currentIndexChanged.connect(self.Cob_parm) 


        """ TAB; CLEAR; PROCESS; VIEW """
        self.dlg.tabWidget.currentChanged.connect(self.ontabChange)               
        self.dlg.pb_view.clicked.connect(self.viewData)
        self.dlg.clear_terminal.clicked.connect(self.clear_log)
        self.dlg.pb_process.clicked.connect(self.startProcess)
        self.dlg.help_btn.clicked.connect(lambda: webbrowser.open('https://sar-tools.readthedocs.io/en/latest/'))
        self.dlg.close_btn.clicked.connect(self.closeui_fn)
        # self.dlg.close_btn.clicked.connect(self.cancel_fn)
        # self.dlg.pb_cancel.clicked.connect(self.cancel_fn)
        # self.dlg.pb_cancel.clicked.connect(lambda: self.worker.stop())
        return action
    
    #@pyqtSlot()  
    # Print the tab/polarimetric mode update to the logger 
    def ontabChange(self,i): #changed!
        if i==0:
            # logger.append("->> Full-pol")
            pol_mode = i            
            # logger.append(str(pol_mode))
        if i==1:
            # logger.append("->> Compact-pol")
            pol_mode = i
            # logger.append(str(pol_mode))
        if i==2:
            # logger.append("->> Dual-pol")
            pol_mode = i
            # logger.append(str(pol_mode))
            
   

    def closeui_fn(self):

        self.dlg.close()
            
        
    def psi_update(self):
        self.psi_val = float(self.dlg.cp_sb_psi.value())

    def chi_update(self):
        self.chi_val = float(self.dlg.cp_sb_chi.value())
            
    def ws_update(self):
        
        if self.dlg.tabWidget.currentIndex()==0:
            self.ws = int(self.dlg.fp_ws.value())
        if self.dlg.tabWidget.currentIndex()==1:
            self.ws = int(self.dlg.cp_ws.value())
        if self.dlg.tabWidget.currentIndex()==2:
            self.ws = int(self.dlg.dp_ws.value())
        if self.ws%2==0:
            self.ws+=1
        # logger = self.dlg.terminal
        # logger.append('->> Window size: '+str(self.ws))
        
 
    def dtype_error(self):
        logger.append('->> Error!! Invalid data folder.')
                 
    def startProcess(self):
        
        
        if self.dlg.tabWidget.currentIndex() == 0:
            # self.inFolder = str(QFileDialog.getExistingDirectory(self.dlg, 
                                                            # "Select T3/C3/T2/C2 Folder"))

            # if(self.fp_cb_C3.isChecked()):
            indX =self.dlg.fp_parm.currentIndex()          
            
            if indX==1:
                try:
                    logger.append('->> --------------------')
                    self.startGRVI()
                except:
                    self.dtype_error()
                    
            if indX==2:
                try:
                    logger.append('->> --------------------')
                    self.startNM3CF()
                except:
                    self.dtype_error()
                    
            if indX==3:
                try:
                    logger.append('->> --------------------')
                    self.startPRVI()
                except:
                    self.dtype_error()
                    
            if indX==4:
                try:
                    logger.append('->> --------------------')
                    self.startDOPfp()
                except:
                    self.dtype_error()
            if indX==5:
                try:
                    logger.append('->> --------------------')
                    self.startRVIFP()
                except:
                    self.dtype_error()
            if indX==6:
                try:
                    logger.append('->> --------------------')
                    self.startMF4CF()
                except:
                    self.dtype_error()
            else:
                pass
            
            
        if self.dlg.tabWidget.currentIndex() == 1:
            # if(self.fp_cb_C3.isChecked()):
            indX =self.dlg.cp_parm.currentIndex()  
            if indX==1:
                try:
                    logger.append('->> --------------------')
                    self.startNM3CC()
                except:
                    self.dtype_error()
            
            if indX==2:
                try:
                    logger.append('->> --------------------')
                    self.startDOPCP()
                except:
                    self.dtype_error()
                    
            if indX==3:
                try:
                    logger.append('->> --------------------')
                    self.startCPRVI()
                except:
                    self.dtype_error()

            if indX==4:
                try:
                    logger.append('->> --------------------')
                    self.startiSOmega()
                except:
                    self.dtype_error()

            else:
                pass
        
        if self.dlg.tabWidget.currentIndex() == 2:
            # if(self.fp_cb_C3.isChecked()):
            indX =self.dlg.dp_parm.currentIndex()   
            if indX==1:
                try:
                    logger.append('->> --------------------')
                    self.startDpRVI()
                except:
                    self.dtype_error()

            if indX==2:
                try:
                    logger.append('->> --------------------')
                    self.startPRVIdp()
                except:
                    self.dtype_error()

            if indX==3:
                try:                    
                    logger.append('->> --------------------')
                    self.startDOPdp()
                except:
                    self.dtype_error()
            if indX==4:
                try:                    
                    logger.append('->> --------------------')
                    self.startRVIdp()
                except:
                    self.dtype_error()


            else:
                pass
            
    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/SAR_Tools/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Process'),
            callback=self.run,
            parent=self.iface.mainWindow())

    def Cob_parm(self):
        # For terminal and UI update
        
        if self.dlg.tabWidget.currentIndex() == 0:
            parm =self.dlg.fp_parm.currentIndex()
            if parm == 1:
                # logger.append('->>      GRVI')
                self.dlg.inFolder_fp.setEnabled(True)
                self.dlg.fp_browse.setEnabled(True)
                # self.dlg.fp_cb_T3.setChecked(True)
                # self.dlg.fp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            elif parm == 2:
                # logger.append('->>      MF3CF')
                self.dlg.inFolder_fp.setEnabled(True)
                self.dlg.fp_browse.setEnabled(True)
                # self.dlg.fp_cb_T3.setChecked(True)
                # self.dlg.fp_browse.setEnabled(True)
                # self.dlg.fp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            elif parm == 3:
                # logger.append('->>      PRVI')
                self.dlg.inFolder_fp.setEnabled(True)
                self.dlg.fp_browse.setEnabled(True)
                # self.dlg.fp_cb_T3.setChecked(True)
                # self.dlg.fp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            elif parm == 4:
                # logger.append('->>      DOP')
                self.dlg.inFolder_fp.setEnabled(True)
                self.dlg.fp_browse.setEnabled(True)
                # self.dlg.fp_cb_T3.setChecked(True)
                # self.dlg.fp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            elif parm == 5:
                # logger.append('->>      RVI')
                self.dlg.inFolder_fp.setEnabled(True)
                self.dlg.fp_browse.setEnabled(True)
                # self.dlg.fp_cb_T3.setChecked(True)
                # self.dlg.fp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            elif parm == 6:
                # logger.append('->>      MF4CF')
                self.dlg.inFolder_fp.setEnabled(True)
                self.dlg.fp_browse.setEnabled(True)
                # self.dlg.fp_cb_T3.setChecked(True)
                # self.dlg.fp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            

            elif parm==0:
                self.dlg.inFolder_fp.setEnabled(False)
                self.dlg.pb_process.setEnabled(False)
                self.dlg.fp_browse.setEnabled(False)
                # self.dlg.fp_ws.setEnabled(False)

        if self.dlg.tabWidget.currentIndex() == 1:
            parm =self.dlg.cp_parm.currentIndex()
            # tau = self.dlg.cp_cb_tau.currentIndex()
            if parm == 1:
                # logger.append('->>     MF3CC')
                self.dlg.inFolder_cp.setEnabled(True)
                self.dlg.cp_browse.setEnabled(True)
                # self.dlg.cp_cb_C2.setChecked(True)
                # self.dlg.cp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            
            if parm == 2:
                # logger.append('->>     DOP')
                self.dlg.inFolder_cp.setEnabled(True)
                self.dlg.cp_browse.setEnabled(True)
                # self.dlg.cp_cb_C2.setChecked(True)
                # self.dlg.cp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            
            if parm == 3:
                # logger.append('->>    CpRVI')
                self.dlg.inFolder_cp.setEnabled(True)
                self.dlg.cp_browse.setEnabled(True)
                # self.dlg.cp_cb_C2.setChecked(True)
                # self.dlg.cp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            
            if parm == 4:
                # logger.append('->>    CpRVI')
                self.dlg.inFolder_cp.setEnabled(True)
                self.dlg.cp_browse.setEnabled(True)
                # self.dlg.cp_cb_C2.setChecked(True)
                # self.dlg.cp_ws.setEnabled(True)
                self.dlg.cp_sb_psi.setEnabled(True)
                self.dlg.cp_sb_chi.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            
            elif parm==0:
                self.dlg.inFolder_cp.setEnabled(False)
                self.dlg.pb_process.setEnabled(False)
                self.dlg.cp_browse.setEnabled(False)
                self.dlg.cp_sb_psi.setEnabled(False)
                self.dlg.cp_sb_chi.setEnabled(False)
                # self.dlg.fp_ws.setEnabled(False)
  
        if self.dlg.tabWidget.currentIndex() == 2:
            parm =self.dlg.dp_parm.currentIndex()
            
            if parm == 1:
                # logger.append('->>      DpRVI')
                # self.dlg.dp_cb_C2.setChecked(True)
                self.dlg.inFolder_dp.setEnabled(True)
                self.dlg.dp_browse.setEnabled(True)
                # self.dlg.dp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            if parm == 2:
                # logger.append('->>      PRVI')
                # self.dlg.dp_cb_C2.setChecked(True)
                self.dlg.inFolder_dp.setEnabled(True)
                self.dlg.dp_browse.setEnabled(True)
                # self.dlg.dp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            if parm == 3:
                # logger.append('->>      DOP')
                # self.dlg.dp_cb_C2.setChecked(True)
                self.dlg.inFolder_dp.setEnabled(True)
                self.dlg.dp_browse.setEnabled(True)
                # self.dlg.dp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            if parm == 4:
                # logger.append('->>      RVI')
                # self.dlg.dp_cb_C2.setChecked(True)
                self.dlg.inFolder_dp.setEnabled(True)
                self.dlg.dp_browse.setEnabled(True)
                # self.dlg.dp_ws.setEnabled(True)
                self.dlg.pb_process.setEnabled(True)
            elif parm==0:
                self.dlg.inFolder_dp.setEnabled(False)
                self.dlg.dp_browse.setEnabled(False)
                self.dlg.pb_process.setEnabled(False)
                # self.dlg.dp_ws.setEnabled(False)
           
            
            
            
    def viewData(self):
        # log_text = self.dlg.terminal
        # log_text.append('->> Data loaded in to QGIS\n')
        
        file_filter = "bin (*.bin);;GeoTiFF (*.tif);;All (*.*)"
                                          
        if self.inFolder:
            f_path = self.inFolder
        else:
            f_path = os.path.dirname(__file__)
        names = QFileDialog.getOpenFileNames(self.dlg, 
                                            "Select files to view/import into QGIS",
                                            f_path,
                                            file_filter       
                                                      )
        
        if names is not None:
            for i in np.arange(0,np.size(list(names[0][0:])),1):
                try:

                    self.iface.addRasterLayer(str(names[0][i]))   
                    logger.append(str(names[0][i]))

                except:
                    logger.append("->> invalid file type!!")

        # logger.append(str(np.size(list(names[0][0:]))))
        # logger.append(str(f_path))    
        
        
  
            
    def clear_log(self):
        self.dlg.terminal.clear()
        self.Startup()
        # self.dlg.cb_mat_type.setCurrentIndex(0)
        self.dlg.inFolder_fp.clear()
        self.dlg.inFolder_cp.clear()
        self.dlg.inFolder_dp.clear()
        # self.dlg.inFolder_fp.setEnabled(False)
        # self.dlg.fp_browse.setEnabled(False)
        self.dlg.progressBar.setValue(0)
        self.dlg.fp_ws.setValue(5)
        self.dlg.cp_sb_psi.setValue(0)
        self.dlg.cp_sb_chi.setValue(45)
        # self.dlg.fp_ws.setEnabled(False)
        self.dlg.fp_parm.setCurrentIndex(0)
        self.dlg.cp_ws.setValue(5)
        # self.dlg.cp_ws.setEnabled(False)
        self.dlg.cp_parm.setCurrentIndex(0)
        self.dlg.dp_ws.setValue(5)
        # self.dlg.dp_ws.setEnabled(False)
        self.dlg.dp_parm.setCurrentIndex(0)
        
        self.dlg.pb_process.setEnabled(False)
        
        # log = self.dlg.terminal
        # log.append('MRS Lab\n')
        
    def showmsg(self, signal):
        log = self.dlg.terminal
        log.append(str(signal))  
        
    def T3_C3(self,T3_stack):
        nrows = np.size(T3_stack,0)
        ncols = np.size(T3_stack,1)
        C3_stack = np.zeros(np.shape(T3_stack),dtype=np.complex64)
        "Special Unitary Matrix"
        D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
        for i in range(nrows):
            # self.dlg.terminal.append('>>> '+str(i)+'/'+str(nrows))
            self.dlg.progressBar.setValue(int((i/nrows)*100))
            for j in range(ncols):
                T3 = T3_stack[i,j,:]
                T3 = np.reshape(T3,(3,3))
                C3 = np.matmul(np.matmul((D.T),T3),D);
                C3_stack[i,j,:] = C3.flatten()
                
        self.dlg.progressBar.setValue(0)
        return C3_stack
    
    def C3_T3(self,C3_stack):
        nrows = np.size(C3_stack,0)
        ncols = np.size(C3_stack,1)
        T3_stack = np.zeros(np.shape(C3_stack),dtype=np.complex64)
        "Special Unitary Matrix"
        D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
        for i in range(nrows):
            # self.dlg.terminal.append('>>> '+str(i)+'/'+str(nrows))
            self.dlg.progressBar.setValue(int((i/nrows)*100))
            for j in range(ncols):
                C3 = C3_stack[i,j,:]
                C3 = np.reshape(C3,(3,3))
                T3 = np.matmul(np.matmul((D),C3),D.T);
                T3_stack[i,j,:] = T3.flatten()
        self.dlg.progressBar.setValue(0)
        return T3_stack
    
    def showTip(self):
       msgBox = QMessageBox()
       msgBox.setIcon(QMessageBox.Information)
       msgBox.setText("Please select a valid matrix folder \
                      \n generated from PolSARpro \
                      \n file format: *.bin, *.hdr")
       msgBox.setWindowTitle("Tip!")
       msgBox.setStandardButtons(QMessageBox.Ok)
       # msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
       # msgBox.buttonClicked.connect(msgButtonClick)
    
       returnValue = msgBox.exec()
       # if returnValue == QMessageBox.Ok:
          # print('OK clicked')

    def showError3(self):
       msgBox = QMessageBox()
       msgBox.setIcon(QMessageBox.Information)
       msgBox.setText("Please select a valid C3/T3 matrix folder \
                      \n generated from PolSARpro \
                      \n file format: *.bin, *.hdr")
       msgBox.setWindowTitle("Error!")
       msgBox.setStandardButtons(QMessageBox.Ok)
       # msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
       # msgBox.buttonClicked.connect(msgButtonClick)
    
       returnValue = msgBox.exec()
       # if returnValue == QMessageBox.Ok:
          # print('OK clicked')

    def showError2(self):
       msgBox = QMessageBox()
       msgBox.setIcon(QMessageBox.Information)
       msgBox.setText("Please select a valid C2 matrix folder \
                      \n generated from PolSARpro \
                      \n file format: *.bin, *.hdr")
       msgBox.setWindowTitle("Error!")
       msgBox.setStandardButtons(QMessageBox.Ok)
       # msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
       # msgBox.buttonClicked.connect(msgButtonClick)
    
       returnValue = msgBox.exec()
       # if returnValue == QMessageBox.Ok:
          # print('OK clicked')


    def openRaster(self):
        """Open raster from file dialog"""
        # logger.append(str(self.dlg.tabWidget.currentIndex()))
        self.showTip() # pop-up tip
        pB = self.dlg.progressBar
        pB.setValue(0)
        if self.dlg.tabWidget.currentIndex() == 0:
            self.inFolder = str(QFileDialog.getExistingDirectory(
                            self.dlg, "Select T3/C3 Folder"))                   
            self.dlg.inFolder_fp.setText(self.inFolder)
            
            if self.inFolder:
                try:
                    self.T3_stack = self.load_T3(self.inFolder)
                    logger.append('->> Ready to process. Click the "process" button!')
                    
                except:
                    try:
                        self.C3_stack = self.load_C3(self.inFolder)
                        logger.append('->> C3 Loaded \n->> Converting C3 to T3...')

                        self.T3_stack  = self.C3_T3(self.C3_stack)
                        
                        logger.append('->> Ready to process. Click the "process" button!')
                        self.dlg.fp_ws.setEnabled(True)
                        self.dlg.fp_parm.setEnabled(True)
                    except:
                        logger.append('->> Error! \n->> Please select a valid C3/T3 folder')
                        self.showError3()
            
            # if self.inFolder:
            #     self.dlg.fp_ws.setEnabled(True)
            #     self.dlg.fp_parm.setEnabled(True)
                
                
        if self.dlg.tabWidget.currentIndex() == 1:
            self.inFolder = str(QFileDialog.getExistingDirectory(
                            self.dlg, "Select C2 Folder"))
            self.dlg.inFolder_cp.setText(self.inFolder)
            
            if self.inFolder:
                try:
                    # logger.append('->> C2 selected')
                    self.C2_stack = self.load_C2(self.inFolder)
                    logger.append('->> Ready to process.')
                    self.dlg.cp_ws.setEnabled(True)
                    self.dlg.cp_parm.setEnabled(True)
                except:
                    logger.append('->> Error! \n->> Please select a valid C2 folder')
                    self.showError2()
            
            # if self.inFolder:


                
        if self.dlg.tabWidget.currentIndex() == 2:
            self.inFolder = str(QFileDialog.getExistingDirectory(
                            self.dlg, "Select C2 Folder"))                   
            self.dlg.inFolder_dp.setText(self.inFolder)
            
            if self.inFolder:
                try:
                    logger.append('->> C2 selected')
                    self.C2_stack = self.load_C2(self.inFolder)
                    logger.append('->> Ready to process.')
                    self.dlg.dp_ws.setEnabled(True)
                    self.dlg.dp_parm.setEnabled(True)
                except:
                    logger.append('->> Error! \n->> Please select a valid C2 folder')
                    self.showError2()

        pB.setValue(0)

            
###############################################################
    def load_C2(self,folder):
    
        C11 = self.read_bin(folder+"/C11.bin")
        C22 = self.read_bin(folder+"/C22.bin")
    
        C12_i = self.read_bin(folder+'/C12_imag.bin')
        C12_r = self.read_bin(folder+'/C12_real.bin')
    
        C12 = C12_r + 1j*C12_i
    
        return np.dstack((C11,C12,np.conj(C12),C22))

    def load_C3(self,folder):
        
        C11 = self.read_bin(folder+"/C11.bin")
        C22 = self.read_bin(folder+"/C22.bin")
        C33 = self.read_bin(folder+"/C33.bin")

        C12_i = self.read_bin(folder+'/C12_imag.bin')
        C12_r = self.read_bin(folder+'/C12_real.bin')
        C13_i = self.read_bin(folder+'/C13_imag.bin')
        C13_r = self.read_bin(folder+'/C13_real.bin')
        C23_i = self.read_bin(folder+'/C23_imag.bin')
        C23_r = self.read_bin(folder+'/C23_real.bin')
            
        C12 = C12_r + 1j*C12_i
        C13 = C13_r + 1j*C13_i
        C23 = C23_r + 1j*C23_i
        
        return np.dstack((C11,C12,C13,np.conj(C12),C22,C23,np.conj(C13),np.conj(C23),C33))
    
    
    def load_T3(self,folder):
        
        T11 = self.read_bin(folder+"/T11.bin")
        T22 = self.read_bin(folder+"/T22.bin")
        T33 = self.read_bin(folder+"/T33.bin")

        T12_i = self.read_bin(folder+'/T12_imag.bin')
        T12_r = self.read_bin(folder+'/T12_real.bin')
        T13_i = self.read_bin(folder+'/T13_imag.bin')
        T13_r = self.read_bin(folder+'/T13_real.bin')
        T23_i = self.read_bin(folder+'/T23_imag.bin')
        T23_r = self.read_bin(folder+'/T23_real.bin')
            
        T12 = T12_r + 1j*T12_i
        T13 = T13_r + 1j*T13_i
        T23 = T23_r + 1j*T23_i
        
        return np.dstack((T11,T12,T13,np.conj(T12),T22,T23,np.conj(T13),np.conj(T23),T33))
    
    def load_T2(self,folder):
    
        T11 = self.read_bin(folder+"/T11.bin")
        T22 = self.read_bin(folder+"/T22.bin")
    
        T12_i = self.read_bin(folder+'/T12_imag.bin')
        T12_r = self.read_bin(folder+'/T12_real.bin')
    
        T12 = T12_r + 1j*T12_i
    
        return np.dstack((T11,T12,np.conj(T12),T22))

    def read_bin(self,file):
        ds = gdal.Open(file)
        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray()
   
        return arr 

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&PolSAR tools'),
                action)
            self.iface.removeToolBarIcon(action)


    def run(self):
        """Run method that performs all the real work"""

        # Create the dialog with elements (after translation) and keep reference
        # Only create GUI ONCE in callback, so that it will only load when the plugin is started
        if self.first_start == True:
            self.first_start = False
            self.dlg = MRSLabDialog()
        
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            # Do something useful here - delete the line containing pass and
            # substitute with your code.
            pass
        
    def Startup(self):
        # For terminal outputs
        logger = self.dlg.terminal
        logger.append('\t                       Welcome to PolSAR tools.'+
                      '\n\t This plugin generates derived SAR parameters')
        logger.append('\t     SAR indices | Decomposition parameters')        
        logger.append('---------------------------------------------------------------------------------------------')
        logger.append('Note: Start by selecting a parameter from "select parameter" dropdown\n')
        """ Process button calls"""

    def startPRVIdp(self):
        
        self.dlg.terminal.append('->> Calculating PRVI... ')
        worker = PRVI_dp(self.inFolder,self.C2_stack,self.ws)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.kill

    def startDOPdp(self):
        
        self.dlg.terminal.append('->> Calculating DOP... ')
        worker = dop_dp(self.inFolder,self.C2_stack,self.ws)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.kill
            
    def startDOPfp(self):  
        self.dlg.terminal.append('->> Calculating DOP...')
        worker = dop_FP(self.inFolder,self.T3_stack,self.ws)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.


    def startPRVI(self):  
        self.dlg.terminal.append('->> Calculating PRVI...')
        worker = PRVI(self.inFolder,self.T3_stack,self.ws)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.

    def startiSOmega(self):   
        
        self.dlg.terminal.append('->> Calculating iS-Omega powers...')
        tau = self.dlg.cp_cb_tau.currentIndex()
            
        worker = iS_Omega(self.inFolder,self.C2_stack,self.ws,tau,self.psi_val,self.chi_val)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.

    def startCPRVI(self):        
        self.dlg.terminal.append('->> Calculating CpRVI...')
        tau = self.dlg.cp_cb_tau.currentIndex()
            
        worker = CpRVI(self.inFolder,self.C2_stack,self.ws,tau)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.


    def startDOPCP(self):        
        self.dlg.terminal.append('->> Calculating DOP...')
        tau = self.dlg.cp_cb_tau.currentIndex()
            
        worker = dop_cp(self.inFolder,self.C2_stack,self.ws,tau)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.
        
    def startNM3CC(self):
        
        self.dlg.terminal.append('->> Calculating MF3CC...')
        tau = self.dlg.cp_cb_tau.currentIndex()
            
        worker = NM3CC(self.inFolder,self.C2_stack,self.ws,tau)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.
    
    def startMF4CF(self):
        
        self.dlg.terminal.append('->> Calculating MF4CF...')
        tau = self.dlg.cp_cb_tau.currentIndex()
            
        worker = MF4CF(self.inFolder,self.T3_stack,self.ws)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.
        
    def startDpRVI(self):
        
        self.dlg.terminal.append('->> Calculating DpRVI... ')
        worker = DpRVI(self.inFolder,self.C2_stack,self.ws)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.kill

    def startRVIdp(self):
        
        self.dlg.terminal.append('->> Calculating RVI... ')
        worker = RVIdp(self.inFolder,self.C2_stack,self.ws)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.kill
            
    def startNM3CF(self):
        
        self.dlg.terminal.append('->> Calculating MF3CF...')
        worker = NM3CF(self.inFolder,self.T3_stack,self.ws)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.kill

    def startRVIFP(self):
        
        self.dlg.terminal.append('->> Calculating RVI...')
        worker = RVI_FP(self.inFolder,self.T3_stack,self.ws)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.kill
        
    def startGRVI(self):
        
        self.dlg.terminal.append('->> Calculating GRVI...')
        
        self.__threads = []
        worker = GRVI(self.inFolder,self.T3_stack,self.ws)

        # messageBar = self.iface.messageBar().createMessage('Doing something time consuming...', )
        # progressBar = QtGui.QProgressBar()
        # progressBar.setAlignment(QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        # cancelButton = QtGui.QPushButton()
        # cancelButton.setText('Cancel')
        # cancelButton.clicked.connect(worker.kill)
        # messageBar.layout().addWidget(progressBar)
        # messageBar.layout().addWidget(cancelButton)
        # self.iface.messageBar().pushWidget(messageBar, self.iface.messageBar().INFO)
        # self.messageBar = messageBar

        # start the worker in a new thread
        thread = QtCore.QThread()
        self.__threads.append((thread, worker))
        worker.moveToThread(thread)
        # self.workerFinished =1


        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)

        
        thread.started.connect(worker.run)
        thread.start()
        # self.sig_abort_workers.connect(worker.abort)

        self.thread = thread
        self.worker = worker

        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)
        # time.sleep(0.1)
        # worker.kill

    def startMF3CD(self):
        
        self.dlg.terminal.append('->> Calculating MF3CD...')
        worker = MF3CD(self.inFolder,self.T2_stack,self.ws)

        # start the worker in a new thread
        thread = QtCore.QThread()
        worker.moveToThread(thread)
        # self.workerFinished =1
        worker.finished.connect(self.workerFinished)
        worker.error.connect(self.workerError)

        worker.progress.connect(self.showmsg)
        worker.pBar.connect(self.pBarupdate)
        thread.started.connect(worker.run)
        thread.start()
        
        self.thread = thread
        self.worker = worker
        # time.sleep(0.1)
        # worker.kill

    def cancel_fn(self):
        # self.sig_abort_workers.emit()
        self.dlg.terminal.append('->> cancelling...')
        for thread, worker in self.__threads:  # note nice unpacking by Python, avoids indexing
            thread.quit()  # this will quit **as soon as thread event loop unblocks**
            thread.wait()  # <- so you need to wait for it to *actually* quit
            # worker.kill()


    def pBarupdate(self, signal):
        pB = self.dlg.progressBar
        pB.setValue(int(signal))
        # log.append(str(signal))  

    def workerFinished(self,finish_cond):

        # if finish_cond:
        logger = self.dlg.terminal
        logger.append('->> Process completed with ' +str(self.ws)+' x ' +str(self.ws)+' window ')
        # clean up the worker and thread
    
        # self.viewData() # Load data into QGIS
        #Open output folder after finishing the process
        path = os.path.realpath(self.inFolder)
        os.startfile(path)

        #set progress bar to Zero
        pB = self.dlg.progressBar
        pB.setValue(100)

        self.worker.deleteLater()
        self.thread.quit()
        self.thread.wait()
        self.thread.deleteLater()

        if finish_cond == 0:
            # self.worke
            logger.append('->> Process stopped in between ! You are good to go again.')
            pB.setValue(0)

    def workerError(self, e, exception_string):
        logger = self.dlg.terminal
        logger.append('->> :-( Error:\n\n %s' %str(exception_string))
    
class UserAbortedNotification(Exception):
    pass 
        
        
        
        
        
        
        
        
