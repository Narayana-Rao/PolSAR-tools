# -*- coding: utf-8 -*-
"""
/***************************************************************************
 MRSLab
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

from .resources import *
# Import the code for the dialog
from .SAR_Tools_dialog import MRSLabDialog
import os.path
from osgeo import gdal
import time

#################
from .mod_DpRVI import DpRVI
from .mod_NM3CF import NM3CF
from .mod_NM3CC import NM3CC
from .mod_GRVI import GRVI
from .mod_PRVI import PRVI
#############################

# Create a lock for multiprocess
p_lock = multiprocessing.Lock()

############################################################################################################################################
############################################################################################################################################

class MRSLab(object):
    """QGIS Plugin Implementation."""

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
        self.menu = self.tr(u'&SAR tools')

        # Check if plugin was started the first time in current QGIS session
        # Must be set in initGui() to survive plugin reloads
        self.first_start = None
        ##################################################################
        # USER VARIABLES
        
        self.inFolder=''
        # self.ws = 5
        
        self.toolbar = self.iface.addToolBar(u'SAR Tools')
        self.toolbar.setObjectName(u'SAR Tools')
        
        
        self.dlg.pb_browse.setEnabled(False)
        self.dlg.le_infolder.setEnabled(False)
        self.dlg.sb_ws.setEnabled(False)
        # self.dlg.lineEdit.clear()

        
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
        global  logger
               
        
        self.dlg.pb_browse.clicked.connect(self.openRaster)
        self.dlg.pb_view.clicked.connect(self.viewData)
        self.dlg.clear_terminal.clicked.connect(self.clear_log)
        
        self.dlg.cb_mat_type.currentIndexChanged.connect(self.Cob_mode)
        
        self.ws = int(self.dlg.sb_ws.value())
        
        self.dlg.sb_ws.valueChanged.connect(self.ws_update)

        
        # self.dlg.cb_parm.currentIndexChanged.connect(self.Cob_mode)
        self.dlg.cb_parm.currentIndexChanged.connect(self.Cob_parm)            
        self.dlg.pb_process.clicked.connect(self.startProcess)
        
        return action
            
               
    def ws_update(self):
        self.ws = int(self.dlg.sb_ws.value())
        if self.ws%2==0:
            self.ws+=1
        logger = self.dlg.terminal
        logger.append('Window size: '+str(self.ws))
        
        
    def startProcess(self):
        logger = self.dlg.terminal
        indX =self.dlg.cb_parm.currentIndex()
        if indX==1:
            logger.append('--------------------')
            self.startGRVI()
        if indX==2:
            logger.append('--------------------')
            self.startNM3CF()
        if indX==3:
            logger.append('--------------------')
            self.startNM3CC()
        if indX==4:
            logger.append('--------------------')
            self.startDpRVI()
        if indX==5:
            logger.append('--------------------')
            self.startPRVI()
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
        # For terminal outputs
        logger = self.dlg.terminal
        # global mat_type
        global parm
        parm =self.dlg.cb_parm.currentIndex()
        if parm == 1:
            logger.append('\n     GRVI\n')
            # self.dlg.le_infolder.setEnabled(True)
            # self.dlg.pb_browse.setEnabled(True)
            # self.dlg.sb_ws.setEnabled(True)
            self.dlg.pb_process.setEnabled(True)
        elif parm == 2:
            logger.append('\n     NM3CF\n')
            # self.dlg.le_infolder.setEnabled(True)
            # self.dlg.pb_browse.setEnabled(True)
            # self.dlg.sb_ws.setEnabled(True)
            self.dlg.pb_process.setEnabled(True)
        elif parm == 3:
            logger.append('\n     NM3CC\n')
            # self.dlg.le_infolder.setEnabled(True)
            # self.dlg.pb_browse.setEnabled(True)
            # self.dlg.sb_ws.setEnabled(True)
            self.dlg.pb_process.setEnabled(True)
        elif parm == 4:
            logger.append('\n     DpRVI\n')
            # self.dlg.le_infolder.setEnabled(True)
            # self.dlg.pb_browse.setEnabled(True)
            # self.dlg.sb_ws.setEnabled(True)
            self.dlg.pb_process.setEnabled(True)
        elif parm == 5:
            logger.append('\n     PRVI\n')
            # self.dlg.le_infolder.setEnabled(True)
            # self.dlg.pb_browse.setEnabled(True)
            # self.dlg.sb_ws.setEnabled(True)
            self.dlg.pb_process.setEnabled(True)
        else:
            # self.dlg.le_infolder.setEnabled(False)
            self.dlg.pb_process.setEnabled(False)
            # self.dlg.sb_ws.setEnabled(False)
            
  
        

    def Cob_mode(self):
        # For terminal outputs
        logger = self.dlg.terminal
        global mat_type

        mat_type=self.dlg.cb_mat_type.currentIndex()
        if mat_type == 1:
            logger.append('\n     Selected Matrix type: T3\n')
            self.dlg.le_infolder.setEnabled(True)
            self.dlg.pb_browse.setEnabled(True)
            self.dlg.sb_ws.setEnabled(True)
        elif mat_type == 2:
            logger.append('\n     Selected Matrix type: C3\n')
            self.dlg.le_infolder.setEnabled(True)
            self.dlg.pb_browse.setEnabled(True)
            self.dlg.sb_ws.setEnabled(True)
        elif mat_type == 3:
            logger.append('\n     Selected Matrix type: C2\n')
            self.dlg.le_infolder.setEnabled(True)
            self.dlg.pb_browse.setEnabled(True)
            self.dlg.sb_ws.setEnabled(True)
        else:
            self.dlg.le_infolder.setEnabled(False)
            self.dlg.pb_browse.setEnabled(False)
            self.dlg.sb_ws.setEnabled(False)
            
            
            
            
    def viewData(self):
        log_text = self.dlg.terminal
        log_text.append('Data loaded in to QGIS\n')
        
        mat_type = self.dlg.cb_mat_type.currentIndex()
        
        if self.inFolder is not None and mat_type==1:
            
            # if os.path.isfile(self.inFolder+'\T11.bin'):
            #     self.iface.addRasterLayer(self.inFolder+'\T11.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\T22.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\T33.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\T12_imag.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\T12_real.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\T13_imag.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\T13_real.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\T23_imag.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\T23_real.bin')

            if os.path.isfile(self.inFolder+'\RVI.bin'):
                self.iface.addRasterLayer(self.inFolder+'\RVI.bin')
            if os.path.isfile(self.inFolder+'\GRVI.bin'):
                self.iface.addRasterLayer(self.inFolder+'\GRVI.bin')
                
        elif self.inFolder is not None and mat_type==2:
            # if os.path.isfile(self.inFolder+'\C11.bin'):
            #     self.iface.addRasterLayer(self.inFolder+'\C11.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C22.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C33.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C12_imag.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C12_real.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C13_imag.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C13_real.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C23_imag.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C23_real.bin')
            if os.path.isfile(self.inFolder+'\RVI.bin'):
                self.iface.addRasterLayer(self.inFolder+'\RVI.bin')
            if os.path.isfile(self.inFolder+'\GRVI.bin'):
                self.iface.addRasterLayer(self.inFolder+'\GRVI.bin')
                
        elif self.inFolder is not None and mat_type==3:
            # if os.path.isfile(self.inFolder+'\C11.bin'):
            #     self.iface.addRasterLayer(self.inFolder+'\C11.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C22.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C12_imag.bin')
            #     self.iface.addRasterLayer(self.inFolder+'\C12_real.bin')
                
            if os.path.isfile(self.inFolder+'\DpRVI.bin'):
                self.iface.addRasterLayer(self.inFolder+'\DpRVI.bin')
            if os.path.isfile(self.inFolder+'\RVI_dp.bin'):
                self.iface.addRasterLayer(self.inFolder+'\RVI_dp.bin')
            
    def clear_log(self):
        self.dlg.terminal.clear()
        self.Startup()
        self.dlg.cb_mat_type.setCurrentIndex(0)
        self.dlg.le_infolder.clear()
        self.dlg.le_infolder.setEnabled(False)
        self.dlg.pb_browse.setEnabled(False)
        self.dlg.progressBar.setValue(0)
        self.dlg.sb_ws.setEnabled(False)
        self.dlg.cb_parm.setCurrentIndex(0)
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
                
        self.dlg.progressBar.setValue(100)
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
        self.dlg.progressBar.setValue(100)
        return T3_stack
    
    def openRaster(self):
        """Open raster from file dialog"""
        
        self.inFolder =str(QFileDialog.getExistingDirectory(self.dlg, "Select T3/C3 Folder"))
        self.dlg.le_infolder.setText(self.inFolder)
        # print(self.inFolder)
        mat_type = self.dlg.cb_mat_type.currentIndex()
        logger = self.dlg.terminal
        
        if self.inFolder is not None and mat_type==1:
            
            self.T3_stack = self.load_T3(self.inFolder)
            # logger.append('>>> T3 Loaded \nConverting T3 to C3...')
            # self.C3_stack  = self.T3_C3(self.T3_stack)
            logger.append('>>> Ready to process.')
            
        elif self.inFolder is not None and mat_type==2:

            logger.append('>>> C3 selected')
            self.C3_stack = self.load_C3(self.inFolder)
            logger.append('>>> C3 Loaded \nConverting C3 to T3...')
            self.T3_stack  = self.C3_T3(self.C3_stack)
            logger.append('>>> Ready to process.')
            
        elif self.inFolder is not None and mat_type==3:
            logger.append('>>> C2 selected')
            self.C2_stack = self.load_C2(self.inFolder)
            logger.append('>>> Ready to process.')
            
            
            
            
###############################################################
    def load_C2(self,folder):
    
        C11 = self.read_bin(folder+"\C11.bin")
        C22 = self.read_bin(folder+"\C22.bin")
    
        C12_i = self.read_bin(folder+'\C12_imag.bin')
        C12_r = self.read_bin(folder+'\C12_real.bin')
    
        C12 = C12_r + 1j*C12_i
    
        return np.dstack((C11,C12,np.conj(C12),C22))

    def load_C3(self,folder):
        
        C11 = self.read_bin(folder+"\C11.bin")
        C22 = self.read_bin(folder+"\C22.bin")
        C33 = self.read_bin(folder+"\C33.bin")

        C12_i = self.read_bin(folder+'\C12_imag.bin')
        C12_r = self.read_bin(folder+'\C12_real.bin')
        C13_i = self.read_bin(folder+'\C13_imag.bin')
        C13_r = self.read_bin(folder+'\C13_real.bin')
        C23_i = self.read_bin(folder+'\C23_imag.bin')
        C23_r = self.read_bin(folder+'\C23_real.bin')
            
        C12 = C12_r + 1j*C12_i
        C13 = C13_r + 1j*C13_i
        C23 = C23_r + 1j*C23_i
        
        return np.dstack((C11,C12,C13,np.conj(C12),C22,C23,np.conj(C13),np.conj(C23),C33))
    
    
    def load_T3(self,folder):
        
        T11 = self.read_bin(folder+"\T11.bin")
        T22 = self.read_bin(folder+"\T22.bin")
        T33 = self.read_bin(folder+"\T33.bin")

        T12_i = self.read_bin(folder+'\T12_imag.bin')
        T12_r = self.read_bin(folder+'\T12_real.bin')
        T13_i = self.read_bin(folder+'\T13_imag.bin')
        T13_r = self.read_bin(folder+'\T13_real.bin')
        T23_i = self.read_bin(folder+'\T23_imag.bin')
        T23_r = self.read_bin(folder+'\T23_real.bin')
            
        T12 = T12_r + 1j*T12_i
        T13 = T13_r + 1j*T13_i
        T23 = T23_r + 1j*T23_i
        
        return np.dstack((T11,T12,T13,np.conj(T12),T22,T23,np.conj(T13),np.conj(T23),T33))
    
    def read_bin(self,file):
        ds = gdal.Open(file)
        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray()
   
        return arr 

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&SAR tools'),
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
        logger.append('   \tWelcome to SARtools.\n\
                      A QGIS plugin to calculate SAR derived parameters \n')
        logger.append('       Let\'s get started --> Start with Selecting polarimetric matrix\n')
        

    def startPRVI(self):  
        self.dlg.terminal.append('>>> Calculating PRVI...')
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


    def startNM3CC(self):
        
        self.dlg.terminal.append('>>> Calculating NM3CC...')
        worker = NM3CC(self.inFolder,self.C2_stack,self.ws)

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
        
        self.dlg.terminal.append('>>> Calculating DpRVI... ')
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
            
    def startNM3CF(self):
        
        self.dlg.terminal.append('>>> Calculating NM3CF...')
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
        
    def startGRVI(self):
        
        self.dlg.terminal.append('>>> Calculating GRVI...')
        worker = GRVI(self.inFolder,self.T3_stack,self.ws)

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
            
    def pBarupdate(self, signal):
        pB = self.dlg.progressBar
        pB.setValue(int(signal))
        # log.append(str(signal))  

    def workerFinished(self,finish_cond):
        logger = self.dlg.terminal
        logger.append('>>> Process completed with window size of '+str(self.ws))
        # clean up the worker and thread
        
        # self.viewData() # Load data into QGIS
        #Open output folder after finishing the process
        path = os.path.realpath(self.inFolder)
        os.startfile(path)
        #set progress bar to Zero
        pB = self.dlg.progressBar
        pB.setValue(0)
        self.worker.deleteLater()
        self.thread.quit()
        self.thread.wait()
        self.thread.deleteLater()

        if finish_cond == 0:
            logger.append('>>> Process stopped in between ! You are good to go again.')

    def workerError(self, e, exception_string):
        logger = self.dlg.terminal
        logger.append('>>> :-( Error:\n\n %s' %str(exception_string))
    
        
        
        
        
        
        
        
        
        
