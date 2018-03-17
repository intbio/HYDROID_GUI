#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) Grigoriy A. Armeev, 2015
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as·
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License v2 for more details.

# Cheers, Satary.
#

import sys,os,time,operator
import numpy as np
from PyQt5 import QtGui, QtCore, QtWidgets

#from file_source_widget import FileSourceWidget

import matplotlib
matplotlib.use('QT5Agg')

#### Uncomment these lines if building py2exe binary with window output only
import warnings
warnings.simplefilter('ignore')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.mlab as mlab

from state_widget import StateWidget
from lane_menu import LaneMenu

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s
    
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        
        self.states=['Locating peaks','Assigning peaks','Quantification','Results','Export']
        self.currentState=0       
        self.workDir=unicode(QtCore.QDir.currentPath())
        
        self.initUI()
        
        self.show()
    
    def initUI(self):
        '''
        Setting up Interface

        MainWindow
        --------------------------------
        | mainWidget                   |
        | ---------- ----------------- |
        | |laneMenu| |plotWidget     | |
        | |        | |---------------| |
        | |        | ||stateWidget  || |
        | |        | ||canvas       || |
        | |        | ||toolbarWidget|| |
        | ---------- ----------------- |
        --------------------------------
        
        '''
        mainWidget = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        
        self.laneMenu=LaneMenu(self.workDir)

        mainWidget.addWidget(self.laneMenu)
        plotAndControlsWidget=QtWidgets.QWidget(self)
        plotAndControlsWidgetLayout=QtWidgets.QHBoxLayout(plotAndControlsWidget)
        
        
        highlightColor = str('white')#tabs.palette().color(QtWidgets.QPalette.HighlightedText).name())
       
        plotWidget=QtWidgets.QWidget(self)
        plotWidgetLayout=QtWidgets.QVBoxLayout(plotWidget)
        
        self.stateWidget=StateWidget(self.states,self.currentState)
        self.stateWidget.setFixedSize(600,50)
        plotWidgetLayout.addWidget(self.stateWidget)
        
        self.figure = plt.figure(facecolor=highlightColor)
        self.canvas = FigureCanvas(self.figure)
        plotWidgetLayout.addWidget(self.canvas)
        
        toolbarWidget=QtWidgets.QWidget(self)
        toolbarWidgetLayout=QtWidgets.QHBoxLayout(toolbarWidget)
        self.toolbar = NavigationToolbar(self.canvas, self)
        toolbarWidgetLayout.addWidget(self.toolbar)
        
        self.BackBtn = QtWidgets.QPushButton("Back")
        self.BackBtn.clicked.connect(self.back)     
        self.BackBtn.setEnabled(False)   
        toolbarWidgetLayout.addWidget(self.BackBtn)
        
        self.NextBtn = QtWidgets.QPushButton("Next")
        self.NextBtn.clicked.connect(self.next)
        toolbarWidgetLayout.addWidget(self.NextBtn)
        plotWidgetLayout.addWidget(toolbarWidget)
        
        plotAndControlsWidgetLayout.addWidget(plotWidget)
        
        mainWidget.addWidget(plotAndControlsWidget)
        self.setCentralWidget(mainWidget)
        self.setWindowTitle('HYDROID WIZARD')   
        
        
        #self.assignConnections()
        
    def toggleState(self,direction):
        self.currentState += 1 if direction=='next' else -1
        self.currentState = 0 if self.currentState<0 else self.currentState
        self.currentState = len(self.states)-1 if self.currentState>=len(self.states) else self.currentState
        if self.currentState==0:
            self.BackBtn.setEnabled(False)
        elif self.currentState==(len(self.states)-1):
            self.NextBtn.setEnabled(False)
        else:
            self.BackBtn.setEnabled(True)
            self.NextBtn.setEnabled(True)
            
        self.stateWidget.currentState=self.currentState
        self.stateWidget.update()
        
        # put some meat here
        
    def back(self):
        self.toggleState('back')

    def next(self):
        self.toggleState('next')
        
    def closeEvent(self, event):
        print("Running cleanup...")
        if hasattr(self.laneMenu, 'config'):
            os.remove(self.laneMenu.config.configFile)
        
'''      
    def assignConnections(self):
        self.connect(self.laneMenu,QtCore.SIGNAL("updateFilePreview"),self.set_data)
        self.connect(self.laneMenu,QtCore.SIGNAL("updateUI"),self.repaint)
        self.connect(self.laneMenu,QtCore.SIGNAL("runCalculations"),self.plot_multiple_files)
        self.connect(self.settingsWidget,QtCore.SIGNAL("settingsUpdatedSignal"),self.apply_settings)
'''   


def main():
    
    app = QtWidgets.QApplication(sys.argv)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()    
