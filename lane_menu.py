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
import sys, os
import pandas as pd
from PyQt5 import QtGui, QtCore, QtWidgets
import matplotlib.pyplot as plt
import numpy as np
from time import sleep

class LaneMenu(QtWidgets.QWidget):
    '''
    Provides Widget for opening multiple files
    '''
    def __init__(self,workDir,parent=None):
        super(LaneMenu, self).__init__(parent)
        self.setAcceptDrops(True)
        self.parent=parent
        self.filters="Excel files (*.xls)"
        self.laneList=[]
        self.widgetWithActiveThread=None

        self.foldersScrollArea = QtWidgets.QScrollArea(self)
        self.foldersScrollArea.setWidgetResizable(True)

        self.foldersScrollAreaWidget = QtWidgets.QWidget()
        self.foldersScrollAreaWidget.setGeometry(QtCore.QRect(0, 0, 380, 280))
        self.folderLayout = QtWidgets.QGridLayout(self.foldersScrollAreaWidget)
        self.folderLayout.setAlignment(QtCore.Qt.AlignTop)
        self.foldersScrollArea.setWidget(self.foldersScrollAreaWidget)
        
        headerWidget=QtWidgets.QWidget()
        headerWidgetLayout = QtWidgets.QGridLayout(headerWidget)
        headerWidgetLayout.setSpacing(0)
        
        headerWidgetLayout.setContentsMargins(0,0,0,0)
        headerWidgetLayout.setAlignment(QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        
        hlabel=QtWidgets.QLabel('Label')
        hlabel.setFixedWidth(50)
        hname=QtWidgets.QLabel('Name')
        hname.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        hstrand=QtWidgets.QLabel('Strand')
        hstrand.setFixedWidth(50)
        href=QtWidgets.QLabel('Ref.')
        href.setFixedWidth(50)
        
        headerWidgetLayout.addWidget(hlabel,0,0)
        headerWidgetLayout.addWidget(hname,0,1)
        headerWidgetLayout.addWidget(hstrand,0,2)
        headerWidgetLayout.addWidget(href,0,3)
        
        self.folderLayout.addWidget(headerWidget)
        
        
        openFiles = QtWidgets.QPushButton("Open lane profiles")
        openFiles.clicked.connect(self.openFile)
 
        self.mainLayout = QtWidgets.QVBoxLayout(self)
        self.mainLayout.addWidget(openFiles)
        self.mainLayout.addWidget(self.foldersScrollArea)
        self.setMaximumWidth(300)
        self.setGeometry(300, 200, 300, 400)
    
    def openFile(self):
        #VERRY CRUDE MUST FIX
        
        filename=QtWidgets.QFileDialog.getOpenFileName(self,'Open Lane profiles',filter=self.filters)
        # returns tupple like ('name', 'filter')
        # that means that filename string is not empty as '' == False in boolean conversion
        if filename[0]:
            self.removeAll()
            profdf=pd.read_csv(filename[0],delimiter="\t",engine='python')
            for label in list(profdf.columns.values)[1::2]:
                self.laneList.append(singleLaneWidget(label))
                self.laneList[-1].killAllThreadsSignal.connect(self.runSingleWidgetThread)
                self.folderLayout.addWidget(self.laneList[-1])
    
    def runSingleWidgetThread(self,lanewidget):
        if self.widgetWithActiveThread == None:
            lanewidget.startThread()
        else:
            self.widgetWithActiveThread.stopThread() 
            self.widgetWithActiveThread=lanewidget
            sleep(0.5)
            self.widgetWithActiveThread.startThread() 
            
            
                
    def removeAll(self):
        for i in reversed(range(len(self.laneList))):
                self.laneList[i].setParent(None) 
                self.laneList.pop(i)
 


class singleLaneWidget(QtWidgets.QWidget):
    killAllThreadsSignal = QtCore.pyqtSignal(QtWidgets.QWidget)
    def __init__(self,label,name=None):
        super(singleLaneWidget, self).__init__()
        self.label=label
        if name is not None:
            self.name=name
        else:
            self.name=label
        
        
        self.Layout = QtWidgets.QGridLayout(self)
        self.Layout.setSpacing(0)
        self.Layout.setContentsMargins(0,0,0,0)
        self.Layout.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTop)
        
        
        self.r_button = QtWidgets.QPushButton(self.label)
        #self.r_button.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.r_button.setFixedWidth(50)
        self.r_button.setStyleSheet("text-align: left;padding: 3px")    
        self.r_button.clicked.connect(self.checkOtherThreads)
        
        self.laneNameWidget=QtWidgets.QLineEdit(self.label)
        
        self.strandCB=QtWidgets.QComboBox()
        self.strandCB.addItems(["TS", "BS"])
        self.strandCB.setFixedWidth(50)

        self.refLaneCheckBox = QtWidgets.QCheckBox(self)
        self.refLaneCheckBox.setTristate(False)
        self.refLaneCheckBox.setChecked(False)
        self.refLaneCheckBox.setFixedWidth(50)
        
        
        self.Layout.addWidget(self.r_button,0,0)
        self.Layout.addWidget(self.laneNameWidget,0,1)
        self.Layout.addWidget(self.strandCB,0,2)
        self.Layout.addWidget(self.refLaneCheckBox,0,3)
        self.thread=display_thread()
        
        
        
    def stopThread(self):
        print self.label
        #self.thread.stop()
        if hasattr(self, 'thread'):
            if self.thread.isRunning:
                self.thread.close()
                self.thread.quit()
                self.thread.wait()
    
    def checkOtherThreads(self):
        self.killAllThreadsSignal.emit(self)

    def startThread(self):
    #    #THAT IS HORRIBLE, we need
        sleep(0.1)
        self.thread.start()
                


class display_thread(QtCore.QThread):

    def __init__(self):
        QtCore.QThread.__init__(self)

    def __del__(self):
        self.wait()
    def close(self):
        plt.close()
        
    def run(self):
        
        t = np.arange(0.0, 2.0, 0.01)
        s1 = np.sin(2*np.pi*t)
        s2 = np.sin(4*np.pi*t)

        plt.figure(1)
        plt.subplot(211)
        plt.plot(t, s1)
        plt.subplot(212)
        plt.plot(t, 2*s1)

        plt.show()
        # your logic here
        
def main():
    
    app = QtWidgets.QApplication(sys.argv)
    workDir=unicode(QtCore.QDir.currentPath())
    ex = LaneMenu(workDir)
    ex.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()    

