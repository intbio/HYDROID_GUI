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
import matplotlib
matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt
import numpy as np
from time import sleep
from hydroid.HYDROIDexp import assign_peaks_interactive
from hydroid_wraper import hydroidConfig,PlotProcess
QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)


class LaneMenu(QtWidgets.QWidget):
    '''
    Provides Widget for opening multiple files
    '''
    def __init__(self,workDir,parent=None):
        super(LaneMenu, self).__init__(parent)
        self.setAcceptDrops(True)
        self.parent=parent
        self.filters="Excel files (*.xls)"
        self.laneWidgetList=[]
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
        
        hlabel=QtWidgets.QLabel('Column')
        hlabel.setFixedWidth(50)
        hname=QtWidgets.QLabel('Name')
        hname.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        hstrand=QtWidgets.QLabel('Strand')
        hstrand.setFixedWidth(50)
        labeledEnd=QtWidgets.QLabel('Label')
        labeledEnd.setFixedWidth(50)
        href=QtWidgets.QLabel('Ref.')
        href.setFixedWidth(50)
        
        headerWidgetLayout.addWidget(hlabel,0,0)
        headerWidgetLayout.addWidget(hname,0,1)
        headerWidgetLayout.addWidget(hstrand,0,2)
        headerWidgetLayout.addWidget(labeledEnd,0,3)
        headerWidgetLayout.addWidget(href,0,4)
        
        self.folderLayout.addWidget(headerWidget)
        
        
        openFiles = QtWidgets.QPushButton("Open lane profiles")
        openFiles.clicked.connect(self.openFile)
 
        self.mainLayout = QtWidgets.QVBoxLayout(self)
        self.mainLayout.addWidget(openFiles)
        self.mainLayout.addWidget(self.foldersScrollArea)
        self.setMaximumWidth(350)
        self.setGeometry(300, 200, 300, 400)
    
    def openFile(self):
        filename=QtWidgets.QFileDialog.getOpenFileName(self,'Open Lane profiles',filter=self.filters)
        # returns tupple like ('name', 'filter')
        # that means that filename string is not empty as '' == False in boolean conversion
        if filename[0]:
            self.removeAll()
            profdf=pd.read_csv(filename[0],delimiter="\t",engine='python')
            #VERRY CRUDE MUST FIX
            laneLabels=list(profdf.columns.values)[1::2]
            #!!!!!!!!!!!!!!!!!!!!
            self.config=hydroidConfig(laneLabels)
            for label in laneLabels:
                self.laneWidgetList.append(singleLaneWidget(label,filename[0],self.config.configFile,mainwindow=self.parent))
                #self.laneWidgetList[-1].killAllThreadsSignal.connect(self.runSingleWidget)
                self.folderLayout.addWidget(self.laneWidgetList[-1])
    
    def runSingleWidget(self,lanewidget):
        if self.widgetWithActiveThread != None:
            self.widgetWithActiveThread.stopThread() 
        self.widgetWithActiveThread=lanewidget
        self.widgetWithActiveThread.startThread() 
            
                
    def removeAll(self):
        for i in reversed(range(len(self.laneWidgetList))):
                self.laneWidgetList[i].setParent(None) 
                self.laneWidgetList.pop(i)
                
    def closeEvent(self, event):
        print("Running cleanup...")
        if hasattr(self, 'config'):
            os.remove(self.config.configFile)

 


class singleLaneWidget(QtWidgets.QWidget):
    killAllThreadsSignal = QtCore.pyqtSignal(QtWidgets.QWidget)
    def __init__(self,label,lane_profile_file,lane_config_file,name=None,mainwindow=None):
        super(singleLaneWidget, self).__init__()
        self.mainwindow=mainwindow
        self.label=label
        self.lane_profile_file=lane_profile_file
        self.lane_config_file=lane_config_file
        self.workingProcess=None
        
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
        self.r_button.clicked.connect(self.runTask)
        
        self.laneNameWidget=QtWidgets.QLineEdit(self.label)
        
        self.strandCB=QtWidgets.QComboBox()
        self.strandCB.addItems(["TS", "BS"])
        self.strandCB.setFixedWidth(50)
        
        self.labeledEnd=QtWidgets.QComboBox()
        self.labeledEnd.addItems(["3`", "5`"])
        self.labeledEnd.setFixedWidth(50)

        self.refLaneCheckBox = QtWidgets.QCheckBox(self)
        self.refLaneCheckBox.setTristate(False)
        self.refLaneCheckBox.setChecked(False)
        self.refLaneCheckBox.setFixedWidth(50)
        
        
        self.Layout.addWidget(self.r_button,0,0)
        self.Layout.addWidget(self.laneNameWidget,0,1)
        self.Layout.addWidget(self.strandCB,0,2)
        self.Layout.addWidget(self.labeledEnd,0,3)
        self.Layout.addWidget(self.refLaneCheckBox,0,4)
        
    
    def runTask(self):
        print self.mainwindow.states[self.mainwindow.currentState]
        if self.workingProcess!=None:
            self.workingProcess.plot_process.terminate()
        if self.mainwindow.currentState==0:
            self.workingProcess=PlotProcess(FUNC='assign_peaks_interactive',lane_profile_file=self.lane_profile_file,
                                lane_config_file=self.lane_config_file,
                                lane_name=self.name)
        elif self.mainwindow.currentState==1:
            self.workingProcess=PlotProcess(FUNC='call_peaks_interactive',lane_profile_file=self.lane_profile_file,
                                lane_config_file=self.lane_config_file,
                                lane_name=self.name,DNAseq=s['seq'],
                                labeled_end=s['label'], helper_prof_names=s['helper_profiles'])


'''
from Bio import SeqIO

lane_profile_file="data/lane_profiles.xls"
lane_config_file="data/lane_config.csv"
#Read DNA seqeunce from file via biopython
TS_seq=SeqIO.parse(open("data/DNA_seq.fasta"),'fasta').next().seq
BS_seq=TS_seq.reverse_complement()
lane_sets=[
{'footprinting_profile':'scCSE4_601TA_BS','helper_profiles':['GA_601TA_BS','CT_601TA_BS'],'seq':BS_seq,'label':'three_prime'},
{'footprinting_profile':'scCSE4_601TA_TS','helper_profiles':['GA_601TA_TS','CT_601TA_TS'],'seq':TS_seq,'label':'three_prime'}
]
call_peaks_interactive(lane_profile_file,lane_config_file,DNAseq=s['seq'],labeled_end=s['label'],lane_name=s['footprinting_profile'],helper_prof_names=s['helper_profiles'])
'''
        
def main():
    
    app = QtWidgets.QApplication(sys.argv)
    workDir=unicode(QtCore.QDir.currentPath())
    ex = LaneMenu(workDir)
    ex.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()    

