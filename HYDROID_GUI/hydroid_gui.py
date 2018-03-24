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
import warnings
warnings.filterwarnings("ignore")
import sys,os
from PyQt5 import QtGui, QtCore, QtWidgets
from state_widget import StateWidget
from lane_menu import LaneMenu

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s
    
class HydroidHui(QtWidgets.QMainWindow):
    def __init__(self):
        super(HydroidHui, self).__init__()
        
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
        | |laneMenu| |fasta widget   | |
        | |        | |---------------| |
        | |        | ||stateWidget  || |
        | |        | ||canvas       || |
        | |        | ||buttonWidget|| |
        | ---------- ----------------- |
        --------------------------------
        
        '''
        mainWidget = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        
        self.laneMenu=LaneMenu(self.workDir,parent=self)

        mainWidget.addWidget(self.laneMenu)
        plotAndControlsWidget=QtWidgets.QWidget(self)
        plotAndControlsWidgetLayout=QtWidgets.QHBoxLayout(plotAndControlsWidget)
        
        
        highlightColor = str('white')#tabs.palette().color(QtWidgets.QPalette.HighlightedText).name())
       
        controlWidget=QtWidgets.QWidget(self)
        controlWidgetLayout=QtWidgets.QVBoxLayout(controlWidget)
        
        self.stateWidget=StateWidget(self.states,self.currentState)
        self.stateWidget.setFixedSize(600,50)
        controlWidgetLayout.addWidget(self.stateWidget)
        
        fastaControlWidget=QtWidgets.QWidget(self)
        fastaControlWidgetLayout=QtWidgets.QHBoxLayout(fastaControlWidget)
        
        fastaTextLabel=QtWidgets.QLabel("Top strand sequence 5' - 3':")
        fastaControlWidgetLayout.addWidget(fastaTextLabel)
        
        fastaTextLabel=QtWidgets.QPushButton("Open FASTA")
        fastaTextLabel.clicked.connect(self.openFastaFile)
        fastaControlWidgetLayout.addWidget(fastaTextLabel)
        
        
        controlWidgetLayout.addWidget(fastaControlWidget)
        
        self.fastawidget = QtWidgets.QTextEdit()
        fixed_font = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont)
        self.fastawidget.setCurrentFont(fixed_font)
        controlWidgetLayout.addWidget(self.fastawidget)
        
        buttonWidget=QtWidgets.QWidget(self)
        buttonWidgetLayout=QtWidgets.QHBoxLayout(buttonWidget)
        
        self.BackBtn = QtWidgets.QPushButton("Back")
        self.BackBtn.clicked.connect(self.back)     
        self.BackBtn.setEnabled(False)   
        buttonWidgetLayout.addWidget(self.BackBtn)
        
        self.NextBtn = QtWidgets.QPushButton("Next")
        self.NextBtn.clicked.connect(self.next)
        buttonWidgetLayout.addWidget(self.NextBtn)
        controlWidgetLayout.addWidget(buttonWidget)
        
        plotAndControlsWidgetLayout.addWidget(controlWidget)
        
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
            self.laneMenu.activateRefButtons(False)
        elif self.currentState==(len(self.states)-1):
            self.NextBtn.setEnabled(False)
        else:
            self.laneMenu.activateRefButtons(True)
            self.BackBtn.setEnabled(True)
            self.NextBtn.setEnabled(True)
        
        if self.laneMenu.laneWidgetList:
            for laneWidget in self.laneMenu.laneWidgetList:
                if laneWidget.workingProcess!=None:
                    laneWidget.workingProcess.plot_process.terminate()
                    
        self.stateWidget.currentState=self.currentState
        self.stateWidget.update()
        
    def openFastaFile(self):
        filename=QtWidgets.QFileDialog.getOpenFileName(self,'Open FASTA',filter="FASTA files (*.fasta)")
        if filename[0]:
            with open(filename[0], 'r') as myfile:
                dataString = myfile.read()
            self.fastawidget.setPlainText(dataString)
    
    def back(self):
        self.toggleState('back')

    def next(self):
        self.toggleState('next')
        
    def closeEvent(self, event):
        print("Running cleanup...")
        if hasattr(self.laneMenu, 'config'):
            os.remove(self.laneMenu.config.configFile)
            for laneWidget in self.laneMenu.laneWidgetList:
                try:
                    os.remove(laneWidget.tempFile)
                except:
                    pass
                if laneWidget.workingProcess!=None:
                    laneWidget.workingProcess.plot_process.terminate()
                if laneWidget.intensities != None:
                    try:
                        os.remove(laneWidget.intensities)
                    except:
                        pass
        


def main():
    
    app = QtWidgets.QApplication(sys.argv)
    ex = HydroidHui()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()    
