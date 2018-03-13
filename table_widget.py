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

from PyQt4 import QtGui,QtCore
import sys, csv
import numpy as np
class TableWidget(QtGui.QTableWidget):
    def __init__(self,parent=None):
        super(TableWidget, self).__init__(parent)
# if you want to use parent's methods and be ugly as I do, use parent =)
        self.parent=parent
        self.clip = QtGui.QApplication.clipboard()

        self.horizontalHeader().setDefaultSectionSize(60)
        self.setMinimumWidth(250)
        self.setMinimumHeight(250)
        self.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Minimum)
        self.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
       

        self.rowOrder=[]
        self.columnOrder=[]


    def buildFromList(self,inList,addHeaders=True):
        '''
        This methods builds table from rectangular 2D list or np.array
        where 1 column and 1 line contain names of lines and columns
        '''
        self.setRowCount(0)
        self.setColumnCount(0)
        if addHeaders:
            for row in inList[1:,0]:
                self.insertRow(self.rowCount())
                self.setVerticalHeaderItem(self.rowCount()-1, QtGui.QTableWidgetItem(row))
            for col in inList[0,1:]:
                self.insertColumn(self.columnCount())
                self.setHorizontalHeaderItem(self.columnCount()-1,QtGui.QTableWidgetItem(col))
            inList=inList[1:,1:]
        #asidning values       
        else:
            for row in inList[:,0]:
                self.insertRow(self.rowCount())
            for col in inList[0,:]:
                self.insertColumn(self.columnCount())
        it = np.nditer(inList, flags=['multi_index'])
        while not it.finished:
            self.setItem(it.multi_index[0],it.multi_index[1],QtGui.QTableWidgetItem(str(it[0])))
            it.iternext()
                      
        self.verticalHeader().setDefaultSectionSize(self.verticalHeader().minimumSectionSize())
    
    def addFromList(self,inList,addHeaders=True):
        numCols=self.columnCount()
        for col in inList[:,0]:
            self.insertColumn(numCols)
        it = np.nditer(inList, flags=['multi_index'])
        while not it.finished:
            self.setItem(it.multi_index[1],it.multi_index[0]+numCols,QtGui.QTableWidgetItem(str(it[0])))
            it.iternext()
            
    def keyPressEvent(self, e):
        if (e.modifiers() & QtCore.Qt.ControlModifier):        
            if e.key() == QtCore.Qt.Key_C:
                self.copySelectionToClipboard()
            elif e.key() == QtCore.Qt.Key_A:
                self.selectAll()
        QtGui.QTableWidget.keyPressEvent(self, e)
        self.parent.keyPressEvent(e)
        
    def contextMenuEvent(self, pos):
        menu = QtGui.QMenu()
        copyAction = menu.addAction("Copy")
        copyAllAction = menu.addAction("Copy All")
        selectAll = menu.addAction("Select All")
        action = menu.exec_(QtGui.QCursor.pos())
        if action == copyAction:
            self.copySelectionToClipboard()
        elif action == copyAllAction:
            self.copySelectionToClipboard(True)
        elif action == selectAll:
            self.selectAll()
    def handleSave(self,path):
        rowLog = range(self.rowCount())
        rowIndx = [self.visualRow(i) for i in rowLog]
        rowVis = [x for (y,x) in sorted(zip(rowIndx,rowLog))]
        
        colLog = range(self.columnCount())
        colIndx = [self.visualColumn(i) for i in colLog]
        colVis = [x for (y,x) in sorted(zip(colIndx,colLog))]
        
    
        with open(unicode(path), 'wb') as stream:
            writer = csv.writer(stream)
            rowdata = []
            rowdata.append("")
            for column in colVis:
                rowdata.append(unicode(self.horizontalHeaderItem(column).text()).encode('utf8'))
            writer.writerow(rowdata) 
            for row in rowVis:
               
                rowdata = []
                rowdata.append(unicode(self.verticalHeaderItem(row).text()).encode('utf8'))
                for column in colVis:
                    
                    item = self.item(row, column)
                    if item is not None:
                        rowdata.append(
                            unicode(item.text()).encode('utf8'))
                    else:
                        rowdata.append('')
                writer.writerow(rowdata)    

    def selectAll(self):    
        self.setRangeSelected(QtGui.QTableWidgetSelectionRange( 0,0, self.rowCount()-1, self.columnCount()-1)
                                , True)
            
              
    def copySelectionToClipboard(self,all=False):
        if all:
            selected = [QtGui.QTableWidgetSelectionRange( 0,0, self.rowCount()-1, self.columnCount()-1)]
        else:
            selected = self.selectedRanges()
        s = ""
        for r in xrange(selected[0].topRow(),selected[0].bottomRow()+1):
            for c in xrange(selected[0].leftColumn(),selected[0].rightColumn()+1):
                try:
                    s += str(self.item(r,c).text()) + "\t"
                except AttributeError:
                    s += "\t"
            s = s[:-1] + "\n" #eliminate last '\t'
            self.clip.setText(s)
            
def main():
    
    app = QtGui.QApplication(sys.argv)
    ex = TableWidget()
    ex.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()    
