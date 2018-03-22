#!/usr/bin/env python
#-*- coding:utf-8 -*-
import csv,io
from shutil import copyfile
import sip
sip.setapi('QString', 2)
sip.setapi('QVariant', 2)

from PyQt5 import QtWidgets, QtCore, QtGui

class csvTable(QtWidgets.QWidget):
    def __init__(self, fileName, parent=None):
        super(csvTable, self).__init__(parent)
        self.clip = QtWidgets.QApplication.clipboard()

        
        #self.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Minimum)
        #self.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        
        self.fileName = fileName

        self.model = QtGui.QStandardItemModel(self)

        self.tableView = QtWidgets.QTableView(self)
        self.tableView.setModel(self.model)
        self.tableView.horizontalHeader().setStretchLastSection(True)
        
        self.tableView.horizontalHeader().setDefaultSectionSize(60)
        self.tableView.setMinimumWidth(500)
        self.tableView.setMinimumHeight(250)
        self.tableView.setSizePolicy(QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Minimum)

        self.pushButtonWrite = QtWidgets.QPushButton(self)
        self.pushButtonWrite.setText("Save File")
        self.pushButtonWrite.clicked.connect(self.on_pushButtonWrite_clicked)

        self.layoutVertical = QtWidgets.QVBoxLayout(self)
        self.layoutVertical.addWidget(self.tableView)
        self.layoutVertical.addWidget(self.pushButtonWrite)
        self.setWindowTitle('Export Intensity file')
        self.loadCsv(self.fileName)

    def loadCsv(self, fileName):
        with open(fileName, "rb") as fileInput:
            for row in csv.reader(fileInput): 
                if ''.join(row)[0] != '#':
                    items = [
                        QtGui.QStandardItem(field)
                        for field in row
                    ]
                    self.model.appendRow(items)

    def writeCsv(self, fileName):
        filename=QtWidgets.QFileDialog.getSaveFileName(self,'Save intensity Lane file',filter="csv. files (*.csv)")
        # returns tupple like ('name', 'filter')
        # that means that filename string is not empty as '' == False in boolean conversion
        if filename[0]:
            try:
                copyfile(fileName,filename[0])
            except IOError:
                print 'Can not write to %s'% filename[0]

    def on_pushButtonWrite_clicked(self):
        self.writeCsv(self.fileName)
    
    def eventFilter(self, source, event):
        if (event.type() == QtCore.QEvent.KeyPress and
            event.matches(QtGui.QKeySequence.Copy)):
            self.copySelection()
            return True
        return super(Window, self).eventFilter(source, event)
    
    def keyPressEvent(self, e):
        if (e.modifiers() & QtCore.Qt.ControlModifier):
            if e.key() == QtCore.Qt.Key_C: #copy
                self.copySelection()
                            
    def copySelection(self,all=False):
        if all:
            self.tableView.selectAll()
        selection = self.tableView.selectedIndexes()
        if selection:
            rows = sorted(index.row() for index in selection)
            columns = sorted(index.column() for index in selection)
            rowcount = rows[-1] - rows[0] + 1
            colcount = columns[-1] - columns[0] + 1
            table = [[''] * colcount for _ in range(rowcount)]
            for index in selection:
                row = index.row() - rows[0]
                column = index.column() - columns[0]
                table[row][column] = index.data()
            stream = io.BytesIO()
            csv.writer(stream, delimiter=',').writerows(table)
            QtWidgets.QApplication.clipboard().setText(stream.getvalue())

    def contextMenuEvent(self, pos):
        menu = QtWidgets.QMenu()
        copyAction = menu.addAction("Copy")
        copyAllAction = menu.addAction("Copy All")
        selectAll = menu.addAction("Select All")
        action = menu.exec_(QtGui.QCursor.pos())
        if action == copyAction:
            self.copySelection()
        elif action == copyAllAction:
            self.copySelection(True)
        elif action == selectAll:
            self.tableView.selectAll()


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName('MyWindow')

    main = csvTable("data/results.csv")
    main.show()

    sys.exit(app.exec_())