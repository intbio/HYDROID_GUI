from PyQt5.QtWidgets import QWidget, QApplication
from PyQt5.QtGui import QPainter, QColor, QBrush, QFont
from PyQt5.QtCore import QRect, Qt

import sys

class StateWidget(QWidget):
    
    def __init__(self,states,currentState):
        super(StateWidget,self).__init__()
        self.states=states
        self.currentState=currentState
        self.initUI()
        
        
    def initUI(self):      

       # self.setGeometry(300, 300, 350, 100)
        self.setWindowTitle('Colours')
        self.show()


    def paintEvent(self, e):
        qp = QPainter()
        qp.begin(self)
        qp.setFont(QFont('Decorative', 10))
        self.drawStates(qp)
        qp.end()

        
    def drawStates(self, qp):
        qp.drawLine(10, 25, 100*len(self.states), 25)
        for i,state in enumerate(self.states):
            if i==self.currentState:
                qp.setBrush(QColor(255, 255, 255))
            else:
                qp.setBrush(QColor(150, 150, 150))
            rect = QRect(5+i*120, 10, 110, 30)
            qp.drawRect(rect)
            qp.drawText(rect, Qt.AlignCenter, state)
            
              
        
if __name__ == '__main__':
    states=['Locating peaks','Assigning peaks','Quantification','Results']
    currentState=0
    app = QApplication(sys.argv)
    ex = StateWidget(states,currentState)
    sys.exit(app.exec_())
