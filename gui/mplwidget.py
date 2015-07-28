__author__ = 'akalinin'

from PyQt5 import QtWidgets, QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class MplCanvas(FigureCanvas):
    def __init__(self):
        self.figure = Figure()
        self.figure.add_subplot(111)
        super(MplCanvas, self).__init__(self.figure)

    def plot_waveform(self, data):
        ax = self.figure.add_subplot(111)
        ax.hold(False)
        ax.plot(data)

class MplWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(MplWidget, self).__init__(parent)

        self.canvas = MplCanvas()
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)

    @QtCore.pyqtSlot(list)
    def plot_waveform(self, data):
        self.canvas.plot_waveform(data)
        self.canvas.draw()
