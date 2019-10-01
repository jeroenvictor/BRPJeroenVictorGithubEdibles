#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Requirements: Python3, pandas, matplotlib, PyQt5

import sys
# import numpy as np
import pandas as pd

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt5 import QtCore
from PyQt5.QtCore import QAbstractTableModel, Qt
from PyQt5.QtWidgets import QApplication, QMainWindow
from gui import Ui_MainWindow
from edibles_spectrum import edibles_spectrum as edspec


# Directory spectra are contained in (TODO: sensible way of initialising this)
path2spectra = '../../spectra/'


class PandasModel(QAbstractTableModel):
    """
    Class to populate a table view with a pandas dataframe
    """
    def __init__(self, data, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._data = data

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self._data.iloc[index.row(), index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self._data.columns[col]
        return None

    def flags(self, index):
        flags = super(self.__class__, self).flags(index)
        flags |= QtCore.Qt.ItemIsEditable
        flags |= QtCore.Qt.ItemIsSelectable
        flags |= QtCore.Qt.ItemIsEnabled
        flags |= QtCore.Qt.ItemIsDragEnabled
        flags |= QtCore.Qt.ItemIsDropEnabled
        return flags

    def sort(self, Ncol, order):
        """Sort table by given column number.
        """
        try:
            self.layoutAboutToBeChanged.emit()
            self._data = self._data.sort_values(self._data.columns[Ncol],
                                                ascending=not order)
            self.layoutChanged.emit()
        except Exception as e:
            print(e)


class SelectionModel(QAbstractTableModel):
    """
    Simple class to populate a table view with Selected data
    """
    def __init__(self, list, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.headers = ['Filename']
        self.list = list

    def rowCount(self, parent=None):
        return len(self.list)

    def columnCount(self, parent=None):
        return 1

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                return self.list[index.row()]
        return None

    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[col]
        return None


class MainWindow(QMainWindow, Ui_MainWindow):
    """
    Main GUI window: used to interactively view EDIBLES spectra.

    Command line usage: $ python main.py

    Current functionality:
    - Use filter tab to select EDIBLES spectra
    - Use plot tab to display matplotlib plot of highlighted spectra

    TODO:
    - I/O functions for subselection?
    - Add FITS header/stellar/etc info into new panel below MPL canvas?
    - Implement filter functions for overview table (by starname, etc)
    - Expand overview file with additional parameters? (exptime, S/N, etc)
    - Integrate fitting/science functions into GUI
      (i.e. interactive lambda selection for profile fitting, etc)
    - add more TODO points
    """
    def __init__(self, parent=None):
        # Initialise the Main window
        QMainWindow.__init__(self, parent=parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.selected_data = []

        # Load overview file from /catalog/DR3_obslist.ext.txt
        self.load_overview()
        # Initialise Filter/Data overview
        self.filtertab()
        # Initialise MPL figure and toolbar
        self.add_mpl()

        # Connect Main window buttons to relevant functions
        self.ui.FilterAddButton.clicked.connect(lambda: self.filter_add())
        self.ui.FilterRemButton.clicked.connect(lambda: self.filter_rem())
        self.ui.PlotButton.clicked.connect(lambda: self.update_plot())

    def keyPressEvent(self, event):

        # Delete key used to remove highlighted spectra from selected sidebar
        if event.key() == QtCore.Qt.Key_Delete:
            # Get highlighted rows
            idx = self.ui.SelectedDataTable.selectionModel().selectedRows()

            # Retrieve filename corresponding to highlighted rows
            # ( TODO probably a nicer way to do this?)
            filenames = []
            for idxxx in idx:
                filenames.append(self.selectionmodel.data(
                                 self.selectionmodel.index(idxxx.row(), 0)))

            # Iterate over highlighted files, remove rom selection list
            for filename in filenames:
                if len(self.selected_data):
                    if filename in self.selected_data:
                        self.selected_data.remove(filename)
                        continue

            # Regenerate selection model and tableview
            self.selectionmodel = SelectionModel(self.selected_data)
            self.ui.SelectedDataTable.setModel(self.selectionmodel)

    def load_overview(self):
        # Load obslist overview into pandas frame
        self.overview = pd.read_csv('../catalog/DR3_obslist_ext.txt',
                                    delim_whitespace=True)
        cols = list(self.overview)
        # move the Filename column to end of list using index, pop and insert
        cols.insert(len(cols), cols.pop(cols.index('Filename')))
        # use ix to reorder
        self.overview = self.overview.ix[:, cols]

    def add_mpl(self):
        # Initial MPL figure and toolbar
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self.ui.matplotlib,
                                         coordinates=True)
        # Add MPL canvas and toolbar to Main window
        self.ui.mpl_layout.addWidget(self.toolbar)
        self.ui.mpl_layout.addWidget(self.canvas)
        self.canvas.draw()

    def update_plot(self):
        try:
            # Retrieve highlighted row from tableview
            idx = self.ui.SelectedDataTable.selectionModel().selectedRows()
            self.fig.clf()
            self.ax = self.fig.add_subplot(111)
            for idxxx in idx:
                filename = self.selectionmodel.data(
                           self.selectionmodel.index(idxxx.row(), 0))
                # Load wav and flux of corresponding spectrum
                wav, flux = edspec(path2spectra+filename).GetSpectrum()

                # Refresh and replot figure
                self.ax.plot(wav, flux)
            self.canvas.draw()
        except(AttributeError, IndexError):
            pass
        except(FileNotFoundError):
            print('FileNotFoundError: \"' + filename + '\" not found...')

    def filtertab(self):
        # Initialise overview table from pandas dataframe, and show the dialog
        self.model = PandasModel(self.overview)
        self.ui.FiltertableView.setModel(self.model)

    def filter_add(self):
            idx = self.ui.FiltertableView.selectionModel().selectedRows()
            skiptotal = 0
            for idxxx in idx:
                filename = self.model.data(self.model.index(idxxx.row(), 4))

                if len(self.selected_data):
                    if filename in self.selected_data:
                        skiptotal += 1
                        continue
                self.selected_data.append(filename)
            if len(self.selected_data):
                self.selectionmodel = SelectionModel(self.selected_data)
                self.ui.SelectedDataTable.setModel(self.selectionmodel)

            if skiptotal:
                if skiptotal == len(idx):
                    statustext = str(skiptotal) + ' spectra already added...'
                else:
                    statustext = (str(skiptotal) + ' spectra already added. ' +
                                  str(len(idx)-skiptotal) +
                                  ' new spectra added...')
            elif len(idx):
                statustext = str(len(idx))+' spectra addded...'
            else:
                statustext = ''
            self.ui.statusBar.showMessage(statustext, 3000)

    def filter_rem(self):
            idx = self.ui.FiltertableView.selectionModel().selectedRows()
            for idxxx in idx:
                filename = self.model.data(self.model.index(idxxx.row(), 4))

                if len(self.selected_data):
                    if filename in self.selected_data:
                        self.selected_data.remove(filename)
                        continue

            self.selectionmodel = SelectionModel(self.selected_data)
            self.ui.SelectedDataTable.setModel(self.selectionmodel)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())