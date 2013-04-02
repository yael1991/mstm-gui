#!/usr/bin/python

from mstm_materials import *
from mstm_parameters import *
from mstm_simparser import *
import time
import sys

#PyQt4 libraries
from PyQt4 import QtGui
from PyQt4 import QtCore
from PyQt4 import uic

class GuiWindow(QtGui.QMainWindow):
    
    params = ParameterClass('msinput.inp')
    
    def setParams(self):
        #update the Gui based on values in the parameters structure
        self.ui.spinStartLambda.setValue(self.params.minLambda)
        self.ui.spinEndLambda.setValue(self.params.maxLambda)
        self.ui.spinNumSamples.setValue(self.params.nSamples)
        self.ui.spinNumSpheres.setValue(int(self.params['number_spheres']))
        #self.ui.spinAlpha.setValue(float(self.params['incident_azimuth_angle_deg']))
        #self.ui.spinBeta.setValue(float(self.params['incident_polar_angle_deg']))
        
        fi = QtCore.QFileInfo(self.params.matFilename)
        self.ui.txtMaterial.setText(fi.baseName())
        
        #update global parameters for the dimer simulation
        self.ui.spinSpacing.setValue(d)
        
    def getParams(self):
        self.params.minLambda = self.ui.spinStartLambda.value()
        self.params.maxLambda = self.ui.spinEndLambda.value()
        self.params.nSamples = self.ui.spinNumSamples.value()
        self.params.nSpheres = self.ui.spinNumSpheres.value()
        self.params['incident_azimuth_angle_deg'] = self.ui.spinAlpha.value()
        self.params['incident_polar_angle_deg'] = self.ui.spinBeta.value()
        
      
        #global parameters for dimers
        d = self.ui.spinSpacing.value()
        
        return self.params
    
    def simulate(self):
        self.results = RunSimulation(self.params)
        
    def saveresults(self):
        fileName = QtGui.QFileDialog.getSaveFileName(w, 'Save Spectral Results', '', 'DAT data files (*.dat)')        
        if fileName:
            self.results.saveFile(fileName)
            
    def loadmaterial(self):
        fileName = QtGui.QFileDialog.getOpenFileName(w, 'Load Material Refractive Index', '', 'TXT data files (*.txt)')
        if fileName:
            self.params.matFilename = fileName
            
            fi = QtCore.QFileInfo(fileName)
            self.ui.txtMaterial.setText(fi.baseName())
        
    def __init__(self):
        QtGui.QWidget.__init__(self)        
        
        #dimer-specific settings
        self.params['number_spheres'] = 2
        self.params['sphere_position_file'] = ''
        
        #load the UI window
        self.ui = uic.loadUi('mstm_guiwindow.ui')
        #update the displayed parameters
        self.setParams()
        #display the UI
        self.ui.show()
        
        #simulation button
        self.connect(self.ui.btnSimulate, QtCore.SIGNAL("clicked()"), self.simulate)
        self.connect(self.ui.mnuSaveResults, QtCore.SIGNAL("triggered()"), self.saveresults)
        self.connect(self.ui.mnuLoadMaterial, QtCore.SIGNAL("triggered()"), self.loadmaterial)

class ProgressBar(QtGui.QWidget):
    def __init__(self, parent=None, total=20):
        super(ProgressBar, self).__init__(parent)
        self.name_line = QtGui.QLineEdit()

        self.progressbar = QtGui.QProgressBar()
        self.progressbar.setMinimum(1)
        self.progressbar.setMaximum(total)

        main_layout = QtGui.QGridLayout()
        main_layout.addWidget(self.progressbar, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")

    def update_progressbar(self, val):
        self.progressbar.setValue(val)
        

def RunSimulation(parameters):
    
    #load the material
    material = MaterialClass(parameters.matFilename)

    #add water if necessary
    if parameters.inWater:
        material.addSolution(1.33)
    
    #set the parameters based on the UI
    parameters = w.getParams()
    
    #range for simulation
    minLambda = parameters.minLambda
    maxLambda = parameters.maxLambda
    nSamples = parameters.nSamples

    #store the simulation results
    results = SimParserClass()
    
    #create a progress bar
    pbar = ProgressBar(total=nSamples)
    pbar.show()
    
    #for each wavelength in the material
    for i in range(nSamples):

        l = minLambda + i*(maxLambda - minLambda)/(nSamples - 1)
        
        

        #set the computed parameters
        m = material[l]
        n = m.n
        parameters['real_ref_index_scale_factor'] = n.real
        parameters['imag_ref_index_scale_factor'] = n.imag
        parameters['length_scale_factor'] = (2.0 * 3.14159)/l
        parameters['scattering_plane_angle_deg'] = gamma;
        

        parameters.clearSpheres()
        parameters.addSphere(a, -(d + 2*a)/2, 0, 0)
        parameters.addSphere(a, (d + 2*a)/2, 0, 0)

        #save the scripted input file
        parameters.saveFile('scriptParams.inp')

        #run the binary
        from subprocess import call
        devnull = open('/dev/null', 'w')
        call(["./ms-tmatrix",  "scriptParams.inp"], stdout=devnull)

        results.parseSimFile(l, 'test.dat')

        #update the progress bar
        pbar.update_progressbar(i+1)


    #plot results of interest
    import matplotlib.pyplot as plt
    wl = results['lambda']
    unpol = results['extinction_unpolarized']
    para = results['extinction_parallel']
    perp = results['extinction_perpendicular']
    plt.plot(wl, unpol, 'r-', wl, para, 'g-', wl, perp, 'b-')
    plt.ylabel('Extinction')
    plt.xlabel('Wavelength (um)')
    plt.show()
    
    #return the results
    return results;


        



#input template file name
inpFilename = 'msinput.inp'

#output spectral file name
outFilename = 'spectralOut.txt'

#sphere radii
a = 0.025
#distance between spheres
d = 0.002
#incident light directions
alpha = 0
beta = 0
gamma = 0

#results stored for each spectral sample
resultLabels = {'lambda', 'extinction_unpolarized', 'extinction_parallel', 'extinction_perpendicular'}





outFile = open(outFilename, 'w')

#number of characters in the progress bar
pb_max = 50



#create a Qt window
app = QtGui.QApplication(sys.argv)
w = GuiWindow()
sys.exit(app.exec_())





