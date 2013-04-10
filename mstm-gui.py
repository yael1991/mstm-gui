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

#Matplotlib libraries
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pylab import *

class GuiWindow(QtGui.QMainWindow):
	
	params = ParameterClass('msinput.inp')
	
	def setParams(self):
		#update the Gui based on values in the parameters structure
		self.ui.spinStartLambda.setValue(self.params.minLambda)
		self.ui.spinEndLambda.setValue(self.params.maxLambda)
		self.ui.spinNearFieldLambda.setValue(self.params.snapshotLambda)
		self.ui.spinNumSamples.setValue(self.params.nSamples)
		self.ui.spinNumSpheres.setValue(int(self.params['number_spheres']))
		#near field stuff
		self.ui.cmbPlaneSlice.setCurrentIndex(int(self.params['near_field_plane_coord']) - 1)
		verts = self.params['near_field_plane_vertices']
		self.ui.spinNearFieldWidth.setValue(verts[2] - verts[0])
		self.ui.spinNearFieldHeight.setValue(verts[3] - verts[1])
		self.ui.spinNearFieldSteps.setValue(self.params.nSteps)
		
		fi = QtCore.QFileInfo(self.params.matFilename)
		self.ui.txtMaterial.setText(fi.baseName())
		
		#update global parameters for the dimer simulation
		self.ui.spinSpacing.setValue(self.params.d)
		self.ui.spinRadius.setValue(self.params.a)
		
	def getParams(self):
		self.params.minLambda = self.ui.spinStartLambda.value()
		self.params.maxLambda = self.ui.spinEndLambda.value()
		self.params.snapshotLambda = self.ui.spinNearFieldLambda.value()
		self.params.nSamples = self.ui.spinNumSamples.value()
		self.params['number_spheres'] = self.ui.spinNumSpheres.value()
		
		#incident light properties
		if self.ui.chkRandomOrientation.isChecked():
			self.params['fixed_or_random_orientation'] = 1
		else:
			self.params['fixed_or_random_orientation'] = 0
		self.params['incident_azimuth_angle_deg'] = self.ui.spinAlpha.value()
		self.params['incident_polar_angle_deg'] = self.ui.spinBeta.value()
		self.params['polarization_angle_deg'] = self.ui.spinGamma.value()
		
		self.params.showOutput = self.ui.chkShowOutput.isChecked()
		self.params.inWater = self.ui.chkInWater.isChecked()
		
			
		#near field
		if self.ui.chkNearField.isChecked():
			self.params['calculate_near_field'] = 1
		else:
			self.params['calculate_near_field'] = 0
		self.params['near_field_plane_coord'] = self.ui.cmbPlaneSlice.currentIndex() + 1
		width = (self.ui.spinNearFieldWidth.value()/2)
		height = (self.ui.spinNearFieldHeight.value()/2)
		self.params['near_field_plane_vertices'] = [-width, -height, width, height]
		dx = self.ui.spinNearFieldWidth.value() / (self.ui.spinNearFieldSteps.value() - 1)
		self.params['spacial_step_size'] = dx
	
		#global parameters for dimers
		self.params.d = self.ui.spinSpacing.value()
		self.params.a = self.ui.spinRadius.value()
		
		#get the spheres from the table
		nSpheres = self.ui.tblSpheres.rowCount()
		print("Row count: " + str(nSpheres))
		print("Orientatino: " + str(self.params['fixed_or_random_orientation']))
		
		self.params.sphereList = []
		for s in range(nSpheres):
			a = float(self.ui.tblSpheres.item(s, 0).text())
			x = float(self.ui.tblSpheres.item(s, 1).text())
			y = float(self.ui.tblSpheres.item(s, 2).text())
			z = float(self.ui.tblSpheres.item(s, 3).text())
			self.params.addSphere(a, x, y, z)
	
		return self.params
	
	def simulate(self):
		self.results = RunSimulation(True)
		
		#plot results of interest
		wl = self.results['lambda']
		
		if int(self.params['fixed_or_random_orientation']) == 0:
			unpol = self.results['extinction_unpolarized']
			para = self.results['extinction_parallel']
			perp = self.results['extinction_perpendicular']
			plt.plot(wl, unpol, 'r-', label='unpolarized')
			plt.plot(wl, para, 'g-', label='parallel')
			plt.plot(wl, perp, 'b-', label='perpendicular')
		else:
			total = self.results['extinction_total']
			plt.plot(wl, total, 'r-', label='extinction')
			
		#plot the near field maximum values if available
		
		if self.params['calculate_near_field']:
			maxima = self.results.maxNearField
			print(len(wl))
			print(len(maxima))
			plt.plot(wl, maxima)
		
		
		
		plt.legend(loc = 'upper left')
		plt.ylabel('Extinction')
		plt.xlabel('Wavelength (um)')
		plt.show()
		
	def func3(self, x,y):
			return (1- x/2 + x**5 + y**3)*exp(-x**2-y**2)
			
	def snapshot(self):
	
		self.results = RunSimulation(False)
		
		if self.params['calculate_near_field']:
			#verts = self.params['near_field_plane_vertices']
			#dx = (verts[2] - verts[0])/(self.params.nSteps)
			#x = arange(verts[0], verts[2], dx)
			#print(len(x))
			#y = arange(verts[1], verts[3], dx)
			#X, Y = meshgrid(x, y)
			E = array(self.results.gridNearField)
			#pcolor(X, Y, E, cmap=cm.RdBu)
			#colorbar()
			#axis([verts[0], verts[2], verts[1], verts[3]])
			
			pcolor(E, cmap=cm.RdBu)
			colorbar()
			print("Maximum enhancement: " + str(abs(E).max()))
		
		# make these smaller to increase the resolution
		#dx, dy = 0.05, 0.05

		#x = arange(-3.0, 3.0001, dx)
		#y = arange(-3.0, 3.0001, dy)
		#X,Y = meshgrid(x, y)

		#Z = self.func3(X, Y)
		#pcolor(X, Y, Z, cmap=cm.RdBu, vmax=abs(Z).max(), vmin=-abs(Z).max())
		#colorbar()
		#axis([-3,3,-3,3])

		show()
		
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
			
	def spherenum(self, i):
		self.ui.tblSpheres.setRowCount(i)
		print(i)
		
	def updatedimers(self):
		
		d = self.ui.spinSpacing.value()
		a = self.ui.spinRadius.value()
		
		self.ui.tblSpheres.setItem(0, 0, QtGui.QTableWidgetItem(str(a)))
		self.ui.tblSpheres.setItem(0, 1, QtGui.QTableWidgetItem(str(-(d + 2*a)/2)))
		self.ui.tblSpheres.setItem(0, 2, QtGui.QTableWidgetItem(str(0.0)))
		self.ui.tblSpheres.setItem(0, 3, QtGui.QTableWidgetItem(str(0.0)))
		
		self.ui.tblSpheres.setItem(1, 0, QtGui.QTableWidgetItem(str(a)))
		self.ui.tblSpheres.setItem(1, 1, QtGui.QTableWidgetItem(str((d + 2*a)/2)))
		self.ui.tblSpheres.setItem(1, 2, QtGui.QTableWidgetItem(str(0.0)))
		self.ui.tblSpheres.setItem(1, 3, QtGui.QTableWidgetItem(str(0.0)))
		
		
	def __init__(self):
		QtGui.QWidget.__init__(self)        
		
		#dimer-specific settings
		self.params['number_spheres'] = 2
		self.params['sphere_position_file'] = ''
				
		#load the UI window
		self.ui = uic.loadUi('mstm_guiwindow.ui')
		
		
		
		#controls
		self.connect(self.ui.btnSimulate, QtCore.SIGNAL("clicked()"), self.simulate)
		self.connect(self.ui.btnEvaluateNearField, QtCore.SIGNAL("clicked()"), self.snapshot)
		self.connect(self.ui.mnuSaveResults, QtCore.SIGNAL("triggered()"), self.saveresults)
		self.connect(self.ui.mnuLoadMaterial, QtCore.SIGNAL("triggered()"), self.loadmaterial)
		self.connect(self.ui.spinNumSpheres, QtCore.SIGNAL("valueChanged(int)"), self.spherenum)
		self.connect(self.ui.spinRadius, QtCore.SIGNAL("valueChanged(double)"), self.updatedimers)
		self.connect(self.ui.spinSpacing, QtCore.SIGNAL("valueChanged(double)"), self.updatedimers)
		
		#update the displayed parameters
		self.setParams()
		
		#update the sphere table with the default dimer values
		self.updatedimers()
		
		#display the UI
		self.ui.show()


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
		

def RunSimulation(spectralSim = True):
	
	#set the parameters based on the UI
	parameters = w.getParams()
	
	
	
	#load the material
	material = MaterialClass(parameters.matFilename)

	#add water if necessary
	if parameters.inWater:
		material.addSolution(1.33)

	#for a spectral simulation, set the range and number of samples
	if spectralSim:
		minLambda = parameters.minLambda
		maxLambda = parameters.maxLambda
		nSamples = parameters.nSamples
	else:
		minLambda = parameters.snapshotLambda
		maxLambda = parameters.snapshotLambda
		nSamples = 1

	#store the simulation results
	results = SimParserClass(parameters)
	
	#create a progress bar
	pbar = ProgressBar(total=nSamples)
	pbar.show()
	
	#for each wavelength in the material
	for i in range(nSamples):

		if i == 0:
			l = minLambda
		else:
			l = minLambda + i*(maxLambda - minLambda)/(nSamples - 1)

		#set the computed parameters
		m = material[l]
		n = m.n
		parameters['real_ref_index_scale_factor'] = n.real
		parameters['imag_ref_index_scale_factor'] = n.imag
		parameters['length_scale_factor'] = (2.0 * 3.14159)/l
		parameters['scattering_plane_angle_deg'] = gamma;
		parameters['near_field_output_data'] = 0
		#parameters['number_spheres'] = 1

		#a = parameters.a;
		#d = parameters.d;
		#parameters.clearSpheres()
		#parameters.addSphere(a, -(d + 2*a)/2, 0, 0)
		#parameters.addSphere(a, (d + 2*a)/2, 0, 0)

		#save the scripted input file
		parameters.saveFile(l, 'scriptParams.inp')

		#run the binary
		from subprocess import call
		if parameters.showOutput:
			call(["./ms-tmatrix",  "scriptParams.inp"])
		else:            
			devnull = open('/dev/null', 'w')
			call(["./ms-tmatrix",  "scriptParams.inp"], stdout=devnull)

		#parse the simulation results
		results.parseSimFile(l, 'test.dat')
		
		if parameters['calculate_near_field']:
			results.parseNearField('nf-temp.dat')
		

		#update the progress bar
		pbar.update_progressbar(i+1)
	
	#return the results
	return results;


		




#incident light directions
alpha = 0
beta = 0
gamma = 0

#results stored for each spectral sample
resultLabels = {'lambda', 'extinction_unpolarized', 'extinction_parallel', 'extinction_perpendicular'}

#create a Qt window
app = QtGui.QApplication(sys.argv)
w = GuiWindow()
sys.exit(app.exec_())





