#this code parses the results of a simulation and stores them in a SimParserClass structure
from pylab import *

class SimParserClass:

    simResults = dict()
    
    sxNearField = 0
    syNearField = 0
    intersectedNearField = 0
    
    #near field data
    gridNearField = []
    maxNearField = []
    
    #the stokes matrix is read from the output
    stokesMatrix = []
    scatAmpMatrix = []
    

    def __init__(self, parameters):
        
        self.params = parameters;
        
        self.simResults['lambda'] = list()      

        self.simResults['extinction_unpolarized'] = list()
        self.simResults['extinction_parallel'] = list()
        self.simResults['extinction_perpendicular'] = list()
        self.simResults['extinction_total'] = list()
        self.simResults['detector_field'] = list()
        
        self.gridNearField = []
        self.maxNearField = []
        
        
    def parseSimFile(self, l, fileName):
        self.simResults['lambda'].append(l)
        inFile = open(fileName, 'r')

        
        while True:
            
            line = inFile.readline().strip()
			
			#if the simulation is for a single plane wave
            if int(self.params['fixed_or_random_orientation']) == 0:
                if line == 'scattering matrix elements':
					#empty the stokes matrix
					self.stokesMatrix = []
					inFile.readline()
					for s in range(0, 181):
						values = map(float, inFile.readline().strip().split())
						self.stokesMatrix.append(values)
					break;
                elif line == 'unpolarized total ext, abs, scat efficiencies, w.r.t. xv, and asym. parm':
                    values = inFile.readline().strip().split(' ')
                    self.simResults['extinction_unpolarized'].append(values[0])
                elif line == 'parallel total ext, abs, scat efficiencies':
                    values = inFile.readline().strip().split(' ')
                    self.simResults['extinction_parallel'].append(values[0])
                elif line == 'perpendicular total ext, abs, scat efficiencies':
                    values = inFile.readline().strip().split(' ')
                    self.simResults['extinction_perpendicular'].append(values[0])
                
            #if the simulation is for random orientations
            else:
                if line == 'scattering matrix elements':
                    break
                elif line == 'total ext, abs, scat efficiencies, w.r.t. xv, and asym. parm':
                    values = inFile.readline().strip().split(' ')
                    self.simResults['extinction_total'].append(values[0])
                    
    def parseNearField(self, fileName):
        
        inFile = open(fileName, 'r')
        
        #get the size of the near field grid
        line = inFile.readline().strip()
        self.sxNearField, self.syNearField = map(int, line.split())
        
        #get the number of spheres that are intersected
        line = inFile.readline().strip()
        self.intersectedNearField = int(line)
        
        #process intersections here-----------
        

        #get the field values
        self.gridNearField = []
        for y in range(self.syNearField):
            self.gridNearField.append([])
            for x in range(self.sxNearField):
                line = inFile.readline().strip()
                values = map(float, line.split())
                self.gridNearField[y].append(values[2])
                
        E = array(self.gridNearField)
        self.maxNearField.append(abs(E).max())
        
    #calculate and return the scattering amplitude matrix
    def calcScatteringAmp(self):
		#compute the number of entries in the stokes matrix
		nEntries = len(self.stokesMatrix)
		
		#initialize the scattering amplitude matrix to empty
		self.scatAmpMatrix = []
		
		for s in range(0, nEntries):
			Z = self.stokesMatrix[s]
			
			scatEntry = []
			s11 = complex(sqrt(0.5 * (Z[1] - Z[2] - Z[5] + Z[6])), 0.0)
			scatEntry.append(s11)
			scatEntry.append(complex(-0.5 * (Z[3] + Z[7]) / s11, 0.5 * (Z[4] + Z[8]) / s11))
			scatEntry.append(complex(-0.5 * (Z[9] + Z[10]) / s11, -0.5 * (Z[13] + Z[14]) / s11))
			scatEntry.append(complex(0.5 * (Z[11] + Z[12]) / s11, -0.5 * (Z[12] - Z[15]) / s11))
			
			self.scatAmpMatrix.append(scatEntry)
			
		S = self.scatAmpMatrix[0]
		E = [S[0], S[2]]
		self.simResults['detector_field'].append(E)
		print(E)
               

    def saveFile(self, fileName):
        outFile = open(fileName, 'w')
        outFile.write(str(self))
        outFile.close()

    def __getitem__(self, key):
        return self.simResults[key];

    def __str__(self):
        result = '';

        for i in range(len(self.simResults['lambda'])):
            result += str(self.simResults['lambda'][i])
            result += '\t' + str(self.simResults['extinction_unpolarized'][i])
            result += '\t' + str(self.simResults['extinction_parallel'][i])
            result += '\t' + str(self.simResults['extinction_perpendicular'][i])
            result += '\t' + str(self.simResults['detector_field'][i][0]) + '\t' + str(self.simResults['detector_field'][i][1])
            
            #parse the near field if it is included in the simulation
            #if int(parameters['calculate_near_field']) == 1:
            #    result += '\t' + str(maxNearField)
                
            result += '\n'
        return result
