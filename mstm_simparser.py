from pylab import *

class SimParserClass:

    simResults = dict()
    
    sxNearField = 0
    syNearField = 0
    intersectedNearField = 0
    
    #near field data
    gridNearField = []
    maxNearField = []
    

    def __init__(self, parameters):
        
        self.params = parameters;
        
        self.simResults['lambda'] = list()      

        self.simResults['extinction_unpolarized'] = list()
        self.simResults['extinction_parallel'] = list()
        self.simResults['extinction_perpendicular'] = list()
        self.simResults['extinction_total'] = list()
        
        self.gridNearField = []
        self.maxNearField = []
        
        
    def parseSimFile(self, l, fileName):
        self.simResults['lambda'].append(l)
        inFile = open(fileName, 'r')

        
        while True:
            
            line = inFile.readline().strip()

            if int(self.params['fixed_or_random_orientation']) == 0:
                if line == 'scattering matrix elements':
                    break
                elif line == 'unpolarized total ext, abs, scat efficiencies, w.r.t. xv, and asym. parm':
                    values = inFile.readline().strip().split(' ')
                    self.simResults['extinction_unpolarized'].append(values[0])
                elif line == 'parallel total ext, abs, scat efficiencies':
                    values = inFile.readline().strip().split(' ')
                    self.simResults['extinction_parallel'].append(values[0])
                elif line == 'perpendicular total ext, abs, scat efficiencies':
                    values = inFile.readline().strip().split(' ')
                    self.simResults['extinction_perpendicular'].append(values[0])
            
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
        self.maxNearField.append(pow(E.max(), 2))
               

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
            
            #parse the near field if it is included in the simulation
            if int(parameters['calculate_near_field']) == 1:
                result += '\t' + str(maxNearField)
                
            result += '\n'
        return result
