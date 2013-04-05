class ParameterClass:
    #minimum and maximum wavelengths for the simulation
    minLambda = 0.300
    maxLambda = 0.700
    snapshotLambda = 0.300
    #number of spectral samples
    nSamples = 40
    
    #spatial samples
    nSteps = 100
    
    #sphere size and separation
    a = 0.025
    d = 0.005
    
    #material file name
    matFilename = 'etaSilver.txt'
    #are the sphere's in water?
    inWater = False
    
    #show console output from the MS-TM FORTRAN program
    showOutput = False

    paramDict = {}
    sphereList = []

    sphereParamNames = ['radius', 'X', 'Y', 'Z', 'n', 'k', 'Xr', 'Xi']

    def __init__(self, fileName):
        self.loadFile(fileName)

    def __getitem__(self, key):
        return self.paramDict[key];

    def __setitem__(self, key, value):
        self.paramDict[key] = value;

    def clearSpheres(self):
        self.sphereList = []

    def addSphere(self, a, x, y, z, n = 1.0, k=1.0):
        s = [a, x, y, z, n, k]
        self.sphereList.append(s)

    def loadFile(self, fileName):
        inpFID = open(fileName, 'r')
        selfparamDict = {}

        while 1:
            key = inpFID.readline().strip()
            
            
            #deal with sphere sizes and positions
            if key == 'sphere_sizes_and_positions':
                
                while True:
                    #load the parameters for a sphere
                    value = inpFID.readline().strip()
                    if value == 'end_of_options':
                        break
                    
                    self.sphereList.append(value.split(' '))           
                    

            elif not key:
                break
            elif key == 'end_of_options':
                break
            #deal with the near field plane
            elif key == 'near_field_plane_vertices':
                value = inpFID.readline().strip()
                #self.paramDict[key] = list(map(float, value.split(',')))
                self.paramDict[key] = [-1, -1, 1, 1]
            
            else:                
                value = inpFID.readline().strip()
                self.paramDict[key] = value
                
        #update the length scale factor to deal with the UI
        self.paramDict['length_scale_factor'] = (2.0 * 3.14159)/self.snapshotLambda

        inpFID.close()

    def saveFile(self, l, fileName):
        
        #print(self)

        #open the output file
        outFID = open(fileName, 'w')

        #write the parameters
        for key in self.paramDict.keys():
            outFID.write(key + '\n')
            
            #deal with the near field plane
            if key == 'near_field_plane_vertices':
                #these have to be scaled by the length scale factor
                ls = (2 * 3.14159)/l
                v = self.paramDict[key]
                outFID.write(str(v[0]*ls) + ',' + str(v[1]*ls) + ',' + str(v[2]*ls) + ',' + str(v[3]*ls) + '\n')
            elif key == 'spacial_step_size':
                ls = (2 * 3.14159)/l
                dx = self.paramDict[key] * ls
                outFID.write(str(dx) + '\n')
                
            else:
                outFID.write(str(self.paramDict[key]) + '\n')

        #write the spheres
        outFID.write("sphere_sizes_and_positions\n")
        for s in self.sphereList:
            for p in s:
                outFID.write(str(p) + ' ')
            outFID.write('\n')
        outFID.write("end_of_options")
        

    def __str__(self):
        #print(self.paramDict)
        result = ""
        for key in self.paramDict.keys():
            #deal with the near field plane
            if key == 'near_field_plane_vertices':
                v = map(str, self.paramDict[key])                
                result += key + ": " + v[0] + ',' + v[1] + ',' + v[2] + ',' + v[3] + '\n'
            else:
                result += key + ": " + str(self.paramDict[key]) + '\n'

        result += "\n"
        result += "Spheres:\n"
        #iterate through each sphere
        for s in self.sphereList:
            result += "------------------\n"
            for i in range(len(s)):
                result += self.sphereParamNames[i] + ": " + str(s[i]) + '\n'

        return result
