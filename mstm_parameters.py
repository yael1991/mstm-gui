class ParameterClass:
    #minimum and maximum wavelengths for the simulation
    minLambda = 0.300
    maxLambda = 0.700
    #number of spectral samples
    nSamples = 40
    
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
        self.paramDict[key] = str(value);

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
            else:                
                value = inpFID.readline().strip()
                self.paramDict[key] = value

        inpFID.close()

    def saveFile(self, fileName):

        #open the output file
        outFID = open(fileName, 'w')

        #write the parameters
        for key in self.paramDict.keys():
            outFID.write(key + '\n')
            outFID.write(self.paramDict[key] + '\n')

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
            result += key + ": " + self.paramDict[key] + '\n'

        result += "\n"
        result += "Spheres:\n"
        #iterate through each sphere
        for s in self.sphereList:
            result += "------------------\n"
            for i in range(len(s)):
                result += self.sphereParamNames[i] + ": " + str(s[i]) + '\n'

        return result
