class MaterialSampleClass:
    
    #constructor      
    def __init__(self, l, n):
        self.l = l
        self.n = n

    #string conversion
    def __str__(self):
        result = ""
        result += str(self.l) + 'um: ' + str(self.n)
        return result
        
class MaterialClass:
    materialList = []

    def __init__(self, fileName=""):
        
        self.materialList = []
        
        if fileName != "":
            self.loadFile(fileName)

    #when the material is cast to a string, create the list of refractive indices
    def __str__(self):
        nSamples = len(self.materialList)
        result = ""
        for i in range(nSamples):
            result += str(self.materialList[i]) + '\n'
        return result

    def __len__(self):
        return len(self.materialList)

    def __getitem__(self, l):
        bigI = smallI = 0;
        bigV = 999;
        smallV = 0;
        #find the smallest sample larger than l
        for i in range(len(self.materialList)):
            if self.materialList[i].l > l and self.materialList[i].l < bigV:
                bigI = i
                bigV = self.materialList[i].l
            if self.materialList[i].l < l and self.materialList[i].l > smallV:
                smallI = i
                smallV = self.materialList[i].l

        a = (l - smallV)/(bigV - smallV)

        bigN = self.materialList[bigI].n
        smallN = self.materialList[smallI].n

        n = a * bigN + (1 - a) * smallN

        return MaterialSampleClass(l, n)
        

        #print(str(self.materialList[smallI].l) + "---" + str(self.materialList[bigI].l))
        
        return self.materialList[smallI]

    def add(self, l, n):
        m = MaterialSampleClass(l, n)
        self.materialList.append(m)

    def clip(self, minLambda, maxLambda):
        #this function clips all material samples to the range [minLambda, maxLambda]
        self.materialList = list(filter(lambda m: m.l > minLambda, self.materialList))
        self.materialList = list(filter(lambda m: m.l < maxLambda, self.materialList))


    def addSolution(self, n):
        #places the material in a solution (divide by the solution's n)
        for i in range(len(self.materialList)):
            self.materialList[i].n = self.materialList[i].n / n

        

    def loadFile(self, fileName):
        #open the real refractive index file
        irFID = open(fileName, 'r')
        #read the first line to get the units (wavelength (um) or wavenumber (cm^2))
        lightUnits = irFID.readline().split()[0]

        #load the material
        for line in irFID:
            l, n, k = map(float, line.split())

            #if units are in wavenumber, convert to wavelength
            if lightUnits == "nu":
                l = l/10000

            self.add(l, complex(n, k))

        #close the file
        irFID.close()
