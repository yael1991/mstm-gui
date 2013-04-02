class SimParserClass:

    simResults = dict()

    def __init__(self):
        self.simResults['lambda'] = list()
        self.simResults['extinction_unpolarized'] = list()
        self.simResults['extinction_parallel'] = list()
        self.simResults['extinction_perpendicular'] = list()
    
    def parseSimFile(self, l, fileName):
        self.simResults['lambda'].append(l)
        inFile = open(fileName, 'r')

        while True:
            line = inFile.readline().strip()

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

    def saveFile(self, fileName):
        outFile = open(fileName, 'w')
        outFile.write(str(self))

    def __getitem__(self, key):
        return self.simResults[key];

    def __str__(self):
        result = '';

        for i in range(len(self.simResults['lambda'])):
            result += str(self.simResults['lambda'][i]) + '\t'
            result += str(self.simResults['extinction_unpolarized'][i]) + '\t'
            result += str(self.simResults['extinction_parallel'][i]) + '\t'
            result += str(self.simResults['extinction_perpendicular'][i]) + '\n'

        return result
