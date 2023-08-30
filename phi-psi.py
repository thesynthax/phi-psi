#program to calculate phi and psi dihedral angles

pdbFile = './1asy.pdb';

class atom:
    def __init__(self, line):
        self.atomNumber = int(line[6:11].strip());
        self.atomName = line[13:16].strip();
        self.residueName = line[17:20].strip();
        self.chainIdentifier = line[21];
        self.residueNumber = int(line[22:28].strip());
        self.xPos = float(line[30:38].strip());
        self.yPos = float(line[38:46].strip());
        self.zPos = float(line[46:54].strip());

class residue:
    def __init__(self, atomObjects):
        pass

fptr = open(pdbFile, 'r');
lines = fptr.readlines();
fptr.close();

atomLines = [];
for line in lines:
    if line.startswith('ATOM'):
        atomLines.append(line);

atomObjects = [];
for atomLine in atomLines:
    atomObjects.append(atom(atomLine));

proteinObjects = [];
rnaObjects = [];
for atomObject in atomObjects:
    if len(atomObject.residueName) == 3:
        proteinObjects.append(atomObject);
    else:
        rnaObjects.append(atomObject);

residues = [];
currentResidue = []
for proteinObjectIndex in range(len(proteinObjects)):
    currentResidue.append(proteinObjects);
    if proteinObjects[proteinObjectIndex].residueNumber != proteinObjects[proteinObjectIndex + 1].residueNumber:
        residues.append(residue(currentResidue));
        currentResidue = []
