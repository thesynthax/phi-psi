#program to calculate phi and psi dihedral angles
import math;

pdbFile = './1asy.pdb';

class vector:
    def __init__(self, x, y, z) -> None:
        self.x = x;
        self.y = y; 
        self.z = z;
        self.length = math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2);

class atom:
    def __init__(self, line):
        self.atomNumber = int(line[6:11].strip());
        self.atomName = line[13:16].strip();
        self.residueName = line[17:20].strip();
        self.chainIdentifier = line[21];
        self.residueNumber = int(line[22:28].strip());
        self.x = float(line[30:38].strip());
        self.y = float(line[38:46].strip());
        self.z = float(line[46:54].strip());

class residue:
    def __init__(self, atomObjects):
        self.residueName = atomObjects[0].residueName;
        self.residueNumber = atomObjects[0].residueNumber;
        self.chainIdentifier = atomObjects[0].chainIdentifier;
        self.phi = None;
        self.psi = None;
        self.n = None;
        self.c = None;
        self.ca = None;

        for atomObject in atomObjects:
            if atomObject.atomName == 'N':
                self.n = atomObject;
            elif atomObject.atomName == 'CA':
                self.ca = atomObject;
            elif atomObject.atomName == 'C':
                self.c = atomObject;

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
currentResidue = [];
for proteinObjectIndex in range(len(proteinObjects)):
    currentResidue.append(proteinObjects[proteinObjectIndex]);
    if proteinObjectIndex == len(proteinObjects) - 1 or proteinObjects[proteinObjectIndex].residueNumber != proteinObjects[proteinObjectIndex + 1].residueNumber:
        residues.append(residue(currentResidue));
        currentResidue = [];

def dot(v1, v2):
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

def cross(v1, v2):
    return vector(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);


def dihedral(atom1, atom2, atom3, atom4):
    #angle = "nigga";

    vecAB = vector(atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z);
    vecBC = vector(atom3.x - atom2.x, atom3.y - atom2.y, atom3.z - atom2.z);
    vecCD = vector(atom4.x - atom3.x, atom4.y - atom3.y, atom4.z - atom3.z);

    n1 = cross(vecAB, vecBC);
    n2 = cross(vecBC, vecCD);

    num1 = dot(n1, n2);
    denom1 = n1.length * n2.length;
    angle = math.acos(num1/denom1) * 180 / math.pi;

    n3 = cross(n1, n2);

    num2 = dot(n3, vecBC);
    denom2 = n3.length * vecBC.length;
    costheta = num2/denom2;

    angle = -angle if costheta < 0 else angle;

    return round(angle, 3);

for residueIndex in range(len(residues)):
    if residueIndex == 0 or residues[residueIndex-1].chainIdentifier != residues[residueIndex].chainIdentifier:
        residues[residueIndex].phi = 180;
        residues[residueIndex].psi = dihedral(residues[residueIndex].n, residues[residueIndex].ca, residues[residueIndex].c, residues[residueIndex+1].n);
    elif residueIndex == len(residues)-1 or residues[residueIndex+1].chainIdentifier != residues[residueIndex].chainIdentifier:
        residues[residueIndex].phi = dihedral(residues[residueIndex-1].c, residues[residueIndex].n, residues[residueIndex].ca, residues[residueIndex].c);
        residues[residueIndex].psi = 180;
    else:
        residues[residueIndex].phi = dihedral(residues[residueIndex-1].c, residues[residueIndex].n, residues[residueIndex].ca, residues[residueIndex].c);
        residues[residueIndex].psi = dihedral(residues[residueIndex].n, residues[residueIndex].ca, residues[residueIndex].c, residues[residueIndex+1].n);

phis = []
psis = []

for residue in residues:
    phis.append(residue.phi);
    psis.append(residue.psi);

print('Number\tPhi\tPsi')

for i in range(len(phis)):
    print(str(i) + '\t' + str(phis[i]) + '\t' + str(psis[i]));

from matplotlib import pyplot;

pyplot.title('Ramachandran plot', fontsize=15);
pyplot.xlabel('Phi', fontsize=13);
pyplot.ylabel('Psi', fontsize=13);
pyplot.axis([-180, 180, -180, 180]);
pyplot.axhline(0, color='black', lw=1);
pyplot.axvline(0, color='black', lw=1);
pyplot.grid(True);
pyplot.plot(phis, psis, 'co');
pyplot.savefig('1asy_rmc.png');
pyplot.show();
