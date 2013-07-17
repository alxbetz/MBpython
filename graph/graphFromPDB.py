import Bio
import numpy
import copy
from numpy import matrix
#from operator import itemgetter
from parsers import bondLength


from Bio.PDB.PDBParser import PDBParser

parser=PDBParser()

structure=parser.get_structure("test", "../peptides/antibiotic_peptide.pdb")
model=structure[0]

struct1 = parser.get_structure("rep1", "./test_data/cluster.rep.1.pdb")
m1 = struct1[0]
struct2 = parser.get_structure("rep2", "./test_data/cluster.rep.2.pdb")
m2 = struct2[0]

bLength = bondLength.getBondsMap()



		

def makeGraph(mol,method):
	
	
	gHash = {}
	atoms = list(mol.get_atoms())
	nodes = atoms
	
	size = len(atoms)
	dist = getDist(mol,1)
	edges = [[] for i in range(size)]
	for i,atom1 in enumerate(atoms):
		if atom1.element == "":
				assignElement(atom1)
		
		for j,atom2 in enumerate(atoms[:i]):
			if atom2.element == "":
				assignElement(atom2)
			bond = atom1.element + "-" + atom2.element
			#print atom1.name
			#print atom2.name
			#print bond
			
			
			if dist[i][j] < bLength[bond]:
				 edges[i].append(j)
				 edges[j].append(i)
				
	
	return nodes,edges

			
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def assignElement(atom):
	if not is_number(atom.name[0]):
		atom.element = atom.name[0]
	else:
		atom.element = atom.name[1] 				


def getDist(mol,full):
	atoms = list(mol.get_atoms())
	size = len(atoms)
	dist = numpy.zeros(shape=(size,size))
	for i,at1 in enumerate(atoms):
		for j,at2 in enumerate(atoms):
			if i == j:
				break
			else:
				dist[i][j] = numpy.linalg.norm(at1.get_coord()-at2.get_coord())
				if full:
					dist[j][i] = dist[i][j]
	return dist
