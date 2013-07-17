import Bio
import molGraph
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from optparse import OptionParser

pdb_parser=PDBParser()



parser = OptionParser()
parser.add_option("-f", "--from", dest="fromFile",
                  help="take coordinates from this file", metavar="FROM")
parser.add_option("-t", "--to", dest="toFile",
                  help="inject coordinates into this file", metavar="TO")
parser.add_option("-o", "--output", dest="out",
                  help="output file with injected coordinates", metavar="OUT")


(options, args) = parser.parse_args()
fromFile = options.fromFile
toFile = options.toFile
out = options.out

structFrom = pdb_parser.get_structure("rep1", fromFile)
mFrom = structFrom[0]
atomsFrom = list(mFrom.get_atoms())

structTo = pdb_parser.get_structure("rep1", toFile)
mTo = structTo[0]
atomsTo = list(mTo.get_atoms())

n1 ,e1 = molGraph.makeGraph(mFrom)
n2 ,e2 = molGraph.makeGraph(mTo)

iso = molGraph.getIsomorphism(e1,e2)

for i in range(len(atomsFrom)):
	f = iso[0][i]
	t = iso[1][i]
	atomsTo[t].set_coord(atomsFrom[f].get_coord())
	
w = PDBIO()
w.set_structure(structTo)
w.save(out)
