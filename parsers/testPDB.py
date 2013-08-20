from utilPDB import *
from atomSelect import *

pdb = pdbObject("../../1A6Z.pdb")

s = select("//A/1-50,59,70-75///CA/",pdb.atom)

