#order by residue and reassign ids

import sys

def prep_field(size,value):
	value = str(value)
	l = len(value)
	if l == size:
		return value
	elif l > size:
		print "field larger than exptected"
		print value + " should only be " + size + "characters"
		return None
	else:
		pre = ''.join([' ' for i in range(size-l)])
		return pre + value

def prep_atomid(value):
	return prep_field(5,value)

def string_replace(fr,to,string,replace):
	return string[:fr] + replace + string[to:]

def parseFile(fName):
	fHandle = open(fName)
	out = []
	for line in fHandle:
		if line.startswith("ATOM") or line.startswith("HETATM"):
			out.append(line)
	return out

def compResID(string1,string2):
	return int(string1[22:26]) - int(string2[22:26])

def sortByResidue(atoms):
	return sorted(atoms, cmp=compResID)

fname = sys.argv[1]

data = parseFile(fname)
at_sort = sortByResidue(data)

out = open(fname,"w")

for i, atom in enumerate(at_sort):
	gt = string_replace(6,11,atom,prep_atomid(i+1))
	out.write(gt)


