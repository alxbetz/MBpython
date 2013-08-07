
fNames=["record name","atom serial number","atom name","alternate location indicatior","residue name","chain identifier","residue sequence number","code for insertions of residue","x coord in A","y coord in A","z-coord in A","occupancy","temperature factor","segment identifier","element symbol","charge"]
fIdent=["recName","serialNum","atomName","alternateLoc","resName","chainID","resSeqNum","insCode","X","Y","Z","occupancy","temp","segID","element","charge"]
fType=["string","int","string","character","string","character","integer","character","float","float","float","float","float","string","string","string"]

fSizes=[6,6,4,1,4,1,4,1,3,8,8,8,6,6,6,4,2,2]

cutPos = []
a = 0
for i in fSizes:
	# skip gaps
	if i == 9 or i == 14:
		a+=i
		continue
	cutPos.append([a,a+i])
	a+=i
if len(cutPos) != len(fIdent):
	print "smt wrong here"

def parseAtom(entry):
	res = []
	if len(entry != 80):
		print "length is " + len(entry)
	for cut,i in cutPos:
		dat = entry[cut[0],cut[1]].strip()
		if fType[i]=="int":
			dat = int(dat)
		elif fType[i]=="float":
			dat = float(dat)
		res.append(dat)
	return res
	
def parseConnect(entry):
	dat = line.split()
	return [dat[1],dat[2:len(dat)-1]]
	
