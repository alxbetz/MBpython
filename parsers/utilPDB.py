from collections import defaultdict

fNames=["record name","atom serial number","atom name","alternate location indicatior","residue name","chain identifier","residue sequence number","code for insertions of residue","x coord in A","y coord in A","z-coord in A","occupancy","temperature factor","segment identifier","element symbol","charge"]
fIdent=["recName","serialNum","atomName","alternateLoc","resName","chain","resNum","insCode","X","Y","Z","occupancy","temp","segid","element","charge"]
fType=["string","int","string","character","string","character","int","character","float","float","float","float","float","string","string","string"]
fSizes=[6,6,4,1,4,1,4,1,3,8,8,8,6,6,6,4,2,2]


#~ cutPos = []
#~ a = 0
#~ for i,val in enumerate(fSizes):
	#~ # skip gaps
	#~ if i == 8 or i == 14:
		#~ a+=val
		#~ continue
	#~ else:
		#~ cutPos.append([a,a+val])
		#~ a+=val
		
cutPos = [[0, 6], [6, 12], [12, 16], [16, 17], [17, 21], [21, 22], [22, 26], [26, 27], [30, 38], [38, 46], [46, 54], [54, 60], [60, 66], [72, 76], [76, 78], [78, 80]]


def parseAtomDict(entry):
	global cutPos
	global fIdent
	res = defaultdict(lambda: "N/A")
	if len(entry) != 80:
		print "length is " + str(len(entry))
	for i,cut in enumerate(cutPos):
		dat = entry[cut[0],cut[1]].strip()
		if len(dat) == 0:
			continue
		elif fType[i]=="int":
			dat = int(dat)
		elif fType[i]=="float":
			dat = float(dat)
		res[fIdent[i]] = dat
	return res
	
def parseAtomArr(entry):
	global cutPos
	res = []
	#~ if len(entry) != 80:
		#~ print "length is " + str(len(entry))
	for i,cut in enumerate(cutPos):
		
		dat = entry[cut[0]:cut[1]].strip()
		if len(dat) == 0:
			res.append(None)
			continue
		elif fType[i]=="int":
			dat = int(dat)
		elif fType[i]=="float":
			dat = float(dat)
		res.append(dat)
	return res
	
def parseConnect(entry):
	dat = line.split()
	return [dat[1],dat[2:len(dat)-1]]
	
	


class pdbObject(object):
	
	models = []
	mxyz = []
	mCount = 0 #modelCount
	xyz = []
	atom = []
	'''
	Parse a PDB file
	If it has multiple models, the atomic records (ATOM and HETATM) are stored in 
	self.models, otherwise the atomic records are stored in self.atom 
	'''
	@classmethod
	def parsePDB(self,fName):
		#resList = []
		try:
			multiModel = False
			with open(fName) as f:
				for line in f:
					# check if file contains multiple models
					if line.startswith("Model"):
						mCount += 1
						#mNum = int(line.split()[1])
						if mCount > 1:
							models.append(atom)
							mxyz.append(xyz)
							atom = []
							xyz = []
							multiModel = True
					if line.startswith("ATOM") or line.startswith("HETATM"):
						record = parseAtomArr(line)
						self.atom.append(record)
						self.xyz.append(record[8:11])
				if multiModel:
					#handle last model
					models.append(atom)
					mxyz.append(xyz)
					atom = []
					xyz = []
					
		except IOError:
			print "File does not exist"
			system.exit()
	
	
	def __init__(self,fName):
		self.res = self.parsePDB(fName)



	
