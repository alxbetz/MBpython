import Bio
import numpy
import math
import copy
from numpy import matrix
#from operator import itemgetter
from parsers import bondLength
import recogniseBonds

from Bio.PDB.PDBParser import PDBParser

parser=PDBParser()

#structure=parser.get_structure("test", "../peptides/antibiotic_peptide.pdb")
#model=structure[0]

#struct1 = parser.get_structure("rep1", "./test_data/shake_360_c0.pdb")
#m1 = struct1[0]
#struct2 = parser.get_structure("rep2", "./test_data/shake_360_c1.pdb")
#m2 = struct2[0]

bLength = bondLength.getBondsMap()



		

def makeGraph(mol):
	
	
	gHash = {}
	atoms = list(mol.get_atoms())
	nodes = atoms
	
	size = len(atoms)
	dist = getDist(mol,1)
	edges = [[] for i in range(size)]
	
	for atom in atoms:
		if atom.element == "":
			assignElement(atom)
			#print dir(atom)
			
	"""		
	for i,atom1 in enumerate(atoms):
		
		for j,atom2 in enumerate(atoms[:i]):
			
			bond = atom1.element + "-" + atom2.element
			#print atom1.name
			#print atom2.name
			#print bond
			
			
			if dist[i][j] < bLength[bond]:
				 edges[i].append(j)
				 edges[j].append(i)
				
	"""

	edges = recogniseBonds.rasMol(dist,atoms)
	#edges = getZMatrix(mol,dist)
	print "dist, edges"
	print len(dist)
	print len(edges)
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



def normalise(vec):
	out = []
	l = numpy.linalg.norm(vec)
	for comp in vec:
		out.append(comp/l)
	return out
	
	
def getNormalisedVectors(atoms,curr,candidates):
	out = []
	curr_xyz = atoms[curr].get_coord()
	for ca in candidates:
		#print atoms[curr].get_full_id()
		#print atoms[int(ca)].get_full_id()
		dir_vec = atoms[int(ca)].get_coord() - curr_xyz
		#print dir_vec
		v = normalise(dir_vec)
		#print v
		#print "---"
		out.append(v)
	return out
		
def getAngle(v1,v2):

	dot = numpy.dot(v1,v2)

	return (math.acos(dot) / (2 * math.pi)) * 360
	
def getMinIndex(array):
	minElem = array[0]
	minInd = 0
	for i,elem in enumerate(array):
		if elem < minElem:
			minElem = elem
			minInd = i
			
	return minInd

#put all possible permutations of a binary string of length rem in the list "total"
def getBinPerm(curr,rem,total):
	if rem==0:
		total.append(curr)
		return
		
	for i in range(2):
		new = list(curr)
		new.append(i)
		getBinPerm(new,rem-1,total)

#calculate the power set of inSet
def getPowerSet(inSet):
	ps =  []
	bin = []
	getBinPerm([],len(inSet),bin)
	for perm in bin:
		tmp = []
		for i,b in enumerate(perm):
			if b ==1:
				tmp.append(inSet[i])
		ps.append(tmp)
	return ps
	

def getMeanError(refAngle,subSet,angles):
	error = 0
	for item in subSet:
		error += math.fabs(refAngle-angles[item]) 
	error /= len(subSet)
	return error

	
def selPartners(nVecs,cInd):
	refAngles = [109.47,120,180]
	angles = []
	pairs = []
	for i,v1 in enumerate(nVecs):
		for j,v2 in enumerate(nVecs[0:i]):
			angles.append(getAngle(v1,v2))
			pairs.append([cInd[i],cInd[j]])
	
	subSets = getPowerSet(range(len(nVecs)))
	subSets.remove([])
	mSize = 0
	fAngle = 0
	fError = 0
	partners = []
	for subset in subSets:
		for refA in refAngles:
			
			error = getMeanError(refA,subset,angles)
			if error < 10 and len(subset) > mSize:
				tmp = []
				for elem in subset:
					tmp.append(pairs[elem])
				partners = tmp
				fAngle = refA
				fError = error
				
	#print partners
	#print fAngle
	#print fError
	uniqueP = set()
	for p in partners:
		
		uniqueP.add(p[0])
		uniqueP.add(p[1])
	
		
	return list(uniqueP)
				
			
def allBetween(array,upper,lower):
	outInd = []
	outElem = []
	for i,elem in enumerate(array):
		if elem < upper and elem > lower:
			outInd.append(i)
			outElem.append(elem)
	return outInd,outElem	
	
	
def getZMatrix(mol,dist):
	bonds = [[] for i in range(len(dist))]
	atoms = list(mol.get_atoms())
	geoms = [109.47,120,180]
	
	for i,row in enumerate(dist):
		cInd,cElem = allBetween(row,1.9,0)
		canLen = len(cInd)
		#print cInd
		if canLen == 1:
			bonds[i].append(cInd[0])
		elif canLen == 2:
			vecs = getNormalisedVectors(atoms,i,cInd)
			
			angle = getAngle(vecs[0],vecs[1])
			#print angle
			if angle >= 90:
				#bonds[i].__add__(cInd)
				bonds[i] = cInd
			else:
				minInd = getMinIndex(cElem)
				bonds[i].append(cInd[minInd])
		else:
			vecs = getNormalisedVectors(atoms,i,cInd)
			partners = selPartners(vecs,cInd)
			#print partners
			#bonds[i].__add__(partners)
			if len(partners) == 0:
				minInd = getMinIndex(cElem)
				bonds[i].append(cInd[minInd])
			else:
				bonds[i] = partners
	return bonds
				
				

def getMinDist(dist,flags,unfinished):
	minval = float("inf")
	minind = -1
	for i,node in enumerate(unfinished):
		if flags[i] == 0 and dist[node] < minval:
			minval = dist[node]
			minind = i
	return minind, minval


def BFS(edges,start):
	if edges == []:
		print "Edges empty"
		return None
	gSize = len(edges)
	prev = [-1 for i in range(gSize)]
	q = set()
	q.add(start)
	F = []
	U = range(gSize)
	
	while not len(q) == 0:
		curr = q.pop()
		
		for i in edges[curr]:
			if i not in F:
				q.add(i)
				prev[i] = curr
		
		U.remove(curr)
		F.append(curr)
	if len(U) > 0:	
		for node in U:
			prev[node] = float("inf")
	return prev

def singleBFS(edges,start,stop):
	
	if edges == []:
		print "Edges empty"
		return None
	gSize = len(edges)
	prev = [-1 for i in range(gSize)]
	q = set()
	q.add(start)
	F = []
	U = range(gSize)
	while not len(q) == 0:
		curr = q.pop()
		for i in edges[curr]:
			if i not in F:
				q.add(i)
				prev[i] = curr
		F.append(curr)
		U.remove(curr)
		if(curr == stop):
			break
	
	if len(U) > 0:	
		for node in U:
			prev[node] = float("inf")
	
	return prev


def dijkstra(edges,start):
	if edges == []:
		print "Edges empty"
		return None
	gSize = len(edges)
	
	dist = [float("Inf") for i in range(gSize)]
	prev = [-1 for i in range(gSize)]
	dist[start] = 0
	F = []
	U = range(gSize)
	flags = [0 for i in range(gSize)]
	
	while 0 in flags:
		minindex, minvalue = getMinDist(dist,flags,U)
		v = int(minindex)
		F.append(v)
		flags[v] = 1
		for v2, edge in enumerate(edges[v]):
			if dist[v] + 1 < dist[edge]:
				dist[edge] = dist[v] + 1
				prev[edge] = v
	return dist,prev
	
def dijkstraSingle(edges,start,stop):
	print "dijkastra single started"
	print edges[start]
	if edges == []:
		print "Edges empty"
		return None
	gSize = len(edges)
	dist = [float("Inf") for i in range(gSize)]
	prev = [-1 for i in range(gSize)]
	dist[start] = 0
	F = []
	U = range(gSize)
	flags = [0 for i in range(gSize)]
	
	while 0 in flags:
		minindex, minvalue = getMinDist(dist,flags,U)
		v = int(minindex)
		F.append(v)
		flags[v] = 1
		for v2, edge in enumerate(edges[v]):
			if dist[v] + 1 < dist[edge]:
				dist[edge] = dist[v] + 1
				prev[edge] = v
		#if target node has been reached, stop
		if v == stop:
			break
	print "dijkastra single finished"	
	return dist,prev
'''	
def calcAllShortestPath(edges):
	gSize = len(edges)
	shortestPath = [[[] for i in range(gSize)] for j in range(gSize)]
	for i in range(gSize):
		
		dist, prev = dijkstra(edges,i)
		
		for j in range(i+1,gSize):
			path = traceShortest(j,i,prev)
			shortestPath[i][j] = path
			shortestPath[j][i] = path
			
	return shortestPath
'''
def calcAllShortestPath(edges):
	gSize = len(edges)
	shortestPath = [[[] for i in range(gSize)] for j in range(gSize)]
	for i in range(gSize):
	
	
	#determine nodes that i has no connection to
		validEdges = allAbove(edges[i],i) #get all edges that have not been considered
		diff = set.difference(set(range(i+1,gSize)),set(validEdges)) #get all nodes that have not been considered yet and do not have an edge to i
		#one BFS for all in diff
		prev = BFS(edges,i)
		for j in diff:
			path = traceShortest(j,i,prev)
			shortestPath[i][j] = path
			shortestPath[j][i] = path
	
		#one BFS for each edge that i is connected to
		for j in edges[i]:
			if len(edges[i]) == 1 or len(edges[j]) == 1:
					shortestPath[i][j] = []
					shortestPath[j][i] = []
					continue

			edgeTMP = copy.deepcopy(edges)
			edgeTMP[i].remove(j)
			edgeTMP[j].remove(i)

			prev = singleBFS(edgeTMP,i,j)
			path = traceShortest(j,i,prev)
			shortestPath[i][j] = path
			shortestPath[j][i] = path

	return shortestPath

def allAbove(array,value):
	out = []
	for elem in array:
		if elem > value:
			out.append(elem)
	return out
				
def traceShortest(fr,to,prev):
	if prev[fr] == float("inf"):
		return []
	if prev[fr] == to:
		return [fr,to]
	if fr == to:
		return []
	path = []
	last = fr
	while True:
		path.append(last)
		last = prev[last]
		if last == to:
			break
	path.append(to)
	return path

def countSubgraphEdges(subg,gEdges):
	ecount = 0
	for n1 in subg:
		for n2 in subg:
			if n2 in gEdges[n1]:
				ecount += 1
	return ecount / 2

def calcSignMatrix(edges):
	shortestP = calcAllShortestPath(edges)
	
	gSize = len(edges)

	

	
	signM = [["" for i in range(gSize)] for j in range(gSize)]

	
	for i in range(gSize):
		for j in range(gSize):
			if i == j:
				signM[i][j] = "-0.1.0"
				continue
			#if shortestP[i][j] == None:
			#	print str(i) + " - " + str(j) + " <-"
			#	return 0
				
			if j in edges[i]:
				prefix = "+"
			else:
				prefix = "-"
			duv = len(shortestP[i][j])-1
			nuv = len(shortestP[i][j])
			muv = countSubgraphEdges(shortestP[i][j],edges)
			
			sign = prefix + str(duv) +"." + str(nuv)+"." + str(muv)
			
			
			signM[i][j] = sign
	return signM

#signList: vector of different signs found in alphabetical order
#freqMatrix list of frequency vectors 
def calcSignFrequencies(signM):
	signSet = set()
	freqList = []
	for i, row in enumerate(signM):
		freqList.append(freqInRow(row,signSet))
		
	signList = list(signSet)
	
	signList.sort()
	freqMatrix = []
	
	for freq in freqList:
		freqVector = [0 for i in range(len(signList))]
		
		for sign in freq.keys():
			ind = signList.index(sign)
			freqVector[ind] = freq[sign] #assign frequency value from hash map to correct position in frequency matrix
		freqMatrix.append(freqVector)
			
	
	reorder = sorted(range(len(freqMatrix)), key= lambda k: freqMatrix[k]) # compute reordering of frequency matrix in row major order

	makeOrder(freqMatrix,reorder)
	return signList,freqMatrix,reorder

"""
v:	vector to be ordered
order:	numerical vector that contains the desired position for each element in v
"""	
def makeOrder(v,order):
	if not len(v) == len(order):
		print "reorder was given 2 arguments with different length"
		return None
	
	for i, index in enumerate(order):
		tmp = v[i]
		v[i] = v[index]
		v[index] = tmp
	
#count element e +1 in hash h
def addToHash(h,e):
	if h.has_key(e):
		h[e] += 1
	else:
		h[e] = 1
		

def freqInRow(v,signSet):
	freq = {}
	for elem in v:
		addToHash(freq,elem)
		signSet.add(elem)
	return freq

'''
swaps two rows of a matrix while updating the indices of the reordering
'''
def swapRow(array,ind1,ind2,reorder):
	tmp = array[ind1]
	tmpInd = reorder[ind1]
	array[ind1] = array[ind2]
	array[ind2] = tmp
	reorder[ind1] = reorder[ind2]
	reorder[ind] = tmpInd

def getIsomorphism(edges1,edges2):
	if not len(edges1) == len(edges2):
		return 0
	gSize = len(edges1)
	
	signM1 = calcSignMatrix(edges1)
	signM2 = calcSignMatrix(edges2)
	signList1,signFreq1,reorder1 = calcSignFrequencies(signM1)
	signList2,signFreq2,reorder2 = calcSignFrequencies(signM2)
	if not signList1 == signList2:
		return 0
	
	for i, row1 in enumerate(signFreq1):
		if not signFreq1[i] == signFreq2[i]:
			flag = False
			for j in range(i+1,gSize):
				if signFreq1[i] == signFreq2[j]:
					flag = True
					swapRow(signFreq2,i,j,reorder2)
					break
			if not flag:
				print "There is no isomorphism, aborted at index " + i
				return 0
	return reorder1,reorder2
	
def checkAtomTypes(n1,n2,mapping):

	for i, at1 in enumerate(n1):
		print i
		if not n1[mapping[0][i]].element == n2[mapping[1][i]].element:
			return False
	return True
			

#nodes,edges = makeGraph(model)

#n1 ,e1 = makeGraph(m1)
#n2 ,e2 = makeGraph(m2)


#print e1
#print e2
#iso = getIsomorphism(e1,e2)
#print iso
#print checkAtomTypes(n1,n2,iso)

#dist, prev = dijkstra(nodes,edges,7)
#print dist
#print prev

