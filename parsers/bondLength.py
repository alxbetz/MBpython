	
def getBondsMap():
	f = open("./data/bondlength.txt","r")
	header = []
	l = 0
	bondLength = {}
	for i,line in enumerate(f):
		if i==0:
			header = line.strip("\n").split()
			l = len(header)
			continue
		if line.startswith("#") or line.strip("\n").strip() == "":
			continue
		line = line.strip("\n").split()
		for j,entry in enumerate(line[1:]):
			f = header[l-i]
			t = header[j]
			forw = f + "-" + t
			
			back = t + "-" + f
			entry = (float(entry) / 100) * 1.2 #convert pm to Angstrom and allow 10% deviation from mean
			bondLength[forw] = entry
			bondLength[back] = entry			
	return bondLength

