from collections import defaultdict

#/model/chain/resNum/resName/atomNum/atomName/
fIdent=["recName","serialNum","atomName","alternateLoc","resName","chain","resNum","insCode","X","Y","Z","occupancy","temp","segid","element","charge"]
'''
turn a given range of numbers into a list containing the explicit numbers
e.g. "1-5,10,12-14" -> [1,2,3,4,5,10,12,13,14]
'''
def rangeAsVector(string):
	selection= []
	dat = string.split(",")
	for item in dat:
		if item.__contains__("-"):
			r = map(int,item.split("-"))
			selection.extend(map(str,[x for x in range(r[0],r[1]+1)]))
		else:
			selection.append(item)
	return selection

def getSelection(criteria,recordList):
	sel = []
	for i, record in enumerate(recordList):
		add = True
		for j in criteria.iterkeys():
			if str(record[j]) not in criteria[j]:
				add = False 
				#~ print "now false"
			#~ print j
			#~ print record[j]
			#~ print criteria[j]
		#~ print add
		#~ print record
		
		if add:
			sel.append(record)
	return sel

def select(selString,recordList):
	fieldInd = [5,6,4,1,2]
	criteria = selString.split("/")
	mSelect = criteria[0]
	
	toCheck = defaultdict()
	
	for i, crit in enumerate(criteria[2:7]):
		if len(crit) == 0:
			continue
		if any(c.isalpha() for c in crit):
			toCheck[fieldInd[i]] = crit.split(",")
		else:
			toCheck[fieldInd[i]] = rangeAsVector(crit)
	#~ print toCheck
	return getSelection(toCheck,recordList)
	
#~ def getSelRows(sel,rec):
	#~ 
	#~ res = []	
	#~ for i in sel:
		#~ res.append(rec)
	#~ return res 
