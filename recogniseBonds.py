
#map elements to covalent radii:
"""
covalrad= { \
"H" : 80, \
"He": 400, \
"Li": 170, \
"Be": 88, \
"B" : 208, \
"C" : 180, \
"N" : 170, \
"O" : 170, \
"F" : 160, \
"Ne": 280, \
"Na": 243, \
"Mg": 275, \
"Al": 338, \
"Si": 300, \
"P" : 259, \
"S" : 255, \
"Cl": 250, \
"Ar": 392, \
"K" : 332, \
"Ca": 248, \
"Mn": 338, \
"Fe": 335, \
"Cu": 380, \
"Zn": 362, \
"Se": 305, \
"Br": 302, \
"Ag": 398, \
"I": 350, \
"Au": 375, \
}
"""

covalrad= { \
"H" : 0.31, \
"He": 0.28, \
"Li": 1.28, \
"Be": 0.96, \
"B" : 0.84, \
"C" : 0.76, \
"N" : 0.71, \
"O" : 0.66, \
"F" : 0.57, \
"Ne": 0.58, \
"Na": 1.66, \
"Mg": 1.41, \
"Al": 1.21, \
"Si": 1.11, \
"P" : 1.07, \
"S" : 1.05, \
"Cl": 1.02, \
"Ar": 1.06, \
"K" : 2.03, \
"Ca": 1.76, \
"Mn": 1.39, \
"Fe": 1.32, \
"Cu": 1.32, \
"Zn": 1.22, \
"Se": 1.20, \
"Br": 1.20, \
"Ag": 1.45, \
"I": 1.39, \
"Au": 1.36, \
}
def rasMol(dist,atoms):
    if len(atoms) >255:
        return ras1(dist,atoms)
    else:
        return ras2(dist,atoms)

# for >255 molecules
def ras1(dist,atoms):
    bonds = [[] for i in range(len(dist))]
    for i, line in enumerate(dist):
        for j, entry in enumerate(line[:i]):
            if atoms[i].element == "H" or atoms[j].element == "H":
                if entry > 0.4 and entry < 1.2:
                    bonds[i].append(j)
                    bonds[j].append(i)
            else:
                if entry > 0.4 and entry < 1.9:
                    bonds[i].append(j)
                    bonds[j].append(i)
    return bonds

#for <=255
def ras2(dist,atoms):
    bonds = [[] for i in range(len(dist))]
    covradii= []
    for atom in atoms:
        rad = covalrad[atom.element]
        covradii.append(rad) #convert to angstroem and append
        
    for i, line in enumerate(dist):
        for j, entry in enumerate(line[:i]):
            maxdist = covradii[i] + covradii[j] + 0.56
            if entry > 0.4 and entry < maxdist:
                bonds[i].append(j)
                bonds[j].append(i)
    return bonds       

            
