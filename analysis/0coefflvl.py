import os
import linecache

#os.system('ls Eigenvectors*.dat > EigenvectorsFiles')

os.system('touch 0coefflvl.dat')
os.system('rm 0coefflvl.dat')
os.system('touch 0coefflvl.dat')

os.system('touch 0coeffFiles.dat')
os.system('rm 0coeffFiles.dat')
os.system('touch 0coeffFiles.dat')

#eigenvecFiles = open('EigenvectorsFiles', 'r')
explicitFile = open('../0explicit.inp', 'r')

#get the total number of atoms
dummy = linecache.getline('../0explicit.inp', 3)
dummy = dummy.split()
nAtom = int(dummy[0])

eigenvecFiles=[]
for a in range(1,nAtom+1):
    name="Eigenvectors"+str(a)+".dat"
    eigenvecFiles.append(name) 

#get the total number of bonds
dummy = linecache.getline('../0explicit.inp', 5)
dummy = dummy.split()
nBond = int(dummy[0])

lines_to_read = range(14,14+nBond)

res = open('0coefflvl.dat', 'a')

fileMoreNodes = open('0coeffFiles.dat', 'a')

Bonds = []
for position, line in enumerate(explicitFile):
    if position in lines_to_read:
        dummy=line.split()
        Bonds.append(( int(dummy[0]), int(dummy[1])))

lvl = 0

for fileN in eigenvecFiles:

    lvl = lvl + 1

    fileN = fileN.strip()
    eigenvectorN = open(fileN, 'r')

    ener = linecache.getline(fileN,1)
    ener = ener.split()
    ener = ener[2]

    coeff = []
    lines_to_read = range(1,nAtom+1) 
    for position, line in enumerate(eigenvectorN): 
        if position in lines_to_read:
            coeff.append(float(line.strip()))

    cont=0
    for vec in coeff:
        if abs(vec) < 0.0000000001 :
            cont = cont + 1

    resultado = str(lvl) + ' ' +str(ener) + ' ' + str(cont)
    res.write(resultado)
    res.write('\n')

    if cont > 10 :
        fileMoreNodes.write(fileN)
        fileMoreNodes.write('\n')

    eigenvectorN.close()

#eigenvecFiles.close()

#os.system('rm EigenvectorsFiles')
