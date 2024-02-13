import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import numpy as np
import linecache


dirOut = '../output/'

#opens the file containing the eigenvectors
eigenVec = open(dirOut+'Eigenvectors/Eigenvectors3.dat', 'r')

#open file with the bonds
myfile = open(dirOut+'0explicit.inp', 'r')

#get the total number of bonds
dummy = linecache.getline(dirOut+'0explicit.inp', 5)
dummy = dummy.split()
nBond = int(dummy[0])  

#get the total number of atoms
dummy = linecache.getline(dirOut+'0explicit.inp', 3)
dummy = dummy.split()
nAtom = int(dummy[0])  

#empty image
G = nx.Graph()

#lines_to_read = [14, 208]
lines_to_read = range(14,14+nBond)

#create a list with the bonds readed from the file
for position, line in enumerate(myfile):
        if position in lines_to_read:
            bond=line.split()
            G.add_edge(bond[0], bond[1])

#close the file
myfile.close()

#Generates the positions
pos = nx.kamada_kawai_layout(G)
#pos = nx.spring_layout(G, iterations=1000,dim=2, seed=200)

lines_to_read = range(1,1+nAtom)

#creates a list with the corresponding coefficients
coeff = []
for position, line in enumerate(eigenVec):
        if position in lines_to_read:
            coeff.append(float(line))


eigenVec.close()

#looks for the maximum and minumum values of the coefficients
maxCoeff = max(coeff)
minCoeff = min(coeff)

#assign a color to each node based on its coefficient
color_map = []
for i in coeff:
    if i > 0 :
        #blue if its positive
        color_map.append((0,0,np.divide(i,maxCoeff)))
        #print(0,0,np.divide(i,maxCoeff))
    else:
        #red if its negative
        color_map.append((np.divide(i,minCoeff),0,0))
        #print(np.divide(i,minCoeff),0,0)
       
#options for the nodes and bonds
options = {
    "with_labels": False,
    "node_size": 80,
    "edgecolors": "black",
    "linewidths": 0,
    "width": 1,
}

#create hte image using the nodes and bonds with the options given
nx.draw_networkx(G, pos, node_color=color_map, **options)
#nx.draw_networkx(G, pos, **options)

# Set margins for the axes so that nodes aren't clipped
ax = plt.gca()
ax.margins(0.20)


#show the plot
plt.axis("off")
plt.show()


