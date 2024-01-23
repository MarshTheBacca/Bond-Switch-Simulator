import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

cmap = mpl.colormaps['viridis']


with open('Many_Atom_Test/test_e_compare.out', 'r') as f:
    array = np.genfromtxt(f)


norm = mpl.colors.Normalize(min(array[:,0]), max(array[:,0]))


for i in range(1,array.shape[1]):
    plt.scatter(array[:,0], array[:,i])
plt.show()

lable = ["p6", "Simple Graphene", "Tersoff Graphene", "Triangle Raft", "Bilayer"]

for i in range(1,array.shape[1]):
    for j in range(i+1, array.shape[1]):
#        for k in range(array.shape[0]):
#            if int(array[k,i])!=0 and int(array[k,j]!=0):
        plt.scatter(array[:,i], array[:,j], color=cmap(norm(array[:,0])))
        plt.title(lable[i] +" vs "+lable[j])
        plt.show()



