import numpy as np
import stochpy
from copy import deepcopy

def Set_metapop(taillepop, nbsites): #Initialise une métapopulation de nbsites sites contenant chacun taillepop individus, un site aléatoire contiendra une individu infecté
    Metapop=[]

    for i in range(nbsites):
        if i == 0: Metapop.append([taillepop-1,1])
        else :Metapop.append([taillepop,0])
    Metapoparray = np.array(Metapop) # Les array sont stockés dans des blocs mémoire contigus et ca devrait nous accélérer entre 'un peu' et 'pas mal'
    return Metapoparray

Metapopulation = Set_metapop(100,10)
print(Metapopulation)