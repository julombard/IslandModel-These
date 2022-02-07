# Algorithme de simulation stochastique pour les métapopulations
import os
import main
import numpy as np
import pandas as pd
import stochpy
from copy import deepcopy
from itertools import islice

#Parametres
sim_time = 0 # Temps de simulation (vrai temps, pas nombre iterations)
vectime = [0] # Vecteur de stockage de variable t
tmax = 40 # Temps au bout duquel on s'arrête
nbsite = 2 # Nombre de sites

#Definition des populations comme instances de classe
ListSites=[] # Liste contenant l'ensemble des objets Sites
for i in range(len(nbsite)): # Crée les sites, le premier contient toujours 1 infecté
    if i == 0:
        newsite = Site(effectifS=99, effectifI=1)
    else:newsite = Site(effectifS=100, effectifI=0)
