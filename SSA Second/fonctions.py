# Algorithme de simulation stochastique pour les métapopulations
import os
import classes
import numpy as np
import pandas as pd
import stochpy
from copy import deepcopy
from itertools import islice


def SetMetapop(nbsite, taillepop): #Creates sites objects containing populations
    ListSites=[] # List that will contain all sites
    for i in range(nbsite-1): # Creates sites, the 1st will always contain one infected
        if i == 0:
            newsite = classes.Site(effectifS=taillepop-1, effectifI=1)
            ListSites.append(newsite)
        else:newsite = classes.Site(effectifS=taillepop, effectifI=0)
        ListSites.append(newsite)
    print(ListSites)
    return ListSites