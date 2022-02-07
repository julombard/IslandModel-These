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
    for i in range(nbsite): # Creates sites, the 1st will always contain one infected
        if i == 0:
            newsite = classes.Site(effectifS=taillepop-1, effectifI=1)
            ListSites.append(newsite)
        else:
            newsite = classes.Site(effectifS=taillepop, effectifI=0)
            ListSites.append(newsite)
    print(ListSites)
    return ListSites

def GetPropensites (Sites, Events):
    Propensities = []
    for i in Events: # For each event
        PropEvent =[]
        for j in Sites : # And each site
            S, I = j.effectifS, j.effectifI # Get the xi
            Prop = i.UpdatePropensity(S,I) #Compute propensity
            PropEvent.append(Prop)
        Propensities.append(PropEvent)
    sumpropensities = []
    for i in Propensities :
        sumpropensities.append(sum(i))
    return Propensities, sumpropensities

def SumDensities(Sites) :
    SumS = 0
    SumI = 0
    for i in Sites:
        SumS += i.effectifS
        SumI += i.effectifI
    return SumS, SumI

