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

def GetCriticals(Propensities, Sites, Events):
    Crit_treshold = 11 #Critical number of individuals
    Criticals = []
    for indexi ,i in enumerate(Events) :
        CriticalEvent = []
        for indexj,j in enumerate(Sites) :
            S, I = j.effectifS, j.effectifI  # Get the xi
            print('formule', i.formula)
            if 'epsilon' in i.formula : # Case of extinction which is always critic
                CriticalEvent.append(1)
                continue
            if i.Schange < 0 : # Case where an S individual is depleted
                if S < Crit_treshold : # If S subpop is low
                    if Propensities[indexi][indexj] > 0 : # But not zero
                        CriticalEvent.append(1) # Its critical
                    else: CriticalEvent.append(0) # Its not critical
                else : CriticalEvent.append(0)
            if i.Ichange < 0 : # Case where an I individual is depleted
                if I < Crit_treshold : # If I subpop is low
                    if Propensities[indexi][indexj] > 0 : # But not zero
                        CriticalEvent.append(1) # Its critical
                    else : CriticalEvent.append(0) # Its not critical
                else: CriticalEvent.append(0)
            elif i.Schange >0 and i.Ichange == 0 : CriticalEvent.append(0) # Cas reproduction S
        Criticals.append(CriticalEvent)
    return Criticals

def ComputeMuNSigma(SumPropensities , Events, Sites):
    MuS_vector = []
    MuI_vector = []

    for index, i in enumerate(Events) :
        if i == 'epsilon':
            continue
        elif i.Schange < 0 : # Case where an S individual is depleted
            vij = i.Schange
            PropS = SumPropensities[index]
            MuS_vector.append(vij*PropS)
        elif i.Ichange < 0 :
            vij = i.Schange
            PropI = SumPropensities[index]
            MuI_vector.append(vij * PropI)
        elif i.Schange > 0 and i.Ichange == 0 :
            vij = 1
            PropS = SumPropensities[index]
            MuS_vector.append(vij*PropS)
    MuS = sum(MuS_vector)
    MuI = sum(MuI_vector)
    Mu = np.array(MuS, MuI)
    return Mu

def GetEpsilonI(SS, SI):
    #Hor_S = 3 # Plus grand ordre de réaction consommant S
    #Hor_I = 2 #Idem pour I
    epsilon2 = 0.03
    gi_s = 3/2*(2+1/(SS-1))
    gi_i = 2
    eixi_s = epsilon2/gi_s * SS
    eixi_i = epsilon2/gi_i * SI
    EpsilonIs = np.array(eixi_s, eixi_i)
    return EpsilonIs