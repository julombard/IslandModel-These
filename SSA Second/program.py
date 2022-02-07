# Algorithme de simulation stochastique pour les métapopulations
import os
import classes
import fonctions
import numpy as np
import pandas as pd
import stochpy
from copy import deepcopy
from itertools import islice

#A story inspired by Modified Poisson Tau leap algorithm from cao et. al. (2005)
#Including New features for metapopulations, designed by Massol F., Lion S. and bibi
#Not Seen on TV !

#Simulation parameters
sim_time = 0 # Simulation time (model time, not an iteration number)
vectime = [0] # to keep t variable
tmax = 40 # Ending time
nbsite = 2 # Number de sites
Taillepop = 100 # Initial local population sizes

#Model parameters
beta = 0.005  #Infectious contact rate
r = 1.5  #Per capita growth rate
k = 1000  #Carrying capacity
d = 0.05  #Dispersal propensity
gamma = 1.5  #Parasite Clearance
alpha = 0.10  #Parasite Virulence
rho = 0.1 #Dispersal Cost
epsilon = 0.1 #Extinction rate

#Define population as class instances
ListSites = fonctions.SetMetapop(nbsite, Taillepop)

#Event definition
#Further expansion : build events from a unique .txt file read by the program, in order to simulate whathever you want
ReproductionS = classes.Event(propensity= 'r*self.S', Schange='1', Ichange='0', order=1)
DeathS = classes.Event(propensity='r*self.S*(self.S+self.I)/k', Schange='-1', Ichange='0', order=3)
DispersalS = classes.Event(propensity='d*self.S', Schange='-1', Ichange='0', order=1)
DispersalI = classes.Event(propensity='d*self.I', Schange='0', Ichange='-1', order=1)
Extinction = classes.Event(propensity='epsilon', Schange='-self.S', Ichange='-self.I', order=0)
Infection = classes.Event(propensity='beta*self.S*self.I', Schange='-1', Ichange='1', order=2)
Recovery = classes.Event(propensity='gamma*self.I', Schange='1', Ichange='-1', order=1)
DeathI = classes.Event(propensity='alpha*self.I', Schange='0', Ichange='-1', order=1)

#Event vector, cause tidying up things is nice
Events = [ReproductionS, DeathS, DeathI, DispersalI, DispersalS, Extinction, Infection, Recovery]
rouge = eval(2)
print(rouge)
print('coucou',ReproductionS.propensity)

Propensities = []
for i in Events: # For each event
    PropEvent = []
    for j in ListSites : # And each site
        S, I = j.effectifS, j.effectifI # Get the xi
        Prop = i.UpdatePropensity(S,I) #Compute propensity
        PropEvent.append(Prop)