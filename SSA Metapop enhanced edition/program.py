# Algorithme de simulation stochastique pour les métapopulations
from copy import deepcopy
import os
import classes
import fonctions
import numpy as np
import pandas as pd


#A story inspired by Modified Poisson Tau leap algorithm from cao et. al. (2005)
#Including New features for metapopulations, designed by Massol F., Lion S. and bibi
#Not Seen on TV !

#Simulation parameters
sim_time = 0 # Simulation time (model time, not an iteration number)
vectime = [0] # to keep t variable
tmax = 40 # Ending time
nbsite = 3 # Number de sites
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
ReproductionS = classes.Event(name='Reproduction S',propensity='r*self.S', Schange='1', Ichange='0', order=1)
DeathS = classes.Event(name='Death S',propensity='r*self.S*(self.S+self.I)/k', Schange='-1', Ichange='0', order=3)
DispersalS = classes.Event(name='Dispersal S',propensity='d*self.S', Schange='-1', Ichange='0', order=1)
DispersalI = classes.Event(name='Dispersal I',propensity='d*self.I', Schange='0', Ichange='-1', order=1)
Extinction = classes.Event(name='Extinction',propensity='epsilon', Schange='-self.S', Ichange='-self.I', order=0)
Infection = classes.Event(name='Infection',propensity='beta*self.S*self.I', Schange='-1', Ichange='1', order=2)
Recovery = classes.Event(name='Recovery',propensity='gamma*self.I', Schange='1', Ichange='-1', order=1)
DeathI = classes.Event(name='Death I',propensity='alpha*self.I', Schange='0', Ichange='-1', order=1)

#Event vector, cause tidying up things is nice
Events = [ReproductionS, DeathS, DeathI, DispersalI, DispersalS, Extinction, Infection, Recovery]

#Compute the propensities
Propensities, Sum_propensities = fonctions.GetPropensites(ListSites, Events) # Get a vector of propensities ordered by event and by sites

SumS, SumI = fonctions.SumDensities(ListSites)
print('Les props', Propensities)
print('Les sommes',Sum_propensities)
#Get Critical Reactions (maybe not useful since we use multinomial samples to identifiy sites where reactions occurs)
Criticals = fonctions.GetCriticals(Propensities, ListSites, Events)

#We now can compute vectors mu and sigma using previous shit
Mus = fonctions.ComputeMuNSigma(Sum_propensities, Events, ListSites) # As each statechange is 1 , 0, or -1 we have sigma = mu
print('les mumu', Mus)

#Get epsilon_i
Epsis = fonctions.GetEpsilonI(SumS, SumI)
print('Les ei*xi', Epsis)
#Get Tau prime
TauPrime = fonctions.GetTauPrime(Epsis, Mus)

# Now that main intermediary computations are done, let's get to the main algorithm Decision tree
aox = sum(Sum_propensities)

if TauPrime < 1/aox : # Take 10/aox 1 is left for ignoring this part
    print('Direct Method performed')
    pass # Insert direct method here
else:
    print('Lets leap baby')
    #Here we do not compute TauPrimePrime to determine how much critical reactions occurs
    #We expect that random sample of the place of reactions will be equivalent
    #As critical reactions in critical population will have low occurences
    Tau=TauPrime

    #So we directly sample the kjs (number of a given event during tau) from a poisson distribution
    Poisson_means = np.multiply(Tau,np.array(Sum_propensities)) #Get ajx * tau from which we will sample the kjs
    #print('Poisson', Poisson_means) # it's working we are so glad

    #Now we sample the kjs in poisson law, aka the number of trigger of each event
    triggers = []
    for i in Poisson_means :
        kj = np.random.poisson(i,1)
        triggers.append(kj[0]) # The [0] is due to array structure of kj
    print('Occurrences', triggers)
    #Now we sample the sites where events will occur from multinomial law
    #And apply the effect of event
    for index,event in enumerate(Events) : # For each event
        #This part define the number of occurences per site
        Noccur = triggers[index] #We get the number of times it should trigger during tau
        props = Propensities[index] # We get the propensity per sites
        #print('site propensities',props)
        SumProp = sum(props) # We get the total propensity
        Probas = [i /SumProp for i in props] # We get probability of occurence in each site
        #print('les probas', Probas) #Good job boy
        trigger_persite = np.random.multinomial(Noccur, Probas) # Working 1st try, you're the best coder ever or what ?
        #print('Occurrences per sites', trigger_persite)

        #This part apply the effect of events in site populations
        for index, Site in enumerate(ListSites) :
            if 'Dispersal' in event.name :
                # Multiply the state change in population by the number of triggers
                Site.effectifS += trigger_persite[index] * event.Schange
                Site.effectifI += trigger_persite[index] * event.Ichange
                nbmigrants = max(abs(trigger_persite[index] * event.Schange), abs(trigger_persite[index] * event.Schange))
                #print('Nombre de migrants', nbmigrants)
                for i in range(nbmigrants) :
                    #Determine which site will receive the dispersing individual
                    receiving_sites = deepcopy(ListSites)# Create a working copy of sites
                    #print('receivers', receiving_sites)
                    del receiving_sites[index] # removing departure site from the copy
                    #print('receivers post suppression', receiving_sites)
                    site_destination = np.random.choice(receiving_sites) # destination is a site object
                    #print('The destination is', site_destination)

                    #add individual to destination
                    if abs(event.Schange) > 0 : #if S are dispersers
                        site_destination.effectifS += 1
                    elif abs(event.Ichange) > 0 :
                        site_destination.effectifI += 1
                    else : print('ERROR : disperser is not S or I and that is very curious !')
            else:
                #Multiply the state change in population by the number of triggers
                Site.effectifS += trigger_persite[index]*event.Schange
                Site.effectifI += trigger_persite[index]*event.Ichange


