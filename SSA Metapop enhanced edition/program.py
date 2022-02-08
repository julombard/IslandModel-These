# Stochastic Simulation Algorithm for Metapopulation models
from copy import deepcopy
import os
import classes
import fonctions
import numpy as np
import pandas as pd


#A story inspired by Modified Poisson Tau leap algorithm from cao et. al. (2006)
#Including New features for metapopulations modelling, designed by Massol F., Lion S. and bibi
#Not Seen on TV (and will never be) !
#Damn efficient compared to previous try

#Simulation parameters
sim_time = 0 # Simulation time (model time, not an iteration number)
vectime = [0] # to keep t variable
tmax = 40 # Ending time
Nexactsteps = 20  # Number of steps to do if/when performing direct method
nbsite = 10 # Number de sites
Taillepop = 100 # Initial local population sizes
Densities_out = [] # Collect densities outputs

#Model parameters
beta = 0.05  #Infectious contact rate
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

#Get initial lists and values in outputs
#We want to get one list per Sx(t) and Ix(t) to store them easily in  dataframe at the end
for i in ListSites:
    Densities_out.append([i.effectifS])
    Densities_out.append([i.effectifI])
#print('Bonjour les amish', Densities_out)
#Get ready for the simulation
#Further expansion : import multiprocessing to launch as many run as CPU in the same time

#Main Loop
while sim_time < tmax :
    print('We have currently passed', sim_time,'time in the simulation')
    vectime.append(sim_time) #Update time vector

    #Compute the propensities
    Propensities, Sum_propensities = fonctions.GetPropensites(ListSites, Events) # Get a vector of propensities ordered by event and by sites
    SumS, SumI = fonctions.SumDensities(ListSites)
    #print('Les props', Propensities)
    #print('Les sommes',Sum_propensities)

    #Get Critical Reactions (maybe not useful since we use multinomial samples to identifiy sites where reactions occurs)
    #Criticals = fonctions.GetCriticals(Propensities, ListSites, Events)

    #We now can compute vectors mu and sigma using previous shit
    Mus = fonctions.ComputeMuNSigma(Sum_propensities, Events, ListSites) # As each statechange is 1 , 0, or -1 we have sigma = mu
    #print('les mumu', Mus)

    #Get epsilon_i
    Epsis = fonctions.GetEpsilonI(SumS, SumI)
    #print('Les ei*xi', Epsis)
    #Get Tau prime
    TauPrime = fonctions.GetTauPrime(Epsis, Mus)

    # Now that main intermediary computations are done, let's get to the main algorithm Decision tree
    aox = sum(Sum_propensities)

    if TauPrime < 10/aox : # Take 10/aox 1 is left for ignoring this part
        print('Direct Method performed')
        Tau = fonctions.DoDirectMethod(Propensities,Sum_propensities,Nexactsteps, Events,ListSites)
    else:
        print('Lets leap baby')
        #Here we do not compute TauPrimePrime to determine how much critical reactions occurs
        #We expect that random sample of the place of reactions will be equivalent
        #As critical reactions in critical population will have low occurences
        Tau=TauPrime
        #print('Voici Tau', Tau)
        #So we directly sample the kjs (number of a given event during tau) from a poisson distribution
        Poisson_means = np.multiply(Tau,np.array(Sum_propensities)) #Get ajx * tau from which we will sample the kjs
        #print('Poisson', Poisson_means) # it's working we are so glad

        #Now we sample the kjs in poisson law, aka the number of trigger of each event
        triggers = []
        for i in Poisson_means :
            kj = np.random.poisson(i,1)
            triggers.append(kj[0]) # The [0] is due to array structure of kj
        #print('Occurrences', triggers)
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

            if Noccur == 0 : #Case where the event can't happen
                trigger_persite=[0 for i in range(nbsite)]
            else : # Normal cases
                trigger_persite = np.random.multinomial(Noccur, Probas) # Working 1st try, you're the best coder ever or what ?
            #print('Occurrences per sites', trigger_persite)
            #This part apply the effect of events in site populations
            for index, Site in enumerate(ListSites) :
                if 'Dispersal' in event.name :
                    # Multiply the state change in population by the number of triggers
                    Site.effectifS += trigger_persite[index] * event.Schange
                    Site.effectifI += trigger_persite[index] * event.Ichange
                    nbmigrants = max(abs(trigger_persite[index] * event.Schange), abs(trigger_persite[index] * event.Ichange))
                    #print('Nombre de migrants', nbmigrants)

                    #Here we apply dispersal cost to determine the number of successful migrants, rho is defined at the top
                    SuccessfulMigrants = 0
                    for i in range(nbmigrants):
                        roll4urlife = np.random.uniform(0,1,1)
                        if roll4urlife > rho : SuccessfulMigrants += 1
                    #Here we distribute successful migrants among neighboring sites
                    #This part can be improved as neighboring rules become more complex, using a specific class 'network' to determine the neighbors
                    for i in range(SuccessfulMigrants) :
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
                            print('PYCHARM WAS HERE')
                            site_destination.effectifI += 1
                            print('le migrant a til été receptionné ?',site_destination.effectifI)
                        else : print('ERROR : disperser is neither S nor I and that is very curious !')
                else:
                    #Multiply the state change in population by the number of triggers
                    Site.effectifS += trigger_persite[index]*event.Schange
                    Site.effectifI += trigger_persite[index]*event.Ichange

    #Update time
    sim_time += Tau
    #print('time increment', Tau)
    #print('Le temps passe si vite',sim_time)

    #Update the output tracking
    #1. Densities
    indexlist = 0
    for i in ListSites :
        if i.effectifS < 0 : #Avoid negative population in the "big fat brute" way
            i.effectifS = 0
        Densities_out[indexlist].append(i.effectifS)
        indexlist += 1
        if i.effectifI<0:
            i.effectifI = 0
        Densities_out[indexlist].append(i.effectifI)
        indexlist += 1
    #2. Propensities


#Structuring outputs to get a .csv file

print('boubou',len(Densities_out))
#Creating the dataframe
data = pd.DataFrame(columns=['t'])
for i in range(nbsite):
    colname_s = 'S'+str(i)
    colname_i = 'I'+str(i)
    data[colname_s] = []
    data[colname_i] = []
print('the DF', data)

for index,colname in enumerate(data):
    if index == 0 : data[colname] = vectime # This one is OK
    else : data[colname] = Densities_out[index-1]
print(data)



