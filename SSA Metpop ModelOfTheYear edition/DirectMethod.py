#Direct method from Gillespie D. (2001)
#Outputs used to compare accuracy with Tau-leap adapted to metapopulations

import fonctions
import classes
import pandas as pd
import numpy as np

#Simulation parameters
np.random.seed(2) #Set seed for reproducibility
sim_time = 0 # Simulation time (model time, not an iteration number)
vectime = [0] # to keep t variable
tmax = 40 # Ending time
Nexactsteps = 20  # Number of steps to do if/when performing direct method
nbsite = 20 # Number de sites
Taillepop = 100 # Initial local population sizes

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
#Further expansion idea : build events from a unique model.txt file read by the program, in order to simulate whathever you want
ReproductionS = classes.Event(name='Reproduction S',propensity='r*self.S', Schange='1', Ichange='0', order=1)
DeathS = classes.Event(name='Death S',propensity='r*self.S*(self.S+self.I)/k', Schange='-1', Ichange='0', order=3)
DispersalS = classes.Event(name='Dispersal S',propensity='d*self.S', Schange='-1', Ichange='0', order=1)
DispersalI = classes.Event(name='Dispersal I',propensity='d*self.I', Schange='0', Ichange='-1', order=1)
Extinction = classes.Event(name='Extinction',propensity='epsilon', Schange='-self.S', Ichange='-self.I', order=0)
Infection = classes.Event(name='Infection',propensity='beta*self.S*self.I', Schange='-1', Ichange='1', order=2)
Recovery = classes.Event(name='Recovery',propensity='gamma*self.I', Schange='1', Ichange='-1', order=1)
DeathI = classes.Event(name='Death I',propensity='alpha*self.I', Schange='0', Ichange='-1', order=1)

#Event vector, cause tidying up things is always nice
Events = [ReproductionS, DeathS, DeathI, DispersalI, DispersalS, Extinction, Infection, Recovery]

#Initializing outputs storage
Densities_out = [] # Collect densities outputs
Propensities_out =[] # Collect propensities outputs
IsTrackPropensites = False #Set false/ true to not track/track propensities
for i in range(len(Events)):
    Propensities_out.append([0]) #Store times series of propensities

#Get initial lists and values in outputs
#We want to get one list per Sx(t) and Ix(t) to store them easily in a dataframe at the end
for i in ListSites:
    Densities_out.append([i.effectifS])
    Densities_out.append([i.effectifI])

while sim_time < tmax :
    print('We have currently passed', sim_time,'time in the simulation') # Kind of a loading bar
    vectime.append(sim_time) #Update time vector

    # Creating a vector with sites indexes, used later and put here to not be computed at each subloop iteration
    sites_indexes = []
    for i in range(len(ListSites)):
        sites_indexes.append(i)

    # Compute the propensities
    Propensities, Sum_propensities = fonctions.GetPropensites(ListSites,Events)  # Get a vector of propensities ordered by event and by sites
    SumS, SumI = fonctions.SumDensities(ListSites)  # Get total densities of species
    Tau = fonctions.DoDirectMethod(Propensities, Sum_propensities, Nexactsteps, Events, ListSites)

    #Update time
    sim_time += Tau

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
    if IsTrackPropensites == True :
        Sum_propensities # Propensities of each event in a list sorted by event
        for index,propensitiy in enumerate(Sum_propensities) :
            Propensities_out[index].append(propensitiy)


    #Not done yet but define a boolean IsTrackPropensity in sim options and track ioi = 1

    # Break the main loop if there are no infected remaining ( this happens essentially at start if the 1st infected dies)
    if SumI == 0:
        print('WARNING : ABORTED SIMULATION, No infected remaining')
        break

#Structuring outputs to get a .csv file event if the loop has broken
if IsTrackPropensites == True : #if we track propensities
    #Creating propensities dataframe
    dataprop = pd.DataFrame(columns=['t'])
    for event in Events :
        EventName = event.name
        dataprop[EventName]= []
    #Filling the dataframe
    for index, colname in enumerate(dataprop):
        if index == 0 : dataprop[colname] = vectime
        else : dataprop[colname] = Propensities_out[index-1]
    #Saving into .csv file
    dataprop.to_csv('DirectMethod_Propensities_outputs.csv')


#Creating the time series dataframe
data = pd.DataFrame(columns=['t'])
for i in range(nbsite): # Define a column for each subpop and site index
    colname_s = 'S'+str(i)
    colname_i = 'I'+str(i)
    data[colname_s] = []
    data[colname_i] = []
#Filling the dataframe
for index,colname in enumerate(data):
    if index == 0 : data[colname] = vectime
    else : data[colname] = Densities_out[index-1] # I don't remember why i put -1, but it doesn't work if it's removed
print(data)
#Saving into .csv file
data.to_csv('DirectMethod_Metapop_Outputs.csv')
# Strange, csv is saved in model directory on lab computer but in User/Appdata/local/temp on personal computer...
#Edit : it seems that computer is totally 'en roue libre' concerning where to save the files