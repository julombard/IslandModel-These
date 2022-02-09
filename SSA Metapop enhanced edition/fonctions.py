# Algorithme de simulation stochastique pour les métapopulations
import classes
import numpy as np
from copy import deepcopy
rho = 0.1 # Not very convenient to put it here

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
    Output_vector = []

    for index, i in enumerate(Events) :
        if i.name == 'Extinction':
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
    Output_vector.append(sum(MuS_vector))
    Output_vector.append(sum(MuI_vector))
    out = np.array(Output_vector)

    return out

def GetEpsilonI(SS, SI):
    #Hor_S = 3 # Highest order reaction consuming S
    #Hor_I = 2 #Idem for I
    epsilon2 = 0.03
    Output_vector = []
    gi_s = 3/2*(2+1/(SS-1)) # Given in the paper
    gi_i = 2 # same
    Output_vector.append(epsilon2/gi_s * SS)
    Output_vector.append(epsilon2/gi_i * SI)

    out = np.array(Output_vector)
    return out

def GetTauPrime(Epsis, Mus):
    print('Coucou attention hahah', Epsis, Mus)
    Tau_candidates = []
    Tau_finalists =[]
    for i in range(len(Mus)): # Get all Tau candidates (2 per specie)
        candidates = []
        upperterm = max(Epsis[i], 1) # Numerator
        upperterm_squared = upperterm**2 # Squared

        # Cases where all infected died in the early sim, we use a very small value for mu
        # (which leads to a huge value for tau which will never be selected but allows program to continue)
        #if Mus[i] == 0 :
            #Mus[i]= 1 / pow(10, 6)
            #print('WARNING : Avoiding division by zero by doing Trickster shit')
        candidates.append(upperterm/ abs(Mus[i]))
        candidates.append(upperterm_squared/abs(Mus[i])) # No need to square because mu = sigma
        Tau_candidates.append(candidates)
    print('Tau candidates', Tau_candidates)
    for i in range(len(Tau_candidates)): # get the tau 'finalists' (1 per specie)
        Tau_finalists.append(min(Tau_candidates[i]))
    TauPrime = min(Tau_finalists) # TauPrime is the minimum of all
    #print('Tau candidats', Tau_candidates)
    #print('Tau finalistes', Tau_finalists)
    #print('and the winner is', TauPrime)

    return TauPrime

def DoDirectMethod(Propensities, Sum_propensities, exactsteps, events, sites):
    for i in range(exactsteps):
        r1 = np.random.uniform(0,1,1) #Draw random numbers
        a0 = sum(Sum_propensities) # Overall propensity

        Tau = (1/a0) * np.log(1/r1) #Time increment
        Probas = [i / a0 for i in Sum_propensities] #Individual pro
        NextReaction = list(np.random.multinomial(1, Probas)) # List to get index easily
        NextReactionIndex = NextReaction.index(1)
        #print('The next reaction Index', NextReactionIndex)

        #Determine where the reaction is going to happen
        props = Propensities[NextReactionIndex]  # We get the propensity per sites
        sumprop = sum(props)
        proba_site = [i/sumprop for i in props]
        NextPlace = list(np.random.multinomial(1, proba_site))
        NextPlaceIndex = NextPlace.index(1)
        #print('Next places', NextPlaceIndex)


        # This part apply the effect of events in site populations
        event = events[NextReactionIndex]
        site = sites[NextPlaceIndex]
        if 'Dispersal' in event.name:
            # Multiply the state change in population by the number of triggers
            site.effectifS +=  event.Schange
            site.effectifI +=  event.Ichange
            nbmigrants = 1
            # Here we apply dispersal cost to determine the number of successful migrants, rho is defined at the top
            SuccessfulMigrants = 0

            roll4urlife = np.random.uniform(0, 1, 1)
            if roll4urlife > rho: SuccessfulMigrants += 1
            # Here we distribute successful migrants among neighboring sites
            # This part can be improved as neighboring rules become more complex, using a specific class 'network' to determine the neighbors
            if SuccessfulMigrants == 1:
                # Determine which site will receive the dispersing individual
                receiving_sites = deepcopy(sites)  # Create a working copy of sites
                # print('receivers', receiving_sites)
                del receiving_sites[NextPlaceIndex]  # removing departure site from the copy
                # print('receivers post suppression', receiving_sites)
                site_destination = np.random.choice(receiving_sites)  # destination is a site object
                # print('The destination is', site_destination)

                # add individual to destination
                if abs(event.Schange) > 0:  # if S are dispersers
                    site_destination.effectifS += 1
                elif abs(event.Ichange) > 0:
                    site_destination.effectifI += 1
                else:
                    pass
                    #print('There was only one migrant, but he died. Nothing happend')
        else:
            # Multiply the state change in population by the number of triggers
            site.effectifS +=  event.Schange
            site.effectifI +=  event.Ichange
        return Tau