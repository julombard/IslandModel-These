# Algorithme de simulation stochastique pour les métapopulations
import numpy as np
import stochpy
from copy import deepcopy

stochpy.SSA
# Implémentation du modèle
# Même principe que pour une population isolée MAIS
# On a N populations
# On tire les évènements pour UN individu parmis TOUTES les populations
# Si l'évènement est "migration" on doit définir un site d'arrivée, puis mettre sa population à jour

# Modified Poisson Tau leap algorithm from cao et. al. (2005)


#Définition des paramètres du modèle
# InitPar
beta = 0.005  # Infectious contact rate
r = 1.5  # Per capita growth rate
k = 1000  # Carrying capacity
d = 0.05  # Dispersal propensity
gamma = 1.5  # Clearance
alpha = 0.10  # Virulence
rho = 0.1 # Dispersal Cost



######## Step 1 : Fabriquer une métapop pour un état initial donné
def Set_metapop(taillepop, nbsites): #Initialise une métapopulation de nbsites sites contenant chacun taillepop individus, un site aléatoire contiendra une individu infecté
    Metapop=[]

    for i in range(nbsites):
        if i == 0: Metapop.append([taillepop-1,1])
        else :Metapop.append([taillepop,0])
    Metapoparray = np.array(Metapop) # Les array sont stockés dans des blocs mémoire contigus et ca devrait nous accélérer entre 'un peu' et 'pas mal'
    return Metapoparray

Metapop = Set_metapop(100, 10) # Définir une métapop en moins de temps qu'il n'en faut pour définir la fonction !


def Tauleap4Metapop() : # Will be used someday, or maybe not
    return 0

def IsCritical(xi,propensity) : # Détermine si une réaction est critique (renvoie un booléen)
    # Ncrit est défini un peu à la louche, Cao dit " entre 2 et 20 " , donc on va partir sur 11
    # Dans ce modèle on a toujours Lj = xi donc c'est pratique
    Ncrit = 11
    iscritical = 0
    if propensity > 0 and xi <= Ncrit : iscritical = 1
    else : pass
    return iscritical
def Get_xi(Metapop) : # Get density values of each subpoppulation of each site
    list_xi = []
    for pop in Metapop :
        list_xi.append(pop[0])
        list_xi.append(pop[1])
    xi = np.array(list_xi)
    return xi

def ComputePropensitiesAndCriticals(Metapop): # Compute propensities for all events
    #Metapop doit être un array de taille quelconque qui contient des listes de taille 2 [S,I]

    list_Propensities = []
    list_criticals = []

    for pop in Metapop : # Pour chaque population
        S = pop[0] # ON récupère son nombre de susceptibles
        I = pop[1] # Idem pour les infectés

        # On passe en revue tout ce qui peut se produire et on en détermine les chances
        # prop x = proba, et ensuite on ajoute la propensity à une liste qui les contiendra toutes
        # Avec la valeur de xi et de propensity, on détermine si la réaction est critique

        # Reproduction S
        prop1 = r*S
        list_Propensities.append(prop1)
        crit = IsCritical(S, prop1)
        list_criticals.append(crit)
        # Mort S
        prop2 = r*S*(S+I)/k
        list_Propensities.append(prop2)
        crit = IsCritical(S, prop2)
        list_criticals.append(crit)
        # Infection
        prop3=beta * S * I
        list_Propensities.append(prop3)
        crit = IsCritical(S, prop3)
        list_criticals.append(crit)
        #MigrationS
        prop4=d * S
        list_Propensities.append(prop4)
        crit = IsCritical(S, prop4)
        list_criticals.append(crit)
        # Guerison:
        prop5=gamma * I
        list_Propensities.append(prop5)
        crit = IsCritical(I, prop5)
        list_criticals.append(crit)
        # Mort I
        prop6=alpha * I
        list_Propensities.append(prop6)
        crit = IsCritical(I, prop6)
        list_criticals.append(crit)
        # Migration I
        prop7=d * I
        list_Propensities.append(prop7)
        crit = IsCritical(I, prop7)
        list_criticals.append(crit)
    propensities = np.array(list_Propensities) # On transvase dans un array parce que c'est mieux
    Sum_propensities = np.sum(propensities)
    Criticals = np.array(list_criticals)

    return propensities, Sum_propensities, Criticals

Propensities, SumProp, Criticals = ComputePropensitiesAndCriticals(Metapop)


def GetTauPrime(xi, propensity, criticals):



    return 0
