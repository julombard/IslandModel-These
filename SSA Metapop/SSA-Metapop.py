# Algorithme de simulation stochastique pour les métapopulations
import numpy
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

# La sélection de tau est basée de "Efficient step size selection for tau leaping simulation", cao et al (2006)


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
print('mes couilles', Metapop)

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

#Propensities, SumProp, Criticals = ComputePropensitiesAndCriticals(Metapop)


def GetTauPrime(xi, propensity, criticals):
    # Traduction de l'équation mathématique en français
    # Pour chaque évènement NON critique
    # On calcule pour chaque espèce la dérivée de la fonction de propensity par rapport a l'espèce, qu'on multiplie
    # par le vecteur de changement d'état


    #Le résultat est remultiplié par la propensity pour avoir mu


    #Le carré du résultats est multiplié par la propensity pour avoir sigma


    # Le choix du tau est le minimum entre deux trucs, l'un fait à partir de mu, l'autre de sigma


    return 0



def GetMainMatrix(Metapop) : # Fonction qui récupère la matrice des propensities et celle qui contient les vecteurs de changement d'états

    #Séparation des S et I par convenance
    Slist = []
    Ilist = []

    for pop in Metapop :
        Slist.append(pop[0])
        Ilist.append(pop[1])
    Sarray = np.array(Slist)
    print('Les S', Sarray)
    Iarray= np.array(Ilist)


    StateChangeMatrixS = np.empty(1)
    StateChangeMatrixI = np.empty(1)
    MatrixPropensitiesS = np.empty(1)
    MatrixPropensitiesI = np.empty(1)
    NbPops = len(Metapop)
    NbSpecies = len(Metapop)*len(Metapop[0])
    print('Nombre de sites', NbPops)
    print('Nombre entités', NbSpecies)

    for specie in range(len(Sarray)) :
        propensities_I = np.zeros(NbSpecies)
        StatechangeVector_I = np.zeros(NbSpecies)

        propensities_S = np.zeros(NbSpecies)
        StatechangeVector_S = np.zeros(NbSpecies)

        # On rempli deux array de taille Nbspecie
        #L'un contient les propensity de chaque évènement
        #L'autre contient les changements d'états

        S = Sarray[specie]
        I = Iarray[specie]

        ##### On décline tous les évènements possibles pour remplir les changements d'états associés + propensity
        #Reproduction S
        prop = r* S
        propensities_S[specie] = prop
        np.append(MatrixPropensitiesS, propensities_S)

        changestate = +1
        StatechangeVector_S[specie] = changestate
        np.append(StateChangeMatrixS, StatechangeVector_S)

        #Mort S
        prop = r * S * (S+I) /k
        propensities_S[specie] = prop
        np.append(MatrixPropensitiesS, propensities_S)

        changestate = -1
        StatechangeVector_S[specie] = changestate
        np.append(StateChangeMatrixS, StatechangeVector_S)

        #Migration S
        prop = d * S
        propensities_S[specie] = prop
        np.append(MatrixPropensitiesS, propensities_S)

        changestate = -1
        StatechangeVector_S[specie] = changestate
        np.append(StateChangeMatrixS, StatechangeVector_S)

        #Migration I

        prop = d * I
        propensities_I[specie] = prop
        np.append(MatrixPropensitiesI, propensities_I)

        changestate = -1
        StatechangeVector_I[specie] = changestate
        np.append(StateChangeMatrixI, StatechangeVector_I)

        # Mort I
        prop = alpha * I
        propensities_I[specie] = prop
        np.append(MatrixPropensitiesI, propensities_I)

        changestate = -1
        StatechangeVector_I[specie] = changestate
        np.append(StateChangeMatrixI, StatechangeVector_I)

        #Guérison I

        prop = alpha * I
        propensities_I[specie] = prop
        np.append(MatrixPropensitiesI, propensities_I)

        changestateS = -1
        changestateI = +1
        StatechangeVector_I[specie] = changestateI
        np.append(StateChangeMatrixI, StatechangeVector_I)

        StatechangeVector_S[specie] = changestateS
        np.append(StateChangeMatrixS, StatechangeVector_S)

        #Infection

        prop = beta * S * I
        propensities_S[specie] = prop
        np.append(MatrixPropensitiesS, propensities_I)

        changestateS = -1
        changestateI = +1

        StatechangeVector_I[specie] = changestateI
        np.append(StateChangeMatrixI, StatechangeVector_I)

        StatechangeVector_S[specie] = changestateS
        np.append(StateChangeMatrixS, StatechangeVector_S)

    return StateChangeMatrixS, StateChangeMatrixI, MatrixPropensitiesS, MatrixPropensitiesI

StateMatrixS, StateMatrixI, PropS, PropI = GetMainMatrix(Metapop)

print('Matrice de changement detat S', StateMatrixS)

print('Matrice de proba S', PropS)