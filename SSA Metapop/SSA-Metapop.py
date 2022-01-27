# Algorithme de simulation stochastique pour les métapopulations
import os
import numpy
import numpy as np
import pandas as pd
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
print('La métapop', Metapop)

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
    list_si = []
    list_ii = []
    for pop in Metapop :
        list_si.append(pop[0])
        list_ii.append(pop[1])
    list_xi = list_si+list_ii
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

    #Specifier la taille des matrices de sortie à l'avance
    Nbevents = 7
    NbPops = len(Sarray)
    NbSpecies = len(Metapop) * len(Metapop[0])

    Nrow = Nbevents*NbPops
    Ncol = NbPops

    StateChangeMatrixS = np.zeros((Nrow, Ncol)) # Taille 70*10
    StateChangeMatrixI = np.zeros((Nrow, Ncol))
    MatrixPropensitiesS = np.zeros((Nrow, Ncol))
    MatrixPropensitiesI = np.zeros((Nrow, Ncol))

    print('Nombre de sites', NbPops)
    print('Nombre entités', NbSpecies)
    print('taille de matrice', StateChangeMatrixI.size)
    print('Taille de boucle', range(len(Sarray)) )
    for specie in range(len(Sarray)) :
        propensities_I = np.zeros(NbPops)
        StatechangeVector_I = np.zeros(NbPops)

        propensities_S = np.zeros(NbPops)
        StatechangeVector_S = np.zeros(NbPops)

        # On rempli deux array de taille Nbspecie
        #L'un contient les propensity de chaque évènement
        #L'autre contient les changements d'états

        S = Sarray[specie]
        I = Iarray[specie]

        # Pour que les bons index soient remplis
        Index = Nbevents * specie # Va de Nbevent en Nbevent en partant de zéro

        ##### On décline tous les évènements possibles pour remplir les changements d'états associés + propensity
        #Reproduction S
        prop = r* S
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestate = +1
        StatechangeVector_S[specie] = changestate
        StateChangeMatrixS[Index] = StatechangeVector_S

        Index += 1
        #Mort S
        prop = r * S * (S+I) /k
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestate = -1
        StatechangeVector_S[specie] = changestate
        StateChangeMatrixS[Index] = StatechangeVector_S

        Index += 1
        #Migration S
        prop = d * S
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestate = -1
        StatechangeVector_S[specie] = changestate
        StateChangeMatrixS[Index] = StatechangeVector_S
        Index += 1
        #Migration I

        prop = d * I
        propensities_I[specie] = prop
        MatrixPropensitiesI[Index] = propensities_I

        changestate = -1
        StatechangeVector_I[specie] = changestate
        StateChangeMatrixI[Index]= StatechangeVector_I

        Index += 1
        # Mort I
        prop = alpha * I
        propensities_I[specie] = prop
        MatrixPropensitiesI[Index] = propensities_I

        changestate = -1
        StatechangeVector_I[specie] = changestate
        StateChangeMatrixI[Index]= StatechangeVector_I

        Index += 1
        #Guérison I

        prop = alpha * I
        propensities_I[specie] = prop
        MatrixPropensitiesI[Index] = propensities_I

        changestateS = +1
        changestateI = -1
        StatechangeVector_I[specie] = changestateI
        StateChangeMatrixI[Index]= StatechangeVector_I

        StatechangeVector_S[specie] = changestateS
        StateChangeMatrixS[Index] = StatechangeVector_S

        Index += 1
        #Infection

        prop = beta * S * I
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestateS = -1
        changestateI = +1

        StatechangeVector_I[specie] = changestateI
        StateChangeMatrixI[Index]= StatechangeVector_I

        StatechangeVector_S[specie] = changestateS
        StateChangeMatrixS[Index] = StatechangeVector_S

    StateChangeMatrix = np.concatenate((StateChangeMatrixS, StateChangeMatrixI), axis = 1)
    MatrixPropensities = np.concatenate((MatrixPropensitiesS, MatrixPropensitiesI), axis = 1)

    return StateChangeMatrix, MatrixPropensities

StateMatrix, Props = GetMainMatrix(Metapop)
print('Taille matrice', StateMatrix.size)
print('Matrice de changement detat ', StateMatrix)

def GetCriticals(effectifs, propensities, statechange) : # Prends en entrée respectivement : Vecteur, Matrice, matrice
    ncrit = 11
    list_Indexcriticals = []
    dim_matrix = propensities.shape
    Critical_Matrix = np.zeros(dim_matrix)
    for index, effectif in enumerate(effectifs):
        if effectif < ncrit : list_Indexcriticals.append(index)
    print('BONJOUR',list_Indexcriticals) # OK

    for i in list_Indexcriticals : # Pour chaque pop critique
        print('VOICI LE I',i)
        for j in range(len(propensities)): # On regarde si aj > 0 ET vij < 0
            if propensities[j,i] > 0 and statechange[j,i] <0 : Critical_Matrix[j,i] =1
            else : pass
    print('REBONJOUR', Critical_Matrix)
    return Critical_Matrix
les_xi = Get_xi(Metapop)
Criticals = GetCriticals(les_xi, Props, StateMatrix)
print('Matrice de proba ', Props)
# On vérifie que la matrice à bien la tronche espérée
path = os.getcwd()
print(path)
Data = pd.DataFrame(data= Criticals)
print(Data)
Data.to_csv('Critics.csv')
