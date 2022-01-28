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
#print('La métapop', Metapop)

def Tauleap4Metapop() : # Will be used someday, or maybe not
    return 0

def Get_xi(Metapop) : # Get density values of each subpoppulation of each site
    list_si = []
    list_ii = []
    for pop in Metapop :
        list_si.append(pop[0])
        list_ii.append(pop[1])
    list_xi = list_si+list_ii
    xi = np.array(list_xi)
    return xi

def GetMainMatrix(Metapop) : # Fonction qui récupère la matrice des propensities et celle qui contient les vecteurs de changement d'états

    #Séparation des S et I par convenance
    Slist = []
    Ilist = []

    for pop in Metapop :
        Slist.append(pop[0])
        Ilist.append(pop[1])
    Sarray = np.array(Slist)
    Iarray= np.array(Ilist)

    #Specifier la taille des matrices de sortie à l'avance
    Nbevents = 7
    NbPops = len(Sarray)
    #NbSpecies = len(Metapop) * len(Metapop[0])

    Nrow = Nbevents*NbPops
    Ncol = NbPops

    StateChangeMatrixS = np.zeros((Nrow, Ncol)) # Taille 70*10
    StateChangeMatrixI = np.zeros((Nrow, Ncol))
    MatrixPropensitiesS = np.zeros((Nrow, Ncol))
    MatrixPropensitiesI = np.zeros((Nrow, Ncol))

    # POur récupérer les ordres de réactions (utile après)
    Orders = []

    #print('Nombre de sites', NbPops)
    #print('Nombre entités', NbSpecies)
    #print('taille de matrice', StateChangeMatrixI.size)
    #print('Taille de boucle', range(len(Sarray)) )
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
        order = 1
        Orders.append(order)

        prop = r* S
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestate = +1
        StatechangeVector_S[specie] = changestate
        StateChangeMatrixS[Index] = StatechangeVector_S

        Index += 1

        #Mort S
        order = 3
        Orders.append(order)
        prop = r * S * (S+I) /k
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestate = -1
        StatechangeVector_S[specie] = changestate
        StateChangeMatrixS[Index] = StatechangeVector_S

        Index += 1

        #Migration S
        order = 1
        Orders.append(order)

        prop = d * S
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestate = -1
        StatechangeVector_S[specie] = changestate
        StateChangeMatrixS[Index] = StatechangeVector_S
        Index += 1
        #Migration I
        order = 1
        Orders.append(order)

        prop = d * I
        propensities_I[specie] = prop
        MatrixPropensitiesI[Index] = propensities_I

        changestate = -1
        StatechangeVector_I[specie] = changestate
        StateChangeMatrixI[Index]= StatechangeVector_I

        Index += 1
        # Mort I
        order = 1
        Orders.append(order)
        prop = alpha * I
        propensities_I[specie] = prop
        MatrixPropensitiesI[Index] = propensities_I

        changestate = -1
        StatechangeVector_I[specie] = changestate
        StateChangeMatrixI[Index]= StatechangeVector_I

        Index += 1
        #Guérison I
        order = 1
        Orders.append(order)

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
        order = 2
        Orders.append(order)

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
    Reaction_orders = np.array(Orders)
    return StateChangeMatrix, MatrixPropensities, Orders

StateMatrix, Props, Orders = GetMainMatrix(Metapop)
print('Les ordres de réaction', Orders)
print('Nombre ordres', len(Orders))

def GetCriticals(effectifs, propensities, statechange) : # Prends en entrée respectivement : Vecteur, Matrice, matrice
    ncrit = 11
    list_Indexcriticals = []
    dim_matrix = propensities.shape
    Critical_Matrix = np.zeros(dim_matrix)
    for index, effectif in enumerate(effectifs):
        if effectif < ncrit : list_Indexcriticals.append(index)
    for i in list_Indexcriticals : # Pour chaque pop critique

        for j in range(len(propensities)): # On regarde si aj > 0 ET vij < 0
            if propensities[j,i] > 0 and statechange[j,i] <0 : Critical_Matrix[j,i] =1
            else : pass
    vect_crit = np.sum(Critical_Matrix, axis= 1) # Axe 1 pour sommer sur les lignes !!!! Ne te gourres plus, ca va bien maintenant !
    return vect_crit # Renvoie un vecteur de booléens

def ComputeMuSigma(propensities, statechange, criticals, reaction_orders): # Matrice, matrice matrice vecteur

    ### Chercher les réctions critiques
    Crit_indexes = np.where(criticals==1) # Pour trouver les index des réactions critiques
    Crit_list = list(Crit_indexes[0]) # Pour les traquer et dans une liste les lier
    print('LA LISTE', Crit_list) # Bon ca fonctionne

    # On calcule la matrice produit des aj(x) * vij (et on triera ce qu'on garde après)
    Matrice_produit_mu = abs(statechange * propensities)
    Matrice_produit_sigma = statechange**2 * propensities

    # On retire de ces matrices les lignes correspondant à des réactions critiques
    Mat_mu = np.delete(Matrice_produit_mu, Crit_list, axis = 0)
    Mat_sigma = np.delete(Matrice_produit_sigma, Crit_list, axis = 0) # EN passant une liste en 2e argument, on s'évite une boucle

    # On fait la somme des colonnes des matrices pour récupérer les vecteurs des mus et sigmas
    Vect_mu = np.sum(Mat_mu, axis = 0)
    Vect_Sigma = np.sum(Mat_sigma, axis = 0)

    # On récupère les ordres des réactions non critiques, car c'est utile juste après
    Ncrit_Orders = np.delete(reaction_orders, Crit_list)

    return Vect_mu, Vect_Sigma, Ncrit_Orders

def GetEpsiloni(xi) : # Vecteur des ordre des réactions non-critiques
    epsilon = 0.03 # Valeur donnée dans l'article
    Nbtypes = int(len(xi)/2) # Bouuuuh c'est laid ça
    g_vector = []

    # Pour les i le HOR (higghest order reaction) est toujours 2 -> beta *s *i
    # Pour les s le HOR est touours 3 (mortalité densité dpdt)
    # Definir la sélection des gi en fonction des ordres de réaction
    # Les xi sont ordonnés de S1 à Sn et de I1 à In
    # Donc HOR = 3 pour x1 à xn et =2 pour xn+1 à xm où m est le nombre total d'entités

    for index, specie in enumerate(xi) :
        if index < Nbtypes : # Si on parcours les S
            x = specie
            g = 3/2 * (1/(x-1))
            g_vector.append(g)
        elif index >= Nbtypes : # Si on parcours les I
            g = 2
            g_vector.append(g)

    epsilon_i =[i / epsilon for i in g_vector]
    return epsilon_i

def GetTauPrime(xi, mu, sigma, epsilon): # Que des vecteurs
    # Prendre tous les max entre epsilon*xi et 1
    Vect_xiepsi = xi * epsilon
    print('COUCUUUUUU', Vect_xiepsi[0])
    Upper_term = []
    Tau_candidates = []

    for i in Vect_xiepsi :
        Upper_term.append(max(i, 1))
    # Prendre tous les min entre les deux tau candidats pour chaque espèce
    for i in range(len(xi)) :
        Candidate_mu= Upper_term[i]/ mu[i]
        Candidate_sigma=Upper_term[i]**2/ sigma[i]
        Tau_candidates.append(min([Candidate_mu, Candidate_sigma]))
    # Prendre le min entre tous les candidats qui ont passé les qualifs
    TauPrime = min(Tau_candidates)
    return TauPrime

les_xi = Get_xi(Metapop)
Criticals = GetCriticals(les_xi, Props, StateMatrix)
mus, sigmas, NC_orders = ComputeMuSigma(Props, StateMatrix, Criticals, Orders)
g_vect = GetEpsiloni(les_xi)
TauPrime = GetTauPrime(les_xi, mus, sigmas, g_vect)

#### Regarder si Tau < 10* a0(x) ou a0(x) et la somme de toutes les propensities ####

#### Si non, on fait gillespie de base sur 100 itérations ###

#### Si oui on va chercher TauPrimePrime #####

#### On compare TauPrime avec TauPrimePrime et on avise ####

# On vérifie que la matrice à bien la tronche espérée
path = os.getcwd()
print(path)
Data = pd.DataFrame(data= Criticals)
print(Data)
Data.to_csv('Critics.csv')
