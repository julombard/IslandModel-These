# Algorithme de simulation stochastique pour les métapopulations
import os
import numpy as np
import pandas as pd
import stochpy
import Fonctions
from copy import deepcopy
from itertools import islice

stochpy.SSA
# Implémentation du modèle
# Même principe que pour une population isolée MAIS
# On a N populations
# On tire les évènements pour UN individu parmis TOUTES les populations
# Si l'évènement est "migration" on doit définir un site d'arrivée, puis mettre sa population à jour

# Modified Poisson Tau leap algorithm from cao et. al. (2005)

# La sélection de tau est basée de "Efficient step size selection for tau leaping simulation", cao et al (2006)


Metapop = Fonctions.Set_metapop(100, 10) # Définir une métapop en moins de temps qu'il n'en faut pour définir la fonction !
# print('La métapop', Metapop)


StateMatrix, Props, Orders, aj = Fonctions.GetMainMatrix(Metapop)
print('lézaji', aj)


les_xi = Fonctions.Get_xi(Metapop)
Criticals = Fonctions.GetCriticals(les_xi, Props, StateMatrix)
mus, sigmas, NC_orders, a0nc = Fonctions.ComputeMuSigma(Props, StateMatrix, Criticals, Orders)

g_vect = Fonctions.GetEpsiloni(les_xi)
TauPrime = Fonctions.GetTauPrime(les_xi, mus, sigmas, g_vect)

# Regarder si Tau < 10* a0(x) ou a0(x) et la somme de toutes les propensities ####
a0x = a0 = np.sum(aj)
a0_crit = a0x-a0nc

TauCrit = 1/ a0x
#print('TAU CRITIQUE', TauCrit)
#print('TAU PRIME', TauPrime)
if TauPrime < TauCrit:
    print('Too bad, you are very mal tombé dans une partie non implemented') # Insérer direct Method Ici
else :
    TauPrimePrime = np.random.exponential(1/a0nc, 1)
    #cherche les reactions critiques
    Crit_index = np.where(Criticals == 1)  # Pour trouver les index des reactions critiques
    Critlist = list(Crit_index[0])  # Pour les traquer et dans une liste les lier

    if TauPrime < TauPrimePrime :
        Tau = TauPrime
        # On génère les kj  d'une loi de poisson pour toutes les réactions NC, si Critique -> kj =0
        kjs = []
        for i in range(len(aj)):
            if i in Critlist :
                kjs.append(0)
            else :
                mean = aj[i] * Tau
                kj = np.random.poisson(mean, 1)
                kjs.append(kj[0])
        print('Les KJS 1', kjs)
    else :
        Tau = TauPrimePrime
        kjs = []
        # Identifier la prochaine réaction critique
        jc = []
        for i in range(len(aj)):
            #print('me voilalaaalaalalal',aj[i])
            if i in Critlist :
                kjs.append(0) # On set à 0 par défaut
                jc.append(aj[i]/a0_crit) # On cal
            else :
                mean = aj[i] * Tau
                #print('LA MOYENNE', mean)
                kj=np.random.poisson(mean, 1)
                kjs.append(kj[0])
        # On définit l'évènement critique qui peut se déclencher une fois
        #On crée le vecteur des probas cumulées
        #print('Les probas extremes', jc)
        Index_event = np.random.choice(Critlist, 1, p=jc)
        print('EVENEMENT XTREM est le numéroooooo :',Index_event)
        kjs[Index_event[0]]= 1
        print('Les KJS 2', kjs)
    #print('Taille des Kj, normalement 70', len(kjs))
#print('TAU', Tau)


xi_tplusun = Fonctions.UpdateMetapop(les_xi, kjs, StateMatrix)


print('Pop à t+1', xi_tplusun)

# On vérifie que la matrice à bien la tronche espérée
path = os.getcwd()
print(path)
Data = pd.DataFrame(data= Criticals)
print(Data)
Data.to_csv('Critics.csv')
