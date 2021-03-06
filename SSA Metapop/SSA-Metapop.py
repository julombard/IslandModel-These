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
# Paramètres de smulation
sim_time = 0
vectime = [0]
tmax = 40
nbsite = 2
taillepopinit = 100



Metapop = Fonctions.Set_metapop(taillepopinit, nbsite) # Définir une métapop en moins de temps qu'il n'en faut pour définir la fonction !
Outputs = [Metapop]
# print('La métapop', Metapop)
while sim_time < tmax :
    vectime.append(sim_time)
    StateMatrix, Props, Orders, aj = Fonctions.GetMainMatrix(Metapop)
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
        #print('Les KJS 1', kjs)
    else :
        Tau = TauPrimePrime
        kjs = []
        # Identifier la prochaine réaction critique
        jc = []
        for i in range(len(aj)):
            #print('me voilalaaalaalalal',aj[i])
            if i in Critlist :
                kjs.append(0) # On set à 0 par défaut
                jc.append(aj[i]/a0_crit)
            else :
                mean = aj[i] * Tau
                #print('LA MOYENNE', mean)
                kj=np.random.poisson(mean, 1)
                kjs.append(kj[0])
        # On définit l'évènement critique qui peut se déclencher une fois
        #On crée le vecteur des probas cumulées
        #print('Les probas extremes', jc)
        print('La critliste',Critlist)
        if Critlist  :
            Index_event = np.random.choice(Critlist, 1, p=jc)
            kjs[Index_event[0]] = 0 ### NE DOIT PAS ETRE A ZERO CETAIT POUR UN TEST FAIS GAFFE
        else : pass
        #print('EVENEMENT XTREM est le numéroooooo :',Index_event)

        #print('Les KJS 2', kjs)
    #print('Taille des Kj, normalement 70', len(kjs))
#print('TAU', Tau)


    Newmetapop = Fonctions.UpdateMetapop(les_xi, kjs, StateMatrix)
    Metapop = Newmetapop
    Outputs.append(Metapop)
    sim_time += Tau[0] # MaJ du temps
    print('Nous avons passé ', sim_time, ' temps')
print('VECTEUR TEMPS', vectime)
#Creation et mise en forme du DataFrame des sorties
data = pd.DataFrame(columns=['t'])
for i in range(nbsite):
    colname_s = 'S'+str(i)
    colname_i = 'I'+str(i)
    data[colname_s] = []
    data[colname_i] = []

print(data)

# On crée une liste vide par sous population, qu'on va mettre a jour avec les valeurs d'outputs
#On recupère ainsi toutes les séries dans l'ordre et dans des listes uniques
Si_values = []
Ii_values = []
for i in range(nbsite):
    Si_values.append([])
    Ii_values.append([])

for i in Outputs : # Pour chaque pas de temps
    #On crée une liste vide par sous pop
    for index, j in enumerate(i):
        #print('VOICI j', j)
        #print('VOICI j0', j[0])
        #print('INDEX', index)
        Si_values[index].append(j[0])
        Ii_values[index].append(j[1])
print('séries S' ,Si_values)
datalist= []
for i in range(len(Si_values)):
    datalist.append(Si_values[i])
    datalist.append(Ii_values[i])
print('courage fuyons', len(datalist))
# On remplit le dataframe avec les listes
for index,colname in enumerate(data.iteritems()):
    print('nom col',colname[0])
    print(index)
    if index == 0 : data[colname[0]] = vectime
    else : data[colname[0]] = datalist[index-1]

print(data)
data.to_csv('TimeSeries.csv')
