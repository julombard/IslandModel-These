# Script Python pour modèle en Iles infinies avec dynamique locale

#tout plein d'imports qui font des choses plus ou moins funky
from copy import deepcopy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Essai de code de la dynamique locale pour un seul site


#Définition de variables de simulation
t = np.linspace(0,20, 500) # Donne le tmin, tmax, et le nombre de points qu'on souhaite entre les deux

#Définition des paramètres du modèle local
beta = 0.005 # Taux de contact infectieux
r = 1.5 # Taux de reproduction per capita
k = 1000 # Capacité de charge du site
d = 0.05 # Propension à la migration
gamma = 1.5 # Taux de guérison / clairance
alpha = 0.10 # Virulence
rho = 0.1 # Dispersal cost

#Initialisation d'une population
dict_pop_init = { 'S' : 100, 'I' : 0, 'Ms' : 0, 'Mi' : 0}
Init = dict_pop_init['S'], dict_pop_init['I'], dict_pop_init['Ms'], dict_pop_init['Mi']
Init_local = dict_pop_init['S'], 1, dict_pop_init['Ms'], dict_pop_init['Mi']


def Continuous_LocalDynamics(densite,t) :  # Fonction qui gère la dynamique locale des sites
    S, I , Ds, Di = densite
    N = S + I

    dS = r * S * (1- N / k) - (beta * S * I) - (d * S) + (gamma * I)
    dI = (beta * S * I) - ((alpha + d + gamma) * I)
    dDs = d * S
    dDi = d * I
    return dS, dI, dDs, dDi

# On essaie de faire tourner le modèle

dynamics = odeint(Continuous_LocalDynamics, Init_local, t) # Renvoie un array avec les valeurs de S, I pour chaque point du temps

S, I , Ds, Ds = dynamics.T # Sépare l'array en deux listes distinctes


# Pour tracer la figure de la dynamique du systeme
plt.figure(figsize=(20,10)) # taille de la fenetre graphique
plt.plot(t, S, label = 'susceptibles')
plt.plot(t,I,label = 'Infectes')
plt.xlabel('temps')
plt.ylabel('Nombre d individus')
plt.legend()
#plt.show() # Ne pas oublier, sinon le graphique ne s'affiche pas et on passe deux heures a chercher ce qu'on a fait de mal...

# OK POUR CA, faudra voir si on garde l'autorisation des demi-individus, peut être faire le truc en temps discret ?

#Tentative pour un modèle en îles multi-sites

def Continuous_Globaldynamics(densite,t) :
    S, I, Ds, Di = densite

    dS = Ds
    dI = Di

    return S, I, dS, dI


Nb_sites = 10
tsim = 100
t_petit = np.linspace(0,1,2) # On prends un point sur un pas très petit pour une résolution pas à pas

# On initialise une pop par site
Pops = []
for i in range(Nb_sites-1) :
    Pops.append(Init)
Pops.append((100,1,0,0))
# A ce moment on a besoin de plusieurs choses : Connaître le nombre de migrants
# Connaitre les valeurs de densités des pops pour chaque site au temps t+dt
# Appliquer la mortalité aux migrants
# Les redistribuer équiprobablement dans chaque site
for t in range(tsim):  # A chaque tour du modèle
    print('Metapopulation',Pops)
    # Calcul du nombre de migrants par site
    Migrant_S = []  # Initialisation des listes qui contiendront les migrants de chaque site
    Migrant_I = []
    for index, pop in enumerate(Pops):

        dynamics = odeint(Continuous_LocalDynamics, pop,
                          t_petit)  # On résoud le système pour 1 étape avec un pas de temps tout petit on remplace les valeurs

        # On remplace les valeurs, du coup on récupère la taille de pop au temps suivant, sans compter les immigrants de chaque population
        Migrant_S.append(dynamics[1][2]) # Extraction du nombre de migrants
        #print('Migrant S',Migrant_S)
        #print('Migrant', Migrant_I)
        #print('migrants S', Migrant_S)
        Migrant_I.append(dynamics[1][3])  # pop[x][y] est de type INT et on est contents
        new_pop = (dynamics[1][0],dynamics[1][1], 0, 0) # On remplace les valeurs et on repart sans migrants "non assignés"
        Pops[index] = new_pop
        #print('Population après dyn locale', new_pop)
    Nb_mig_s_perpatch = sum(Migrant_S) / Nb_sites
    Nb_mig_i_perpatch = sum(Migrant_I) / Nb_sites
    #print('Migrant S par site',Nb_mig_s_perpatch)
    for index, pop in enumerate(Pops): # Redistribution des migrants
        #print('population avant dyn globale', pop)
        new_pop = (pop[0]
                ,pop[1], Nb_mig_s_perpatch, Nb_mig_i_perpatch) # On remplace les valeurs et on repart sans migrants "non assignés"
        #print('pop juuuuuuust avant lODE', new_pop )
        final_pop = (new_pop[0] + new_pop[2], new_pop[1]+ new_pop[3], 0, 0 )
        #print('population après dyn globale', final_pop)
        Pops[index] = final_pop