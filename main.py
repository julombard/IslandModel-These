# Script Python pour modèle en Iles infinies avec dynamique locale

#tout plein d'imports qui font des choses plus ou moins funky
from copy import deepcopy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Essai de code de la dynamique locale pour un seul site
#xouxou

#Définition de variables de simulation
t = np.linspace(0,20, 500) # Donne le tmin, tmax, et le nombre de points qu'on souhaite entre les deux

#Définition des paramètres du modèle local
beta = 0.005 # Taux de contact infectieux
r = 1.5 # Taux de reproduction per capita
k = 1000 # Capacité de charge du site
d = 0.05 # Propension à la migration
gamma = 1.5 # Taux de guérison / clairance
alpha = 0.10 # Virulence

#Initialisation d'une population
dict_pop_init = { 'S' : 100, 'I' : 1}
Init = dict_pop_init['S'], dict_pop_init['I']


def Continuous_LocalDynamics(densite,t) :  # Fonction qui gère la dynamique locale des sites
    S, I = densite
    N = S + I

    dS = r * S * (1- N / k) - (beta * S * I) - (d * S) + (gamma * I)
    dI = (beta * S * I) - ((alpha + d + gamma) * I)
    return dS, dI

# On essaie de faire tourner le modèle

dynamics = odeint(Continuous_LocalDynamics, Init, t) # Renvoie un array avec les valeurs de S, I pour chaque point du temps
S, I = dynamics.T # Sépare l'array en deux listes distinctes


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

Nb_sites = 10
tsim = 10
t_petit = np.linspace(0,0.1,2)

# On initialise une pop par site
Pops = []
for i in range(Nb_sites) :
    Pops.append(Init)

# A ce moment on a besoin de plusieurs choses : Connaître le nombre de migrants
# Connaitre les valeurs de densités des pops pour chaque site au temps t+dt
# Appliquer la mortalité aux migrants
# Les redistribuer équiprobablement dans chaque site
for t in range(tsim):  # A chaque tour du modèle
    print('il court il court le furet',Pops)
    # Calcul du nombre de migrants par site
    Migrant_S = []  # Initialisation des listes qui contiendront les migrants de chaque site
    Migrant_I = []
    for pop in Pops:
        print('johnny', pop)
        Migrant_S.append(d * pop[0])
        print('Migrant S', Migrant_S)
        Migrant_I.append(d * pop[1])
        print('olalalalalala', type(pop[0]))
        dynamics = odeint(Continuous_LocalDynamics, pop,
                          t_petit)  # On résoud le système pour 1 étape avec un pas de temps tout petit, et on remplace les valeurs
        print('ehehehe', dynamics[1][0], type(dynamics[1][0]))
        pop[0] = int(dynamics[1][0])# Du coup on récupère la taille de pop au temps suivant, sans compter les immigrants de chaque population
        pop[1] = dynamics[1][1]
    #Calcul du nombre de migrants S qui vont atterir dans chaque patch
    print('coco pops',Pops)
    Nb_mig_s = sum(Migrant_S)
    Nb_mig_s_perpatch = Nb_mig_s / len(Pops)
    # Calcul du nombre de migrants I qui vont atterir dans chaque patch
    Nb_mig_i = sum(Migrant_I)
    Nb_mig_i_perpatch = Nb_mig_i / len(Pops)
    #Distribution de ces migrants là
    for pop in Pops:
        pop[0] = pop[0] + Nb_mig_s_perpatch
        pop[1] = pop[1] + Nb_mig_i_perpatch




