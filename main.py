# Script Python pour modèle en Iles infinies

#Un tas d'imports +/- rock'n'roll
from copy import deepcopy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import fonctions

########################################################################################################################
######################## CODE QUI SIMULE POUR LA DYNAMIQUE LOCAL D'UN SITE UNIQUE#######################################
########################################################################################################################


#Définition de variables de simulation
t = np.linspace(0,20, 500) # Donne le tmin, tmax, et le nombre de points qu'on souhaite entre les deux

#Initialisation d'une population
dict_pop_init = { 'S' : 100, 'I' : 0, 'Ms' : 0, 'Mi' : 0}
Init = dict_pop_init['S'], dict_pop_init['I'], dict_pop_init['Ms'], dict_pop_init['Mi']
#Init pour modèle à une pop
Init_local = dict_pop_init['S'], 1, dict_pop_init['Ms'], dict_pop_init['Mi']

# On fait tourner le modèle

dynamics = odeint(fonctions.Continuous_LocalDynamics, Init_local, t) # Renvoie un array avec les valeurs de S, I pour chaque point du temps
S, I , Ds, Ds = dynamics.T # Sépare l'array en n listes distinctes

# Pour tracer la figure de la dynamique du systeme
plt.figure(figsize=(20,10)) # taille de la fenetre graphique
plt.plot(t, S, label = 'susceptibles')
plt.plot(t,I,label = 'Infectes')
plt.xlabel('temps')
plt.ylabel('Nombre d individus')
plt.legend()
#plt.show() # Passé en commentraire car encombrant lors du run

########################################################################################################################
############################ MODELE VERSION ILES INFINIES (Fonctionnel)#################################################
########################################################################################################################

Nb_sites = 10
tsim = 100
t_petit = np.linspace(0,0.1,2) # On prends un point sur un pas très petit pour une résolution pas à pas
rho = 0.1 # Dispersal cost, utile car paramètres métapop exclusif

# On initialise une pop par site, dont une contient un infecté
Pops = [] # Liste de tuples, un tuple = une pop et contient les valeurs d'effectif
for i in range(Nb_sites-1) :
    Pops.append(Init)
Pops.append((100,1,0,0))

#On décrit le modèle
for t in range(tsim):  # A chaque tour du modèle
    print('Metapopulation',Pops)

    Migrant_S = []  # Initialisation des listes qui contiendront les migrants de chaque site
    Migrant_I = []
    for index, pop in enumerate(Pops): # On parcours chaque population de la métapop

        dynamics = odeint(fonctions.Continuous_LocalDynamics, pop,
                          t_petit)  # On résoud le système pour le pas t+dt

        Migrant_S.append(dynamics[1][2]) # On ajoute le nombre de migrants à la liste
        Migrant_I.append(dynamics[1][3])  # pop[x][y] est de type INT et on est contents

        new_pop = (dynamics[1][0],dynamics[1][1], 0, 0) # On introduit les nouvelles valeurs de S et I
        Pops[index] = new_pop # Et on met a jour la metapop
    #Quand on est passé sur toutes les pops, on calcule le nombre total de migrants
    Nb_mig_s_perpatch = sum(Migrant_S) / Nb_sites # Puis on divise par le nombre de patch pour avoir la qtité reçue pa chacun
    Nb_mig_i_perpatch = sum(Migrant_I) / Nb_sites
    for index, pop in enumerate(Pops): # Boucle de Redistribution des migrants

        new_pop = (pop[0]
                ,pop[1], Nb_mig_s_perpatch, Nb_mig_i_perpatch) # On crée une population avec les valeurs de S, I après dynamique locale, et les futurs immigrants
        final_pop = (new_pop[0] + new_pop[2], new_pop[1]+ new_pop[3], 0, 0 ) # On range les immigrants dans S, I et on repart avec 0 migrants "non assignés"
        Pops[index] = final_pop # Mise à jour de la metapop

