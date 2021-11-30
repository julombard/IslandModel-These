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

#Initialisation d'une population
dict_pop_init = { 'S' : 100, 'I' : 1}
Init = dict_pop_init['S'], dict_pop_init['I']
print('eh coucou',type(Init))


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
plt.show() # Ne pas oublier, sinon le graphique ne s'affiche pas et on passe deux heures a chercher ce qu'on a fait de mal...

# OK POUR CA, faudra voir si on garde l'autorisation des demi-individus, peut être faire le truc en temps discret ?
#Tentative pour un modèle en îles multi-sites

