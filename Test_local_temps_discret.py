#tout plein d'imports qui font des choses plus ou moins funky
from copy import deepcopy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Définition de variables de simulation
tsim = 20 # Donne le tmin, tmax, et le nombre de points qu'on souhaite entre les deux

#Définition des paramètres du modèle local
beta = 0.06 # Taux de contact infectieux
r = 1.5 # Taux de reproduction per capita
k = 1000 # Capacité de charge du site
d = 0.05 # Propension à la migration
gamma = 1.5 # Taux de guérison / clairance
alpha = 0.10 # Virulence

#Initialisation d'une population
dict_pop_init = { 'S' : 99, 'I' : 1}

def Discrete_LocalDynamics (dict_pop) :
    Nb_S = dict_pop['S']
    Nb_I = dict_pop['I']

    S_tplus1 = Nb_S  + r * Nb_S * (1- (Nb_S+Nb_I)/ k) - beta * Nb_S * Nb_I - d * Nb_S + gamma * Nb_I
    I_tplus1 = beta * Nb_S * Nb_I - d * Nb_I - alpha * Nb_I - gamma * Nb_I

    #Gestion des cas négatifs
    if S_tplus1 < 0 : S_tplus1 = 0
    if I_tplus1 < 0 : I_tplus1 = 0

    Pop_tplus1 = {'S' : S_tplus1, 'I' : I_tplus1}

    return Pop_tplus1

Pop = dict_pop_init
for i in range(tsim) :
    Pop = Discrete_LocalDynamics(Pop)
    print(Pop)