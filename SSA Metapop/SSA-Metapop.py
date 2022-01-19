# Algorithme de simulation stochastique pour les métapopulations
import numpy as np
import stochpy

stochpy.SSA
# Implémentation du modèle
# Même principe que pour une population isolée MAIS
# On a N populations
# On tire les évènements pour UN individu parmis TOUTES les populations
# Si l'évènement est "migration" on doit définir un site d'arrivée, puis mettre sa population à jour

# Modified Poisson Tau leap algorithm from cao et. al. (2005)


#Définition des paramètres du modèle
# InitPar
beta = 0.005  # Taux de contact infectieux
r = 1.5  # Taux de reproduction per capita
k = 1000  # Capacité de charge du site
d = 0.05  # Propension à la migration
gamma = 1.5  # Taux de guérison / clairance
alpha = 0.10  # Virulence



######## Step 1 : Fabriquer une métapop pour un état initial donné
def Set_metapop(taillepop, nbsites): #Initialise une métapopulation de nbsites sites contenant chacun taillepop individus, un site aléatoire contiendra une individu infecté
    Metapop=[]

    for i in range(nbsites):
        if i == 0: Metapop.append([taillepop-1,1])
        else :Metapop.append([taillepop,0])
    Metapoparray = np.array(Metapop)
    return Metapoparray
Metapop = Set_metapop(100, 10)
print('La métapop', Metapop)

def Tauleap() :
    return 0

def Compute_propensities(Metapop):
    #Metapop doit être un array de taille quelconque qui contient des listes de taille 2 [S,I]
    #Event est un array qui contient les taux pour chaque évènement

    Propensities = []

    for pop in Metapop :
        S = pop[0]
        I = pop[1]

        # On passe en revue tout ce qui peut se produire et on en détermine les chances
        # Reproduction S
        prop1 = r*S
        Propensities.append(prop1)
        # Mort S
        prop2 = r*S*(S+I)/k
        Propensities.append(prop2)
        # Infection
        prop3=beta * S * I
        Propensities.append(prop3)
        #MigrationS
        prop4=d * S
        Propensities.append(prop4)
        # Guerison:
        prop5=gamma * I
        Propensities.append(prop5)
        # Mort I
        prop6=alpha * I
        Propensities.append(prop6)
        # Migration I
        prop7=d * I
        Propensities.append(prop7)
    Sum_propensities = np.sum(Propensities)

    return Propensities, Sum_propensities

Propensities, SumProp = Compute_propensities(Metapop)
print('coucou',Propensities)
print('recoucou',SumProp)

def define_tau(Propensities, Metapop, Changestatevector) :
