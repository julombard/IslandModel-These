# Stochastic Simulation Algorithm input file
# Doit être un fichier d'entrée de base de donnée (.psc)
# Essai pour un site
# --> S --> I

# Reactions
R1: # Reproduction des individus S
S > {2} S
r*S

R2 : #Mort des individus S
S > $pool
r*S*(S+I)/k

R3: #Infection
S > I
beta * S * I

R4: #Migration S
S > $pool
d * S

R5: #Guérison
I > S
gamma * I

R6: #Mort d'un infecté
I > $pool
alpha * I

R7 : # Migration infecté
I > $pool
d* I
# Fixed species

# Variable species
S = 999
I = 1

# Parameters
beta = 0.005  # Taux de contact infectieux
r = 1.5  # Taux de reproduction per capita
k = 1000  # Capacité de charge du site
d = 0.05  # Propension à la migration
gamma = 1.5  # Taux de guérison / clairance
alpha = 0.10  # Virulence

