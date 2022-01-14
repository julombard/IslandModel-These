# Stochastic Simulation Algorithm input file
# S --> I --> S

# Reactions
ReproductionS: 
S > {2} S
r*S

MortS: 
S > $pool
r*S*(S+I)/k

Infection: 
S > I
beta * S * I

MigrationS: 
S > $pool
d * S

Guerison: 
I > S
gamma * I

MortI: 
I > $pool
alpha * I

MigrationI: # Migration infecté
I > $pool
d* I
# Fixed species

# InitVar
S = 999
I = 1

# InitPar
beta = 0.005  # Taux de contact infectieux
r = 1.5  # Taux de reproduction per capita
k = 1000  # Capacité de charge du site
d = 0.05  # Propension à la migration
gamma = 1.5  # Taux de guérison / clairance
alpha = 0.10  # Virulence

