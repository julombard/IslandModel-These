
def Continuous_LocalDynamics(densite,t) :  # Fonction qui gère la dynamique locale des sites
    # Définition des paramètres du modèle local
    beta = 0.005  # Taux de contact infectieux
    r = 1.5  # Taux de reproduction per capita
    k = 1000  # Capacité de charge du site
    d = 0.05  # Propension à la migration
    gamma = 1.5  # Taux de guérison / clairance
    alpha = 0.10  # Virulence
    rho = 0.1  # Dispersal cost


    S, I , Ds, Di = densite
    N = S + I

    dS = r * S * (1- N / k) - (beta * S * I) - (d * S) + (gamma * I)
    dI = (beta * S * I) - ((alpha + d + gamma) * I)
    dDs = d * S
    dDi = d * I
    return dS, dI, dDs, dDi
