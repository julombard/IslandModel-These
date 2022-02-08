# Model parameters
beta = 0.005  # Infectious contact rate
r = 1.5  # Per capita growth rate
k = 1000  # Carrying capacity
d = 0.05  # Dispersal propensity
gamma = 1.5  # Parasite Clearance
alpha = 0.10  # Parasite Virulence
rho = 0.1  # Dispersal Cost
epsilon = 0.1  # Extinction rate

class Site():
    def __init__(self,effectifS,effectifI):
        self.effectifS = effectifS
        self.effectifI = effectifI
        #self.traitvalue = traitvalue #Affect a trait to each individual (vector) without being individual-based -later-
class Event():
    def __init__(self,name, propensity, Schange, Ichange, order):
        self.name = name
        self.S = 0
        self.I = 0
        self.formula = propensity
        self.propensity = eval(self.formula) #Convert string in maths instruction
        self.Ichange = eval(Ichange)
        self.Schange = eval(Schange)
        self.order = order

    def UpdatePropensity(self, S, I):
        self.S = S
        self.I = I
        self.propensity = eval(self.formula)
        return self.propensity