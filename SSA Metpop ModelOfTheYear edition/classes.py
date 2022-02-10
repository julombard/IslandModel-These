#Model parameters
beta = 0.005  # Infectious contact rate
r = 1.5  # Per capita growth rate
k = 1000  # Carrying capacity
d = 0.05  # Dispersal propensity
gamma = 1.5  # Parasite Clearance
alpha = 0.10  # Parasite Virulence
rho = 0.1  # Dispersal Cost
epsilon = 0.1  # Extinction rate

class Site(): #Site object containing (non explicit) individuals
    def __init__(self,effectifS,effectifI, *args): #First try with arg way to implement feature unsure
        self.effectifS = effectifS #S density
        self.effectifI = effectifI #I density
        #self.traitvalue = traitvalue Affect a trait to each individual (vector) without being individual-based
        #self.pos = pos (tuple) : for future improvements, position of the site on the network grid, as matrix coordinates
        #self.neighbor = [] : for future, maybe including neighbors as an attribute (so it's computed only once)

class Event():
    def __init__(self,name, propensity, Schange, Ichange, order):
        self.name = name # Event name in letter and not in memory address, handful to identify what's happening
        self.S = 0 #Has to take density values to understand the maths
        self.I = 0
        self.formula = propensity # The unique formule (str) given by model construction
        self.propensity = eval(self.formula)#Convert string in maths instruction, very useful to externalise model building
        self.Ichange = eval(Ichange) # State Change due to event, Typically -1, 0, 1 except for extinctions
        self.Schange = eval(Schange)
        self.order = order # Reaction order, not really useful but we never know

    def UpdatePropensity(self, S, I): # Class method to compute propensities without creating new objects
        self.S = S
        self.I = I
        self.propensity = eval(self.formula) # Changes propensity values while keeping formula untouched
        return self.propensity