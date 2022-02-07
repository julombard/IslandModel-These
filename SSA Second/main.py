class Site():
    def __init__(self, effectifS,effectifI, traitvalue):
        self.effectifS = effectifS
        self.effectifI = effectifI
        self.traitvalue = traitvalue #A terme on en fera une liste avec une valeur par individu pour avoir l'evolution
class Event():
    def __init__(self, rate, Schange, Ichange):
        self.rate = rate
        self.Ichange = Ichange
        self.Schange = Schange


#Definition des parametres principaux
