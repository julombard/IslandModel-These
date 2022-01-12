# Algorithme de simulation stochastique
import stochpy
import os

# On récupère le répertoire de travail avec les fichiers d'entrée

stochpy.SSA(method='Direct', model_file='localSIS.psc')
model = stochpy.SSA()
model.DoStochSim(IsTrackPropensities=True)
model.PlotSpeciesTimeSeries()
model.PlotPropensitiesTimeSeries()