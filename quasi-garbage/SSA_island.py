# Algorithme de simulation stochastique
# Fonctionnel sur Jupyter mais n'affiche rien sous Pycharm pour des raisons ma foi bien mystérieuses...
import stochpy
import os

# On récupère le répertoire de travail avec les fichiers d'entrée
path = os.getcwd()

model = stochpy.SSA(IsInteractive=False, method='tauleap', model_file='localSIS.psc')
model.Timesteps(10**5)
model.DoStochSim(IsTrackPropensities=True)
model.PlotSpeciesTimeSeries()
stochpy.plt.show()
stochpy.plt.savefig('timeseries.pdf')
model.PlotPropensitiesTimeSeries()
stochpy.plt.show()
stochpy.plt.savefig('propensities.pdf')

# Quelques données rigolotes
model.PrintSpeciesMeans()
model.PlotSpeciesDistributions()
stochpy.plt.show()