import stochpy

smod = stochpy.SSA()
smod.DoStochSim(IsTrackPropensities=True)
smod.PlotSpeciesTimeSeries()
smod.PlotPropensitiesTimeSeries()