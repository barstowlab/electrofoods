#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# fig-co2fixation.py
# Calculates area and thickness of thin film in hydrogen-transport by diffusion scenario
# 
# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2019-10-25
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Bargraph, Generate_EfficienciesDict_Keys_Sorted_by_Efficiency, \
Export_Efficiency_Bargraph

from rewiredcarbon.utils import ensure_dir




scenarioTableFileName = 'input/fig-cbb_n2_to_amino_acids.csv'
outputFilenameEff = 'output/fig-cbb-n2_to_amino_acids/fig-cbb_n2_to_amino_acids_eff.csv'
outputFilenameFuelMassEff = 'output/fig-cbb-n2_to_amino_acids/fig-cbb_n2_to_amino_acids_fuel_mass.csv'

ensure_dir(outputFilenameEff)
ensure_dir(outputFilenameFuelMassEff)

scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict)


# keysArray = \
# Generate_EfficienciesDict_Keys_Sorted_by_Efficiency(efficienciesDict, 'effTotalElectricalToFuel')


keysArray = list(efficienciesDict.keys())

Plot_Efficiency_Bargraph(efficienciesDict, 'effTotalElectricalToFuel', \
'effTotalElectricalToFuel_lowerError', 'effTotalElectricalToFuel_upperError', keysToPlot=keysArray)


Plot_Efficiency_Bargraph(efficienciesDict, 'effTotalElectricalFuelMassEfficiency', \
'effTotalElectricalFuelMassEfficiency_lowerError', \
'effTotalElectricalFuelMassEfficiency_upperError', keysToPlot=keysArray)


Export_Efficiency_Bargraph(outputFilenameEff, efficienciesDict, scenarioDict, \
'effTotalElectricalToFuel', 'effTotalElectricalToFuel_lowerError', \
'effTotalElectricalToFuel_upperError', keysToPlot=keysArray)

Export_Efficiency_Bargraph(outputFilenameFuelMassEff, efficienciesDict, scenarioDict, \
'effTotalElectricalFuelMassEfficiency', 'effTotalElectricalFuelMassEfficiency_lowerError', \
'effTotalElectricalFuelMassEfficiency_upperError', keysToPlot=keysArray)
