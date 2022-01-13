#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# Last updated by Buz Barstow on 2022-1-13
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Bargraph, Generate_EfficienciesDict_Keys_Sorted_by_Efficiency, \
Export_Efficiency_Bargraph

from rewiredcarbon.utils import ensure_dir




scenarioTableFileName = 'input/fig-cbb_n2_to_glycine_faradaic.csv'
outputFilenameEff = 'output/fig-cbb_h2_butanol_faradaic/fig-cbb_h2_butanol_faradaic.csv'
outputFilenameFuelMassEff = 'output/fig-cbb_h2_butanol_faradaic/fig-cbb_h2_butanol_faradaic_mass.csv'

ensure_dir(outputFilenameEff)
ensure_dir(outputFilenameFuelMassEff)

scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict)


keysArray = \
Generate_EfficienciesDict_Keys_Sorted_by_Efficiency(efficienciesDict, 'effTotalElectricalToFuel')


keysArray = list(efficienciesDict.keys())

# ------------------------------------------------------------------------------------------------ #

# Get data out
electronsPerFuelH2 = efficienciesDict['GLY_CBB_H2']['electronsPerFuel']
electronsDownForFuelH2 = efficienciesDict['GLY_CBB_H2']['electronsDownATPforFuelInH']
h2PerFuel = electronsPerFuelH2/2
o2UsedPerFuelH2 = electronsDownForFuelH2/4
o2MadePerFuelH2 = h2PerFuel/2
o2UsedToMadeRatioH2 = o2UsedPerFuelH2/o2MadePerFuelH2

electronsPerFuelEEU = efficienciesDict['GLY_CBB_EET']['electronsPerFuel']
electronsDownForFuelEEU = efficienciesDict['GLY_CBB_EET']['electronsDownForFuel']
o2UsedPerFuelEEU = electronsDownForFuelEEU/4
o2MadePerFuelEEU = electronsPerFuelEEU/4
o2UsedToMadeRatioEEU = o2UsedPerFuelEEU/o2MadePerFuelEEU





print('#----------------------------------------------------------#')
print('Electron and O2 use in H2-mediated EMP')
print('Total electrons per fuel molecule: ' + str(electronsPerFuelH2))
print('Total electrons downhill per fuel molecule: ' + str(electronsDownForFuelH2))
print('Total H2 per fuel molecule: ' + str(h2PerFuel))
print('Total O2 consumed per fuel molecule: ' + str(o2UsedPerFuelH2))
print('O2 generated to O2 used ratio O2: ' + str(o2UsedToMadeRatioH2))


print('#----------------------------------------------------------#')
print('')
print('#----------------------------------------------------------#')
print('Electron and O2 use in EEU-mediated EMP')
print('Total electrons per fuel molecule: ' + str(electronsPerFuelEEU))
print('Total electrons downhill per fuel molecule: ' + str(electronsDownForFuelEEU))
print('Total O2 made per fuel molecule: ' + str(o2MadePerFuelEEU))
print('Total O2 consumed per fuel molecule: ' + str(o2UsedPerFuelEEU))
print('O2 generated to O2 used ratio O2: ' + str(o2UsedToMadeRatioEEU))


print('#----------------------------------------------------------#')
