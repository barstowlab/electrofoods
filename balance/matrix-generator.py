import re
from numpy import unique, int8
from utils.balanceUtils import SolveFluxBalanceEquation, ConvertIndexedSMatrix, \
ImportReactionFile, GenerateIndexedSMatrixT, PrintStoichiometry, GenerateMergedIOStatusList, \
ImportIOStatus
import os.path
from utils.vectorOutput import generateOutputMatrixWithHeaders, writeOutputMatrix



aaIOStatusDirName = 'input/aa_iostatus'
co2IOStatusDirName = 'input/co2_iostatus'
aaReactionsDirName = 'input/aa_reactions'
co2ReactionsDirName = 'input/co2_reactions'


aaIOStatusFileNames = [\
'ALA_iostatus.csv', \
'ARG_iostatus.csv', \
'ASN_iostatus.csv', \
'ASP_iostatus.csv', \
'CYS_iostatus.csv', \
'GLN_iostatus.csv', \
'GLU_iostatus.csv', \
'GLY_iostatus.csv', \
'ILE_iostatus.csv', \
'LEU_iostatus.csv', \
'LYS_iostatus.csv', \
'MET_iostatus.csv', \
'PHE_iostatus.csv', \
'PRO_iostatus.csv', \
'SER_iostatus.csv', \
'THR_iostatus.csv', \
'TRP_iostatus.csv', \
'TYR_iostatus.csv', \
'VAL_iostatus.csv' \
]

aaReactionsFileNames = [\
'ALA_reactions.txt', \
'ARG_reactions.txt', \
'ASN_reactions.txt', \
'ASP_reactions.txt', \
'CYS_reactions.txt', \
'GLN_reactions.txt', \
'GLU_reactions.txt', \
'GLY_reactions.txt', \
'ILE_reactions.txt', \
'LEU_reactions.txt', \
'LYS_reactions.txt', \
'MET_reactions.txt', \
'PHE_reactions.txt', \
'PRO_reactions.txt', \
'SER_reactions.txt', \
'THR_reactions.txt', \
'TRP_reactions.txt', \
'TYR_reactions.txt', \
'VAL_reactions.txt' \
]

co2ReactionsFileNames = [\
'3HP4HB_reactions.txt', \
'3HP_reactions.txt', \
'4HB_reactions.txt', \
'CBB_reactions.txt', \
'FORM_reactions.txt', \
'WL_reactions.txt', \
'rTCA_reactions.txt' \
]

co2IOStatusFileNames = [\
'3HP4HB_iostatus.csv', \
'3HP_iostatus.csv', \
'4HB_iostatus.csv', \
'CBB_iostatus.csv', \
'FORM_iostatus.csv', \
'WL_iostatus.csv', \
'rTCA_iostatus.csv' \
]


reactantsToGet = ['ATP', 'NADH', 'Fdred', 'Sulfate', 'N2', 'CO2', 'HCO3-', 'HCO2-', 'CO']
reactantOutputFilePath = 'output/amino_acids/stable-nReactants.csv'


# ------------------------------------------------------------------------------------------------ #
# Generate a set of dictionaries for amino acid io status

aa_ioStatusDict = {}
for fileName in aaIOStatusFileNames:
	
	prefix = os.path.splitext(fileName)[0].split('_')[0]
	ioStatusAAFilePath = os.path.join(aaIOStatusDirName, fileName)
	iostatusAA = ImportIOStatus(ioStatusAAFilePath)
	
	aa_ioStatusDict[prefix] = iostatusAA
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Generate a set of dictionaries for co2 io status

co2_ioStatusDict = {}
for fileName in co2IOStatusFileNames:	

	prefix = os.path.splitext(fileName)[0].split('_')[0]
	ioStatusCO2FilePath = os.path.join(co2IOStatusDirName, fileName)
	iostatusCO2 = ImportIOStatus(ioStatusCO2FilePath)
	
	co2_ioStatusDict[prefix] = iostatusCO2
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Generate a set of dictionaries for aa reactions

aa_reactionDict = {}
for fileName in aaReactionsFileNames:	

	prefix = os.path.splitext(fileName)[0].split('_')[0]
	aa_reactionFilePath = os.path.join(aaReactionsDirName, fileName)

	compounds_AA, reactions_AA, sMatrixTKeyIndexed_AA = \
	ImportReactionFile(aa_reactionFilePath, reactionArrow='→')
	
	aa_reactionDict[prefix] = [compounds_AA, reactions_AA, sMatrixTKeyIndexed_AA]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Generate a set of dictionaries for co2 reactions

co2_reactionDict = {}
for fileName in co2ReactionsFileNames:	

	prefix = os.path.splitext(fileName)[0].split('_')[0]
	co2_reactionFilePath = os.path.join(co2ReactionsDirName, fileName)

	compounds_co2, reactions_co2, sMatrixTKeyIndexed_co2 = \
	ImportReactionFile(co2_reactionFilePath, reactionArrow='→')
	
	co2_reactionDict[prefix] = [compounds_co2, reactions_co2, sMatrixTKeyIndexed_co2]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Loop over each amino acid and generate a stoichiometric matrix for each CO2 fixation method

aaKeys = aa_reactionDict.keys()
co2Keys = co2_reactionDict.keys()

# Generate a scenario dict

scenarioDict = {}
scenarioReactantsDict = {}
scenarioTargetDict = {}

for aaKey in aaKeys:

	# Import the unique amino acid synthesis reactions here
	[compounds_AA, reactions_AA, sMatrixTKeyIndexed_AA] = aa_reactionDict[aaKey]
	iostatusAA = aa_ioStatusDict[aaKey]

	# Loop over each carbon-fixation reactions here
	for co2Key in co2Keys:

		scenarioKey = aaKey + '_' + co2Key
		scenarioReactantsDict[scenarioKey] = {}
		scenarioTargetDict[scenarioKey] = {}
		
		print(scenarioKey)
		
		[compounds_CO2, reactions_CO2, sMatrixTKeyIndexed_CO2] = co2_reactionDict[co2Key]
		iostatusCO2 = co2_ioStatusDict[co2Key]

		reactions_full = list(reactions_CO2) + list(reactions_AA)
		compounds_full = list(compounds_CO2) + list(compounds_AA)
		compoundsUnique_full = list(unique(compounds_full))

		sMatrixTKeyIndexed = GenerateIndexedSMatrixT(compoundsUnique_full, reactions_full, \
		reactionArrow='→')

		mergedIOStatus, ioStatusCO2Dict, ioStatusAADict = \
		GenerateMergedIOStatusList(compoundsUnique_full, iostatusCO2, iostatusAA)

		# Convert the stoichiometric matrix to a non-key indexed matrix
		sMatrixT = ConvertIndexedSMatrix(sMatrixTKeyIndexed, compoundsUnique_full)
		sMatrix = sMatrixT.transpose()

		[fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, result] = \
		SolveFluxBalanceEquation(sMatrix, reactions_full, compoundsUnique_full, mergedIOStatus)
	
		scenarioDict[scenarioKey] = [fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, result, \
		compoundsUnique_full, reactions_full, sMatrixTKeyIndexed]
		
		for reactantToGet in reactantsToGet:
			# If the reactant to query isn't in the list of reactants for this matrix, record 0
			if reactantToGet not in compoundsUnique_full:
				scenarioReactantsDict[scenarioKey][reactantToGet] = 0
			
			# On the other hand, if it is, add the number of if to its array
			elif reactantToGet in compoundsUnique_full:
				reactantIndex = compoundsUnique_full.index(reactantToGet)
				nReactant = cDotVectorOptNorm[reactantIndex]
				scenarioReactantsDict[scenarioKey][reactantToGet] = nReactant
	
		i = 0
		while i < len(compoundsUnique_full):
			ioStatusReactant = mergedIOStatus[i]
				
			if ioStatusReactant == 'Target':	
				scenarioTargetDict[scenarioKey][compoundsUnique_full[i]] = cDotVectorOptNorm[i]
				
			i += 1
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Output the results of stoichiometry calculation as a table

scenarioKeys = list(scenarioReactantsDict.keys())
nReactantsDict = {}

for reactantToGet in reactantsToGet:
	nReactantsDict[reactantToGet] = []
	for key in scenarioKeys:
		nReactantsDict[reactantToGet].append(scenarioReactantsDict[key][reactantToGet])

targetsArray = []
nTargetsArray = []

for key in scenarioKeys:
	targetMoleculeKeys = scenarioTargetDict[key].keys()
	targetMolecule = list(targetMoleculeKeys)[0]
	nTargetMolecule = scenarioTargetDict[key][targetMolecule]
	targetsArray.append(targetMolecule)
	nTargetsArray.append(nTargetMolecule)
	

headers = ['Scenario'] + reactantsToGet + ['Target Molecule', 'Target']
vectorList = [scenarioKeys]

for reactantToGet in reactantsToGet:
	vectorList.append(array(nReactantsDict[reactantToGet], float)*-1)

vectorList.append(targetsArray)
vectorList.append(array(nTargetsArray, float))

oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
writeOutputMatrix(reactantOutputFilePath, oMatrix)
# ------------------------------------------------------------------------------------------------ #


		
		

