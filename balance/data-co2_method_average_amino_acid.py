from utils.vectorOutput import generateOutputMatrixWithHeaders, writeOutputMatrix
from utils.specutils12 import GenerateFileList, ensure_dir
import pdb
from utils.balanceUtils import ImportStoichiometricMatrix, BalanceStoichiometricMatrix, \
PrintStoichiometry
from numpy import array


#directory = 'input/Amino Acid Synthesis Pathways/'
directory = 'input/Dataset S3 - Amino Acid Synthesis Pathways Draft 2/'


outputFileName = 'output/data-co2_method_average_amino_acid_2.csv'

ensure_dir(outputFileName)

startIndex = 2
endIndexOffset = 3
printIntermediates = True


reactantsToGet = ['ATP', 'NADH', 'Fdred', 'CO2', 'HCO3-', 'HCO2-']
co2FixationMethods = ['CBB', 'WL', 'RTCA', '3HP', '4HB', '3HP4HB', 'FORM']
# co2FixationMethods = ['4HB']


# Initialize dictionary for storing calculation results
co2FixationMethodsDict = {}
for method in co2FixationMethods:
	co2FixationMethodsDict[method] = {}
	for reactant in reactantsToGet:
		co2FixationMethodsDict[method][reactant] = []
	
co2FixationMethodsAverageDict = {}
for method in co2FixationMethods:
	co2FixationMethodsAverageDict[method] = {}
	for reactant in reactantsToGet:
		co2FixationMethodsAverageDict[method][reactant] = 0


# Calculate stoichiometry for each amino acid for each CO2 fixation method, and average across 
# amino acids
for method in co2FixationMethods:

	regText = r'.*' + method + r'\.csv'

	fileList = GenerateFileList(directory=directory, regex=regText, ignoreCase=True)

	for fileName in fileList:

		dictKey = fileName.split('-')[1]
		dictKey = dictKey.split('.')[0]
	
		print(dictKey)

		matrix = ImportStoichiometricMatrix(directory + '/' + fileName)
		endIndex = len(matrix[0]) - endIndexOffset
	
		[sMatrix, reactions, fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, reactants, ioStatus, \
		result] = \
		BalanceStoichiometricMatrix(matrix, startIndex, endIndex)

		PrintStoichiometry(cDotVectorOptNorm, reactants, ioStatus, \
		printIntermediates=printIntermediates)
		print()
		
		# Initialize storage for number of reactants 
		nReactantsDict = {}
		for reactant in reactantsToGet:
			nReactantsDict[reactant] = []
	
	
		for reactantToGet in reactantsToGet:
			# If the reactant to query isn't in the list of reactants for this matrix, record 0
			if reactantToGet not in reactants:
				nReactantsDict[reactantToGet].append(0)
				co2FixationMethodsDict[method][reactantToGet].append(0)
				
			# On the other hand, if it is, add the number of if to its array
			elif reactantToGet in reactants:
				reactantIndex = reactants.index(reactantToGet)
				nReactant = cDotVectorOptNorm[reactantIndex]
				co2FixationMethodsDict[method][reactantToGet].append(int(nReactant))
		
	# Finally, average across all amino acids
	for reactantToGet in reactantsToGet:
		co2FixationMethodsAverageDict[method][reactantToGet] \
		= average(co2FixationMethodsDict[method][reactantToGet])


	
# Compile together average reactant data
averageReactantsDict = {}
for reactant in reactantsToGet:
	averageReactantsDict[reactant] = []
	for method in co2FixationMethods:
		averageReactantsDict[reactant].append(co2FixationMethodsAverageDict[method][reactant])


# Output averaged data
headers = ['CO2 Fixation Method'] + reactantsToGet

vectorList = [co2FixationMethods]

for reactant in reactantsToGet:
	vectorList.append(array(averageReactantsDict[reactant], float)*-1)

 
oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
writeOutputMatrix(outputFileName, oMatrix)




