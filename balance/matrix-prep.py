
from utils.balanceUtils import ImportReactionFile, ExportUniqueCompoundsWithIOStatus, \
GenerateIOStatusList
import os.path
from utils.specutils12 import ensure_dir

# dirName = 'input/co2_reactions/'
# exportDirName = 'output/co2_iostatus/'

# fileNames = [\
# '3HP4HB_reactions.txt', \
# '3HP_reactions.txt', \
# '4HB_reactions.txt', \
# 'CBB_reactions.txt', \
# 'FORM_reactions.txt', \
# 'WL_reactions.txt', \
# 'rTCA_reactions.txt' \
# ]

dirName = 'input/aa_reactions/'
exportDirName = 'output/aa_iostatus/'

fileNames = [\
# 'ALA_reactions.txt', \
# 'ARG_reactions.txt', \
# 'ASN_reactions.txt', \
# 'ASP_reactions.txt', \
# 'CYS_reactions.txt', \
# 'GLN_reactions.txt', \
# 'GLU_reactions.txt', \
# 'GLY_reactions.txt', \
# 'ILE_reactions.txt', \
# 'LEU_reactions.txt', \
# 'LYS_reactions.txt', \
'MET_reactions.txt', \
'PHE_reactions.txt', \
# 'PRO_reactions.txt', \
# 'SER_reactions.txt', \
# 'THR_reactions.txt', \
'TRP_reactions.txt', \
'TYR_reactions.txt' \
# 'VAL_reactions.txt' \
]

for fileName in fileNames:
	
	print(fileName)
	
	filePath = os.path.join(dirName, fileName)
	
	uniqueCompounds, reactions, sMatrixTKeyIndexed = \
	ImportReactionFile(filePath, reactionArrow='→')

	uniqueCompoundsIOStatus = GenerateIOStatusList(uniqueCompounds)

	baseName = os.path.basename(fileName)

	exportFileName = os.path.splitext(baseName)[0].split('_')[0] + '_iostatus.csv'
	exportFilePath = os.path.join(exportDirName, exportFileName)
	
	ensure_dir(exportFilePath)

	ExportUniqueCompoundsWithIOStatus(exportFilePath, uniqueCompoundsIOStatus)