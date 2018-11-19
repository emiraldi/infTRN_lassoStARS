#!/usr/bin/python

Usage = """
tfTargetGSEA_sortPvals.py
Perform GSEA on the targets of each TF in a network
Assumes the Gene names / IDs used in the network file are already matched to gene set file

USAGE:
	python ${scriptHome}/tfTargetGSEA_sortPvals.py networkFile geneSetFile backgroundGenes
		rawPcut signOpt interactionWeightCut minSetSize outDirBase
INPUTS:
	networkFile -- three column tab-separated format: TF, target genes, interaction weight (w/ header)
	geneSetFile -- three column space-seaparated format: set ID, set name, gene names separated with a "|" (w/o header)
	backgroundGenes -- a one-column list of genes to be considered in the background set,
		NOTE: The best background is the set of all genes with an annotation in database of gene sets
	rawPcut -- only output sets that have an adjusted p-value of enrichment below this cutoff
	signOpt -- 1 to take sign of TF-target interaction into account, 0 otherwise
	interactionWeightCut -- magnitude of interaction weight must be greater than this cutoff
	minSetSize -- minimum size of gene sets to be considered, e.g., for TF targets or
		gene sets, only those sets equivalent to or larger than this number will be tested
		for enrichment
	outDirBase -- basename of directory for the output
OUTPUTS:
	Each TF with enrichments will have it's own tabular output:
		col 0 = TF name
		col 1 = Gene set ID
		col 2 = Gene set name
		col 3 = # genes in set
		col 4 = # TF targets
		col 5 = # overlapping set and TF target genes
		col 6 = raw p-value
		col 7 = adjusted p-value
		col 8 = names of genes in the overlap, separated by spaces	
"""	

import os
import sys
import re
import itertools
import cPickle as pickle
import numpy
import errno
from scipy.stats import hypergeom
import glob

# sys.argv = ["/Users/emiraldi/erm/MariaP/bin/tfTargetGSEA.py",
# 	"/Users/emiraldi/erm/Shared/Jason-Maria-Emily/fiveILCquantAtacNetworks/ILCs150523_qc_max4_body_bp10000_pFDR10_FC1_gFDR10_FC1/ILCnets_pFDR10_FC1_gFDR10_FC1_rawp0001_hyg001_famShare/CCR6mILC3_SI_sp.tsv", #/Users/emiraldi/erm/Shared/Jason-Maria-Emily/ILCquantAtacNetworks/Genesets_FDR10/ILCs150523_qc_max4_body_bp10000_2comb/ILCnets_pFDR10_FC1_rawp0001_hyg001_famShare/ILC3_sp.tsv",
# 	"/Users/emiraldi/erm/MATLAB/Genesets/GWAS_160724_mu.txt",#"/Users/emiraldi/erm/MATLAB/Genesets/awongGOslim/go-mouse/gene_association.exp.slim.symbol",
# 	"/Users/emiraldi/erm/MATLAB/Genesets/awongGOslim/go-mouse/mm10_geneSymbolsUsed.txt", #"/Users/emiraldi/erm/tmpRfolder/ILC_RNA_1607/ILC_1607_genes.txt",
# 	".1",		# raw p-value cutoff
# 	"0",		# take sign of interaction into account
# 	"0",		# cutoff for TF-gene interaction weight strength
# 	"5",		# minimum set size 
# 	"/Users/emiraldi/erm/Shared/Jason-Maria-Emily/fiveILCquantAtacNetworks/ILCs150523_qc_max4_body_bp10000_pFDR10_FC1_gFDR10_FC1/ILCnets_pFDR10_FC1_gFDR10_FC1_rawp0001_hyg001_famShare/CCR6mILC3_SI_GWAS"]

# some parameters I set
minOverlapPrint = 2 # want to see at least 2 genes in enrichment

if len(sys.argv) < 8:
	print Usage
	sys.exit(1)

networkFile = sys.argv[1]
geneSetFile = sys.argv[2]
backgroundGenes = sys.argv[3]
rawPcut = float(sys.argv[4])
signOpt = int(sys.argv[5])
interactionWeightCut = float(sys.argv[6])
minSetSize = int(sys.argv[7])
outDirBase = sys.argv[8]

signName = 'abs'
if signOpt:
	signName = 'dir'
	
# make output directory
outDir = outDirBase + '_Praw' + str(rawPcut).replace('.','p') + '_' + signName + '_wCut' + \
	str(interactionWeightCut).replace('.','p') + '_minSet' + str(minSetSize)	
try:
    os.mkdir(outDir)
    print outDir
except OSError:
#     if e.errno != errno.EEXIST:
#         raise e
    pass

## Strategy
# 0. Get a list of background genes
# 1. Parse Gene Set file
#		geneSetDic : key = [gene set name, GO ID, geneSet size], items = genes
# 2. Parse Network file
#		tfSetDic : key = [TF name + sign of interaction], item = list of target genes
# 3. Calculate p-values,
#		pvalDic : key = [TF name + sign of interaction, gene set name, GO ID],
#			item = [overlap genes, pval]
#		pvalList : list of p-values, will be used for adj p-value estimates
#		raw2adjDic : key - raw p-value, item - adjusted p-value
# 4. Output

# 0. Get a list of background genes
backIn = open(backgroundGenes,'r')
backGenes = list()
for line in backIn:
	backGenes.append(line.strip('\n'))
backIn.close()
backGenes = set(backGenes)
totGenes = len(backGenes)

# 1. Parse Gene Set file
setIn = open(geneSetFile,'r')
geneSetDic = dict()
for line in setIn:
	# geneSetFile -- three column space-seaparated format: set ID, set name, gene names separated with a "|" (w/o header)
	lineInf = line.strip('\n').split(' ')
	setNames = lineInf[0:2]
# 	print setNames
	setGenes = set(lineInf[2].split('|'))
# 	print setGenes
	# limit the genes to those in the background set	
	setGenes = setGenes.intersection(backGenes)
# 	print setGenes
	geneSetSize = len(setGenes)
# 	gompers
	if geneSetSize >= minSetSize:
		geneSetDic[(setNames[0], setNames[1],str(geneSetSize))] = setGenes
setIn.close()
	
	
# 2. Parse Network file
netIn = open(networkFile,'r')
tfSetDic = dict()
for line in itertools.islice(netIn,1,None): # skip line 1, as it's a header
	# networkFile -- three column tab-separated format: TF, target genes, interaction weight (w/ header)
	lineInf = line.strip('\n').split('\t')
	tfName = lineInf[0]
	target = lineInf[1]
	if target in backGenes:	# limit to target genes in background set
		weight = float(lineInf[2])
		if abs(weight) >= interactionWeightCut:
			if not signOpt:
				weight = 0
			tfSetName = (tfName,str(int(numpy.sign(weight))).replace('-1','down').replace('1','up'))
			try:
				tfSetDic[tfSetName].append(target)
			except KeyError:
				tfSetDic[tfSetName] = [target,]
	else: 
		print target + " not in " + backgroundGenes
netIn.close()

# 3. Calculate p-values
pvalDic = dict()
pvalList = []
enrichDict = dict() # key = pval, item = list of sets with that enrichment
tfSetsUsed = dict() # track the tf Sets that were big enough
for tfSet in tfSetDic.keys():
	tfTargets = set(tfSetDic[tfSet])
	tfSetSize = len(tfTargets)
	if tfSetSize >= minSetSize:	# calculate enrichment
		tfName = tfSet[0]
		tfSetSign = tfSet[1]
		try:
			tfSetsUsed[tfName].append(tfSetSign)
		except KeyError:
			tfSetsUsed[tfName] = [tfSetSign,]
		for geneSet in geneSetDic.keys():
			setGenes = geneSetDic[geneSet]
			geneSetSize = len(setGenes)
			enTest = (tfSet, geneSet)
			overlap = sorted(tfTargets.intersection(setGenes))
			ovSize = len(overlap)
			enrich = 1 - hypergeom.cdf(max(ovSize-1,0),totGenes,tfSetSize,geneSetSize)
			pvalDic[enTest] = (overlap,enrich)
			pvalList.append(enrich)
			try: 
				enrichDict[enrich].append(enTest)
			except KeyError:
				enrichDict[enrich] = [enTest,]

# Calculate adjusted p-values, Benjamini Hochberg
totTests = len(pvalList)
sortedPs = sorted(pvalList)
currMult = totTests
raw2adjDic = dict()
for currP in sortedPs:
	currAdj = min(1,float(currP)*currMult)
	currMult += -1
	raw2adjDic[currP] = currAdj

# 4. Output TF enrichments
outFile = outDir + '.txt'
output = open(outFile,'w')
output.write('TF\tGeneSet_ID\tGeneSetName\t#Set genes\tDirection\t#TF targets\t# in overlap\tP_raw\tP_adj\tEnrichedGenes\n')
tfOutputDic = dict() # keys = tf, items = list of enrichments
for pval in sorted(set(sortedPs)):
	enTests = enrichDict[pval]
	for enTest in enTests:
		tfSet = enTest[0]
		tfName = tfSet[0]
		tfSetSign = tfSet[1]
		tfSetSize = str(len(tfSetDic[tfSet]))
		geneSet = enTest[1]
		overlap = pvalDic[enTest][0]
		pval = pvalDic[enTest][1]
		adjp = str(raw2adjDic[pval])
		if len(overlap) >= minOverlapPrint and pval < rawPcut:			
			output.write(tfName + '\t' + '\t'.join(geneSet) + '\t' + tfSetSign + '\t' + tfSetSize + '\t' + str(len(overlap)) + '\t' + str(pval) + '\t' + adjp + '\t' + ' '.join(overlap) + '\n')
			stuff2print = '\t'.join(geneSet) + '\t' + tfSetSign + '\t' + tfSetSize + '\t' + str(len(overlap)) + '\t' + str(pval) + '\t' + adjp + '\t' + ' '.join(overlap)
			try:
				tfOutputDic[tfName].append(stuff2print)
			except KeyError:
				tfOutputDic[tfName] = [stuff2print,]
output.close()
print outFile	

# output individual TF files -- (comment out)
# for tfName in tfOutputDic.keys():
# 	outFileTF = outDir + '/' + tfName + '.txt'
# 	outputTF = open(outFileTF,'w')
# 	outputTF.write('GeneSet_ID\tGeneSetName\t#Set genes\tDirection\t#TF targets\t# in overlap\tP_raw\tP_adj\tEnrichedGenes\n')
# 	outputTF.write('\n'.join(tfOutputDic[tfName]))
# 	outputTF.close()
# 	print outFileTF	
