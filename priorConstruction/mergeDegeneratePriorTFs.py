#!/usr/bin/python

Usage = """
mergeDegeneratePriorTFs.py
Given degeneracy of TF motifs, many priors based on ATAC-seq data contain TFs with identical
target genes and interaction strengths.  This function creates meta-TFs (whose name will 
be a compound of individual TFs in the merge (e.g., TF1_TF2...)).  Note: because mergers
often contained many TF names, we cut off at the first two, and then user can track down
the rest of the group in the merger file below  
USAGE:
	python ${scriptHome}/mergeDegeneratePriorTFs.py networkFile outFileBase
INPUTS:
	networkFile -- three column tab-separated format: TF, target genes, interaction weight (w/ header)
	outFileBase -- output file base name
OUTPUTS:
	0. new networkFile, where meta-TF names replace individual TFs with identical targets
	netOutFile = outFileBase + '_merged_sp.tsv'
	1. a target overlap table (TF X TF) w/ # of of overlapping targets
	overlapsOutFile = outFileBase + '_overlaps.txt'
	2. a table, containing number of targets per TF	
	targetTotalsFile = outFileBase + '_targetTotals.txt'
	3. a table, column 1 = abbreviated merged TF name, column 2 = all TFs included in merger
	mergedTfsFile = outFileBased + '_mergedTfs.txt'

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

## Debugging Inputs:
# sys.argv = ["/Users/emiraldi/erm/MariaP/bin/mergeDegeneratePriorTFs.py",
# 	"/Users/emiraldi/erm/Shared/Jason-Maria-Emily/fiveILCquantAtacNetworks/ILCs150523_qc_max4_body_bp10000_pFDR10_FC1_gFDR10_FC1/ILCnets_pFDR10_FC1_gFDR10_FC1_rawp0001_hyg001_famShare/CCR6pILC3_SI_sp.tsv", #/Users/emiraldi/erm/Shared/Jason-Maria-Emily/ILCquantAtacNetworks/Genesets_FDR10/ILCs150523_qc_max4_body_bp10000_2comb/ILCnets_pFDR10_FC1_rawp0001_hyg001_famShare/ILC3_sp.tsv",
# 	"/Users/emiraldi/erm/Shared/Jason-Maria-Emily/fiveILCquantAtacNetworks/ILCs150523_qc_max4_body_bp10000_pFDR10_FC1_gFDR10_FC1/ILCnets_pFDR10_FC1_gFDR10_FC1_rawp0001_hyg001_famShare/CCR6pILC3_SI"]
# sys.argv = ['/Users/emiraldi/erm/MariaP/bin/mergeDegeneratePriorTFs.py',
# 	'/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/mm10_merged/ilc/amitILC_cut4_bp10000_sATAC_p1Em4_huA_sp.tsv',
# 	'/Users/emiraldi/Desktop/test']


if len(sys.argv) < 3:
	print Usage
	sys.exit(1)

networkFile = sys.argv[1]
outFileBase = sys.argv[2]

## Strategy
# 1. Parse Network file
#		tfTargDic : key = TF name, item = list of target genes
# 2. Determine overlaps among each TF pair.  Make a dictionary of TFs to merge:
#		key = TF, item = merge group name
# 3. Output

# 1. Parse Network file
netIn = open(networkFile,'r')
tfTargDic = dict()
for line in itertools.islice(netIn,1,None): # skip line 1, as it's a header
	# networkFile -- three column tab-separated format: TF, target genes, interaction weight (w/ header)
	lineInf = line.strip('\n').split('\t')
	tfName = lineInf[0]
	target = lineInf[1]
	weight = lineInf[2]
	targetWeight = target + '*&$%||' + weight
	try:
		tfTargDic[tfName].append(targetWeight)
	except KeyError:
		tfTargDic[tfName] = [targetWeight,]
netIn.close()

# 2. Determine TF overlaps
tfNames = sorted(tfTargDic.keys())
totTfs = len(tfNames)
tfMergers = dict() # key = TF, item = list of TFs that should be merged
overlaps = dict() # key = TF1-TF2, item = overlap
tfTargNums = dict() # key = TF, item = number of targets
for tf1ind in range(0,totTfs):
	tf1 = tfNames[tf1ind]
	tf1targets = set(tfTargDic[tf1])
	numTf1Targs = len(tf1targets)
	tfTargNums[tf1] = numTf1Targs
	for tf2ind in range(tf1ind+1,totTfs):
		tf2 = tfNames[tf2ind]
		tf2targets = set(tfTargDic[tf2])
		numTf2Targs = len(tf2targets)	
		tfTargNums[tf2] = numTf2Targs
		overlap = tf2targets.intersection(tf1targets)
		overlapSize = len(overlap)
		overlaps[tf1 + '_' + tf2] = overlapSize
		overlaps[tf2 + '_' + tf1] = overlapSize
		# check to see if the targets/interaction signs are identical
		if numTf1Targs == numTf2Targs and numTf1Targs == overlapSize:
			if tf1 in tfMergers.keys():
				tfMergers[tf1].append(tf2)
			else:
				tfMergers[tf1] = [tf1,tf2]
			if tf2 in tfMergers.keys():
				tfMergers[tf2].append(tf1)
			else: 
				tfMergers[tf2] = [tf1,tf2]

# 3. Output
allMergedTfs = set(tfMergers.keys())
usedMergedTfs = list()	# keep tracks of used TFs, so we don't output mergers twice		
printedTfs = list()
overlapsToPrint = dict() # key = printTfName, item = overlaps with other TFs
netOutFile = outFileBase + '_merged_sp.tsv'
overlapsOutFile = outFileBase + '_overlaps.txt'
targetTotalsFile = outFileBase + '_targetTotals.txt'
mergedTfsFile = outFileBase + '_mergedTfs.txt'
netOut = open(netOutFile,'w')
netOut.write('Regulator\tTarget\tWeight\n')
overlapsOut = open(overlapsOutFile,'w')
targetTotalsOut = open(targetTotalsFile,'w')
mergedTfsOut = open(mergedTfsFile,'w')

for tf in tfNames:
	# if a TF has been merged and the compound TF has not already been printed,
	# add change its name to the compound TF name (and print it)
	printIt = 0
	if tf in allMergedTfs and tf not in usedMergedTfs:
		mergedTfs = tfMergers[tf]
		for mTf in mergedTfs:
			usedMergedTfs.append(mTf)
		if len(mergedTfs) > 2:	
			tfPrint = '_'.join(mergedTfs[0:2]) + '...'
		else:
			tfPrint = '_'.join(mergedTfs)
		mergedTfsOut.write(tfPrint + '\t' + ', '.join(mergedTfs) + '\n')
		printIt = 1
	elif tf not in allMergedTfs:
		tfPrint = tf
		printIt = 1
	if printIt:
		targetTotalsOut.write(tfPrint + '\t' + str(tfTargNums[tf]) + '\n')
		for targ in tfTargDic[tf]:
			netOut.write(tfPrint + '\t' + '\t'.join(targ.split('*&$%||')) + '\n')
		overlapsToPrint[tfPrint] = []
		for tf2 in printedTfs:
			# if either name has been merged, use only the first
			indTf1 = tfPrint.split('_')
			indTf2 = tf2.split('_')
			overlapsToPrint[tfPrint].append(str(overlaps[indTf1[0] + '_' + indTf2[0]]))
		printedTfs.append(tfPrint)
mergedTfsOut.close()
targetTotalsOut.close()
netOut.close()
print mergedTfsFile
print targetTotalsFile
print netOutFile
# print overlaps
overlapsOut.write('\t' + '\t'.join(printedTfs) + '\n')
for tf in printedTfs:
	printStuff = overlapsToPrint[tf]
	totSims = len(printStuff)
	zerosToPad = len(printedTfs) - totSims - 1
	if len(printStuff) > 0:
		overlapsOut.write(tf + '\t' + '\t'.join(printStuff) + '\t' + str(tfTargNums[tf.split('_')[0]]) + '\t0'*zerosToPad + '\n')
	else:
		overlapsOut.write(tf + '\t' + str(tfTargNums[tf.split('_')[0]]) + '\t0'*zerosToPad + '\n')
overlapsOut.close()
print overlapsOutFile		
