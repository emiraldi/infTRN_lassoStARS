#!/usr/bin/python

Usage = """
priorsTable2Sparse.py
convert a prior in table format (columns = regulators, rows = gene targets) into a sparse,
3-column format (regulator, row, interaction strength)
Usage: python ${scriptHome}/priorsTable2Sparse.py
INPUTS:
	priorTable -- tab-delimited file, columns = regulators, rows = gene targets, values = 
		interaction strengths 
	outFileName -- name for output file
	quantCut -- value s.t. interactions with |quantile rank| > quantCut will be included
		in resulting interaction file		
OUTPUTS:
	priorTable_sp -- sparse, 3-column output: 
		regulator, target, and interaction strength
"""	

import sys
import os
import re
import itertools
import cPickle as pickle
import numpy

## Debugging inputs

# currDir = '/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/DCs'
# #'/Users/emiraldi/erm/MariaP/Inferelator/input/Tr1_Th17_microarray_k17'
# # priors = ['encode_footprinting_priors_th17_allsamples.tsv',
# # 	'trrust_priors_th17.tsv',
# priors = ['Hh_DCs_ATAC_simple_deTF.tsv',
# 	'SFB_DCs_ATAC_simple_deTF.tsv']
# 
# #'DCs_ATAC_simple.tsv']
# 
# currCut = 0 
# 
# for prior in priors:
# 
# 	currPrior = currDir + '/' + prior	
# 	sys.argv = ['/Users/emiraldi/erm/MariaP/bin/priorsTable2Sparse.py',
# 		currPrior,
# 		currCut]

# /Users/emiraldi/erm/Shared/DCproject/RNAseq/DESeq2/DCs_exp40.txt

# 	sys.argv = ['/Users/emiraldi/erm/MariaP/bin/priorsTable2Sparse.py',
# 		'/Users/emiraldi/erm/MariaP/Inferelator/input/Tr1_Th17_microarray_k17/Tr1_Th17_noBatch_cvCut0025_k17.tsv',
# 		'/Users/emiraldi/erm/MariaP/Inferelator/input/Tr1_Th17_microarray_k17/Tr1_Th17_noBatch_cvCut0025_TFs_k17.tsv',
# 		currPrior,
# 		currCut]

	## END INPUTS

if len(sys.argv)<4:
	print Usage
	sys.exit(1)

priorGs=sys.argv[1]
quantCut=float(sys.argv[3])
outFileSparse=sys.argv[2]
# if quantCut:
# 	cutInf = "_cut" + str(quantCut)
# else:
# 	cutInf = ""
# 
# fileBase = os.path.splitext(priorGs)
# outFileSparse = fileBase[0] + cutInf + '_sp.tsv'
sparseOut = open(outFileSparse,'w')
sparseOut.write('regulator\ttarget\tweight\n')

# parse table of prior network
priorIn = open(priorGs,'r')
header = 1
totPairs = 0
totInts = 0
geneTot = 0

for line in priorIn:
	lineInf = list()
	lineInf = line.strip('\n').replace('"','').split('\t')
	if header:
		header = 0
		if len(lineInf[0]) < 1:
			# input table might have a tab in the header line
			lineInf = lineInf[1:]
# 			print str(len(lineInf)) + ' TFs possible in prior:'	
# 			for tf in sorted(lineInf):
# 				print tf.capitalize()
# 			print '\n'
# 		gompers
		keepInds = list()
		keepTfs = list()
		count = -1
		for tf in lineInf:
			count += 1
			tf = tf.capitalize()
			keepTfs.append(tf)
			keepInds.append(count)
		print str(len(keepTfs)) + ' TFs kept in prior.'
		intersPerTf = dict()
		for tf in keepTfs:
			intersPerTf[tf] = 0
	else:
		# we're dealing with interaction terms
		gene = lineInf[0].capitalize()
		geneTot += 1
		lineVals = lineInf[1:]
		sparseList = list()
		kind = -1
		for ind in keepInds:
			kind += 1
			currVal = float(lineVals[ind])

			if abs(currVal) > quantCut:
				sparseList.append(keepTfs[kind] + '\t' + gene + '\t' + str(currVal))
				tf = keepTfs[kind]
				intersPerTf[tf] += abs(currVal)
		if len(set(sparseList)) > 1:
			# if there are any non-zero values, print that line
			sparseOut.write('\n'.join(sparseList) + '\n')
			totPairs += 1
			totInts += len(set(sparseList)) - 1
sortedTfs = sorted(keepTfs)
print '\n'
print 'TF\t# Targets'
for tf in sortedTfs:
	print tf + '\t' + str(intersPerTf[tf])
print '\n'

print str(totInts) + ' / ' + str(geneTot*len(keepTfs)) + ' interactions.'
print str(totPairs) + ' / ' + str(geneTot) + ' genes with at least one TF interaction.'
print str(len(keepTfs)) + ' total TFs in prior.'

sparseOut.close()
print outFileSparse + ' generated.'
