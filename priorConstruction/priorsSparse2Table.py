#!/usr/bin/python

Usage = """
priorsSparse2Table.py
convert a prior in a sparse, 3-column format (regulator, row, interaction strength) to a 
table format (columns = regulators, rows = gene targets)
Usage: python ${scriptHome}/priorsSparse2Table.py priorSparseFile priorTableOutFile quantCut
INPUTS:
	priorSparseFile -- sparse, 3-column output: 
		regulator, target, and interaction strength
	priorTableOutFile -- name for the tabular prior
	quantCut -- value s.t. interactions with |quantile rank| > quantCut will be included
		in resulting interaction file		
OUTPUTS:
	priorTableOutFile -- tab-delimited file, columns = regulators, rows = gene targets, 
		values = interaction strengths 	
"""	

import sys
import os
import re
import itertools
import cPickle as pickle
import numpy

## Debugging inputs

# currDir = '/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/Dayanne260'
# # priors = ['encode_footprinting_priors_th17_allsamples.tsv',
# # 	'trrust_priors_th17.tsv',
# priors = ['th17PN.f1.InsertionPrior.unsigned.cut1.numbers_sp.tsv',
# 	'th17PN.f1_simpleATAC.unsigned.cut1.numbers_sp.tsv',
# 	'th17PN.mm9_hg19_InsertionPrior.unsigned.cut1.numbers_sp.tsv',
# 	'th17PN.mm9_hg19_simpleATAC.unsigned.cut1.numbers_sp.tsv']
# currCut = 0 

# currDir = '/Users/emiraldi/erm/Shared/Jason-Maria-Emily/ILCquantAtacNetworks/Genesets_FDR10/ILCs150523_qc_max4_body_bp10000_2comb/DCnets_pFDR10_FC1_rawp0001_hyg001_famShare'
# priors = ['allILCints_sp.tsv']
# currCut = 0
# 
# for prior in priors:
# 
# 	currPrior = currDir + '/' + prior	
# 	sys.argv = ['/Users/emiraldi/erm/MariaP/bin/priorsTable2Sparse.py',
# 		currPrior,
# 		currCut]
## END INPUTS

if len(sys.argv)<4:
	print Usage
	sys.exit(1)
priorSparseFile=sys.argv[1]
outFileTab=sys.argv[2]
quantCut=float(sys.argv[3])

# if quantCut:
# 	cutInf = "_cut" + str(quantCut)
# else:
# 	cutInf = ""

# fileBase = os.path.splitext(priorSparseFile)
# outFileTab = fileBase[0] + cutInf + '_tab.tsv'
# outFileTab = outFileTab.replace('_sp','')

# parse table of prior network
priorIn = open(priorSparseFile,'r')
header = 1
totPairs = 0
totInts = 0
geneTot = 0

inters = dict()
regList = []
targList = []
for line in priorIn:
	lineInf = list()
	lineInf = line.strip('\n').replace('"','').split('\t')
	if header:
		header = 0
	else:
		# we're dealing with interaction terms
		reg = lineInf[0]
		targ = lineInf[1]
		inter = lineInf[2]
		if reg in regList:
			inters[reg][targ] = inter
			targList.append(targ)
		else:
			inters[reg] = dict()
			inters[reg][targ] = inter
			regList.append(reg)
			targList.append(targ)
		targList = list(set(targList))

regList = sorted(regList)
targList = sorted(targList)

# output the table
tabOut = open(outFileTab,'w')				
tabOut.write('\t' + '\t'.join(regList) + '\n')

for targ in targList:
	tabOut.write(targ)
	for reg in regList:
		try:
			inter = inters[reg][targ]				
		except KeyError:
			inter = '0'
		tabOut.write( '\t' + inter)
	tabOut.write('\n')

tabOut.close()
print str(len(targList)) + ' targets in the prior.'
print str(len(regList)) + ' regulators in the prior.'
print outFileTab + ' generated.'
