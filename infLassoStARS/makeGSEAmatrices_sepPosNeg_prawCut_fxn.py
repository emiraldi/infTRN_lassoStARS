#!/usr/bin/python

Usage = """
makeGSEAmatrices_sepPosNeg_fxn.py
Output Summary file from tfTargetGSEA_sortPvals.py and generates tables [GeneSets X TFs] 
for ease of visualization in Matlab
Positive and negative TF edges become separate files, also works for unsigned analyses
USAGE:
	python ${scriptHome}/makeGSEAmatrices_sepPosNeg_fxn.py
INPUTS:
	outDir1 -- output file name and location
	prawCut -- geneset AND TF must have one enrichment at this praw cutoff in order to be 
		included in the final data matrix
	inputDir -- input directory for GSEAs	
	gseaFile -- name of GSEA file as generated from, e.g., tfTargetGSEA_sortPvals.py	
	gseaSetsName -- name of GSEA sets
	prawMax -- if non-zero, then a sparse network output will be made (3 column format:
		TFs, gene, interaction weight = -log10(max(adjp,prawMax)) * sign (up or down) )
OUTPUTS:
	Tables (gene sets X TFs, where TFs are separated by + and - edges, if needed)
	outDir1_${gseaSetsName}_praw${prawCut}_pvals_${sign} -- raw pvalues [sets X TFs]
	outDir1_${gseaSetsName}_praw${prawCut}_adjps_${sign} -- adj. pvalues [sets X TFs]
	outDir1_${gseaSetsName}_praw${prawCut}_numGenes_${sign} -- # genes in overlap [sets X TFs]
	outDir1_${gseaSetsName}_praw${prawCut}_genes_${sign} -- space-separated genes in overlap [sets X TFs]
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

## Debugging
# prawCut = 1
# outDir1 = '/Users/emiraldi/erm/Shared/Maria-Emily/1712_Th17networks/sAaTh_bias50_combined'
# inputDir = '/Users/emiraldi/erm/Shared/Maria-Emily/1712_Th17networks/sAaTh_bias50_combined'
# gseaFile = 'GSEA/sAaTh_bias50_combined_cut01_aGOslim_Praw0p1_dir_wCut0p0_minSet5'
# gseaSetsName = 'aGOslimXXX'
# sys.argv = ["/Users/emiraldi/erm/MariaP/bin/makeGSEAmatrices_sepPosNeg.py",
# 	outDir1,
# 	prawCut,
# 	inputDir,
# 	gseaFile,
# 	gseaSetsName]

if len(sys.argv) < 5:
	print Usage
	sys.exit(1)

outDir1 = sys.argv[1]
prawCut = float(sys.argv[2])
inputDir = sys.argv[3]
gseaFile = sys.argv[4]
gseaSetsName = sys.argv[5]

outFolderBase = outDir1 + '/' + gseaFile 

# STRATEGY:
# construct dictionary: key = set, item = dict(), key = TF, item = rest of line (pvals, etc.)
# construct dictionary: key = set, item = source
# keep track of sets that reached prawCut in a list; do the same for TFs
enrichDic = dict()
setDic = dict()
sigTfs = dict() # -- keys = direction, item = TF sig
sigSets = dict() # -- keys = direction, item = sets sig
gseaTypes = list()


directionSymbolReplacements = {'down':'-','up':'+','0':'0'}
for direction in directionSymbolReplacements.keys():
	sigTfs[direction] = list()
	sigSets[direction] = list()
	enrichDic[direction] = dict()

# parse input file
resultsIn = open(inputDir + '/' + gseaFile + '.txt','r')
for line in itertools.islice(resultsIn,1,None):		
	tf,geneSet_ID,geneSetsName,setNum,direction,tfNum,numOverlap,praw,padj,geneNames = line.strip('\n').split('\t')
	tfName = tf # directionSymbolReplacements[direction] + tf
	if geneSetsName in enrichDic[direction].keys():	
		enrichDic[direction][geneSetsName][tfName] = line.strip('\n').split('\t')
	else:
		enrichDic[direction][geneSetsName] = dict()
		enrichDic[direction][geneSetsName][tfName] = dict()
		enrichDic[direction][geneSetsName][tfName] = line.strip('\n').split('\t')
	if float(praw) < prawCut:
		sigTfs[direction].append(tfName)
		sigSets[direction].append(geneSetsName)
resultsIn.close()

sortedSigTfs = dict()
sortedSigSets = dict()

# limit each list to unique TFs and sets
for direction in directionSymbolReplacements.keys():
	currTfs = list(set(sigTfs[direction]))
	sortedSigTfs[direction] = sorted(currTfs)
	currSets = list(set(sigSets[direction]))
	sortedSigSets[direction] = sorted(currSets)

# compose output tables
for direction in directionSymbolReplacements.keys():
	outDir1 = outFolderBase + '/' + gseaSetsName + '_praw' + str(int(100*prawCut)) + '_' + direction
	rawpOutFile = outDir1 + '_pval.txt'
	adjpOutFile = outDir1 + '_adjp.txt'
	ovTotsOutFile = outDir1 + '_ovTots.txt'
	ovGenesFile = outDir1 + '_ovGenes.txt'

	rawpOut = open(rawpOutFile,'w')
	adjpOut = open(adjpOutFile,'w')
	ovTotsOut = open(ovTotsOutFile,'w')
	ovGenesOut = open (ovGenesFile,'w')

	rawpOut.write('\t' + '\t'.join(sortedSigTfs[direction]) + '\n')
	adjpOut.write('\t' + '\t'.join(sortedSigTfs[direction]) + '\n')
	ovTotsOut.write('\t' + '\t'.join(sortedSigTfs[direction]) + '\n')
	ovGenesOut.write('\t' + '\t'.join(sortedSigTfs[direction]) + '\n')

	for geneSet in sortedSigSets[direction]:
		rawLine = geneSet
		adjLine = geneSet
		otLine = geneSet
		ogLine = geneSet
		for tf in sortedSigTfs[direction]:
			try:
				lineInf = enrichDic[direction][geneSet][tf]
				numOv = lineInf[6]
				praw = lineInf[7]
				padj = lineInf[8]
				genes = lineInf[9]
			except KeyError:
				numOv = '0'
				praw = '1'
				padj = '1'
				genes = ''
			rawLine += '\t' + praw
			adjLine += '\t' + padj
			otLine += '\t' + numOv
			ogLine += '\t' + genes
		rawpOut.write(rawLine + '\n')
		adjpOut.write(adjLine + '\n')
		ovTotsOut.write(otLine + '\n')
		ovGenesOut.write(ogLine + '\n')
	rawpOut.close()
	adjpOut.close()
	ovTotsOut.close()
	ovGenesOut.close()	

	print str(len(sortedSigSets[direction])) + ' enriched gene sets'
	print str(len(sortedSigTfs[direction])) + ' enriched TFs'
	print rawpOutFile # = outDir1 + '_pval.txt'
	print adjpOutFile # = outDir1 + '_adjp.txt'
	print ovTotsOutFile #= outDir1 + '_ovTots.txt'
	print ovGenesFile #= outDir1 + '_ovGenes.txt'

print "Finished"