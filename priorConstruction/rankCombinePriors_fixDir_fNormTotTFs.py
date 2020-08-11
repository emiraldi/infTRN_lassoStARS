#!/usr/bin/python

Usage = """
rankCombinePriors_fixDir_fNormTotTFs.py
Priors can take any set of values.  Here each prior is treated equally by normalizing 
edge weights to the Frobenius norm and multiplying by # of TFs in the prior before combining, 
then prior is re-scaled so that the  max edge length is "Total number of input networks"
rank combine prior networks by ADDING |edge weights| and using the sum of signed 
interactions to choose a sign
This might be more appropriate for genome-scale prior information sources (if you have information
about 10 TFs then the total edge weight contributed will be proportial to 10 TFs)
input matrices are in 3-column format (regulator, target, edge weight)
Provide:
	outDir -- folder for output matrices
	list of desired networks / priors:
		each network input is tuple (Network Name,[list of signed priors], [list of unsigned priors])	

"""	

import sys
import os
import re
import itertools
import cPickle as pickle
import numpy
import subprocess as sp
import math

sys.argv = ['/Users/emiraldi/erm/MariaP/bin/rankCombinePriors_fixDir_fNormTotTFs.py',
 '/Users/way5hl/Desktop/tmp/combined_priors',
 [('ATAC_Miraldi2019_TRRUST2hm',
   ['/Users/way5hl/Desktop/CCHMC/research/databank/priors/TRRUST2_Audrey/prior_TRRUST2_h_m_sp.tsv'],
   ['/Users/way5hl/Desktop/CCHMC/research/databank/ATAC_Miraldi2019/priors/prior_miraldi_ATAC_Th17_sp.tsv'])]]

# debugging
# sys.argv = ['/Users/emiraldi/erm/MariaP/bin/rankCombinePriors_fixDir_fNorm.py',
#  '/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/mm10_merged/ilc',
#  [('aC_lC_oC',
#    [],
#    ['/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/mm10_merged/ilc/osheaILC_cut4_bp10000_sATAC_p1Em5_huA_sp.tsv',
#     '/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/mm10_merged/ilc/amitILC_cut4_bp10000_sATAC_p1Em5_huA_sp.tsv',
#     '/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/mm10_merged/ilc/littmanILC_cut4_bp10000_sATAC_p1Em5_huA_sp.tsv'])]]
# sys.argv = ['/Users/emiraldi/erm/MariaP/bin/rankCombinePriors_fixDir_fNorm.py',
#  	'/Users/emiraldi/erm/MariaP/Inferelator/input/priorLists/priorList_3col_priors_fNorm',
# 	[('T_qAm_sAt17',
# 	['/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/mm10_merged/TRRUST_sp.tsv',
#     	'/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/mm10_merged/th17/invitroTh_1503_ATAC_max4_body_bp10000_pFDR25_FC0p58_gFDR25_FC0p58_allNets_sp.tsv'],
#    	['/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/mm10_merged/th17/Th17_48h_cut4_bp10000_sATAC_p1Em5_huA_sp.tsv'])]]

# each network input is tuple (Network Name,[list of signed priors], [list of unsigned priors])	

if len(sys.argv)<3:
	print Usage
	sys.exit(1)

outDir=sys.argv[1]
netPairs=sys.argv[2]

print outDir
print netPairs

for currNet in netPairs:
	netName = currNet[0]
	koFiles = currNet[1]
	nets2comb = currNet[2]
	
	totPriors = len(koFiles) + len(nets2comb)

# 	outFileTabUS = outDir + '/' + netName + '_us.tsv'	# output as a table UNSIGNED
# 	outFileSpUS = outDir + '/' + netName + '_us_sp.tsv'	# output as a sparse, 3-column matrix UNSIGNED
	outFileTab = outDir + '/' + netName + '.tsv'			# output as a table, SIGNED
	outFileSp = outDir + '/' + netName + '_sp.tsv'		# output as a sparse, 3-column matrix, SIGNED

	regs = list()
	targs = list()
	koInts = dict() # keys  key = reg\ttarg, item = weight
	inters = dict() # will store ALL interactions
	for koFile in koFiles:
		print koFile
		# calculate Frobenius norm
		# use awk to get squared Frobenius norm
		args = ['awk','{sum+=($3)^2} END {print sum}',koFile]
		p = sp.Popen(args, stdin = sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE)
		sqFrob = float(p.stdout.readline().strip('\n'))
		frobNorm = math.sqrt(sqFrob) # get Frobenius norm
		args = ['cut','-f1',koFile,'|','sort','|','uniq','|','wc','-l']
		cmd = " ".join(args)
		ps = sp.Popen(cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
		output = ps.communicate()[0]
		totTFs = float(output.strip('\n').replace(' ',''))-1
		frobNorm = frobNorm/totTFs # now there's some normalization for total number of TFs, I think this is okay for genome-scale assays		
		fin = open(koFile,'r')
		for line in itertools.islice(fin,1,None):
			lineInf = list()
			lineInf = line.strip('\n').split('\t')
			reg = lineInf[0]
			targ = lineInf[1]
			weight = float(lineInf[2])
			if reg in regs:
				try: 
					koInts[reg][targ] += weight/frobNorm
					inters[reg][targ] += abs(weight)/frobNorm
				except KeyError:
					koInts[reg][targ] = weight/frobNorm
					inters[reg][targ] = abs(weight)/frobNorm
			else:
				koInts[reg] = dict()
				koInts[reg][targ] = weight/frobNorm
				inters[reg] = dict()
				inters[reg][targ] = abs(weight)/frobNorm
				regs.append(reg)
			targs.append(targ)
		fin.close()

	# koIntsList = sorted(koInts.keys())

	# regs = list()
	# targs = list()
	# inters = dict() # keys = reg-targ, items = value associated with interaction
	# get regulators and targets from the positive interaction networks
	for inFile in nets2comb:
		print inFile
		# calculate Frobenius norm
		# use awk to get squared Frobenius norm
		args = ['awk','NR>1{sum+=($3)^2} END {print sum}',inFile]
		p = sp.Popen(args, stdin = sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE )
		sqFrob = float(p.stdout.readline().strip('\n'))
		frobNorm = math.sqrt(sqFrob) # get Frobenius norm	
		args = ['cut','-f1',inFile,'|','sort','|','uniq','|','wc','-l']
		cmd = " ".join(args)
		ps = sp.Popen(cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
		output = ps.communicate()[0]
		totTFs = float(output.strip('\n').replace(' ',''))-1
		frobNorm = frobNorm/totTFs # now there's some normalization for total number of TFs, I think this is okay for genome-scale assays					
		fin = open(inFile,'r')
		for line in itertools.islice(fin,1,None):
			lineInf = list()
			lineInf = line.strip('\n').split('\t')
			reg = lineInf[0]
			targ = lineInf[1]
			weight = abs(float(lineInf[2])) # this input is treated as unsigned
			if reg in regs:
				try:
					inters[reg][targ] += weight/frobNorm
				except KeyError:
					inters[reg][targ] = weight/frobNorm
			else:
				inters[reg] = dict()
				inters[reg][targ] = weight/frobNorm
				regs.append(reg)
			targs.append(targ)	
		fin.close()
		
	regs = sorted(list(set(regs)))
	targs = sorted(list(set(targs)))
	
	# Find the absolute max interaction value, and use it to divide
	maxVal = 0
	for targ in targs:
		for reg in regs:
			try:
				currVal = abs(inters[reg][targ])
				if currVal > maxVal:
					maxVal = currVal
			except KeyError:
				dog = 'no interaction'

	outTab = open(outFileTab,'w')
	outSp = open(outFileSp,'w')
# 	outTabUS = open(outFileTabUS,'w') 	# US = unsigned
# 	outSpUS = open(outFileSpUS,'w')		# US = unsigend

	outTab.write('\t' + '\t'.join(regs) + '\n')
	outSp.write('Regulator\tTarget\tWeight\n')
# 	outTabUS.write('\t' + '\t'.join(regs) + '\n')
# 	outSpUS.write('Regulator\tTarget\tWeight\n')

	intCount = 0
	for targ in targs:
# 		print targ
		outTab.write(targ)
# 		outTabUS.write(targ)
		for reg in regs:
			try:
				# is there an interaction
				intTot = float(totPriors)*float(inters[reg][targ])/maxVal
				intCount += 1
				intVal = str(intTot)
				# is it signed?
				try:
					intSig = numpy.sign(koInts[reg][targ])
					intSig = str(intSig*intTot)
					outTab.write('\t' + intSig)
					outSp.write(reg + '\t' + targ + '\t' + intSig + '\n')
				except KeyError:			
					outTab.write('\t' + intVal)
					outSp.write(reg + '\t' + targ + '\t' + intVal + '\n')
# 				outTabUS.write('\t' + intVal)
# 				outSpUS.write(reg + '\t' + targ + '\t' + intVal + '\n')
			except KeyError:
# 				outTabUS.write('\t0')
				outTab.write('\t0')
		outTab.write('\n')
# 		outTabUS.write('\n')

	outTab.close()
	outSp.close()
# 	outTabUS.close()
# 	outSpUS.close()

	print outFileTab
	print outFileSp
	print str(len(regs)) + ' regs'		
	print str(len(targs)) + ' targs'
	print str(intCount) + ' interactions'
	print ' ' 