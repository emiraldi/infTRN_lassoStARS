#!/bin/bash
# tfTarget_GSEA.sh
# given a database of gene sets, test for enrichment of each set in the targets of each TF
# in the network. Positive and negative edges can be considered separately (see signOpts 
# below). GSEA results are further organized into tables (gene sets X TFs) of raw and 
# adjusted p-values, as well as other statistics, for downstream visualization (e.g., in 
# MATLAB)

geneSetAbbrev="geneSets_hMP_LPS"
geneSetFile="/Users/way5hl/Desktop/CCHMC/research/projects/scTRN/hMP/analysis/GRN/geneSets_GRN_hMP_LPS_DE_log2FC0p58_FDR10.txt"
backgroundGeneFile="/Users/way5hl/Desktop/CCHMC/research/projects/scTRN/hMP/inputs/gene_scRNA_Bryson_hMP_10X_QC_scater_N500M1_sct1000.txt"
	
# Network files (NOTE: for signed-analysis of TF targets, we work with TRNs where edges 
# have been filtered based on partial correlation, e.g., as in filter_TRNs_by_pcorr.sh
declare -a networkFiles=("/Users/way5hl/Desktop/CCHMC/research/projects/scTRN/hMP/analysis/GRN/outputs/TRN_scRNA_Bryson_hMP_10X_QC_scater_N500M1_Norm_tpm_Dist_CIDR_Sub_sct1k_Clust_NN_size20_adapt0_sigma1_Norm_tpm_pseudo1_log2_prior_ATAC_Bryson_MP_MACS_MOODS_p5_sGB10kb_bias50_sp_pCorrCut0p01.tsv" \
	"/Users/way5hl/Desktop/CCHMC/research/projects/scTRN/hMP/analysis/GRN/outputs/TRN_scRNA_Bryson_hMP_10X_QC_scater_N500M1_Norm_tpm_Dist_CIDR_Sub_sct1k_Clust_NN_size20_adapt0_sigma1_Norm_tpm_pseudo1_log2_prior_ATAC_Bryson_MP_MACS_MOODS_p5_sGB10kb_bias50_TFmRNA_sp_pCorrCut0p01.tsv" \
	"/Users/way5hl/Desktop/CCHMC/research/projects/scTRN/hMP/analysis/GRN/outputs/TRN_scRNA_Bryson_hMP_10X_QC_scater_N500M1_Norm_tpm_pseudo1_log2_prior_ATAC_Bryson_MP_MACS_MOODS_p5_sGB10kb_bias50_sp_pCorrCut0p01.tsv" \
	"/Users/way5hl/Desktop/CCHMC/research/projects/scTRN/hMP/analysis/GRN/outputs/TRN_scRNA_Bryson_hMP_10X_QC_scater_N500M1_Norm_tpm_pseudo1_log2_prior_ATAC_Bryson_MP_MACS_MOODS_p5_sGB10kb_bias50_TFmRNA_sp_pCorrCut0p01.tsv")

scriptHome="../infLassoStARS"

# Parameters for GSEA -- see tfTargetGSEA_sortPvals.py for definitions
interactionWeightCut="0"
# In initial GSEA, include an enrichment with a raw p-value <= to this cutoff
rawPcut=".1"
# Only test gene sets that overlap with our background gene set at this level:
minSetSize=5
# TF-gene interactions are signed, and so a TF's postive and negative targets can be 
# analyzed separately (signed or "dir" for directionally, signOpt = 1) or together 
# (unsigned, "abs" for absolute value of interaction, signOpt = 0). Signed analysis makes
# sense in terms of finding pathways that are turned on or off (note this could be 
# imperfect if some genes in the pathway are negative regulators). Unsigned analysis makes
# sense for GWAS-derived gene sets, where the genetic effect on gene function is often unknown.
declare -a signOpts=("1" "0")
# For putting results in table format for downstream visualization
fdrCut="1"
rawPcutTable=".1"

function log {
	echo $(date +%F_%T) $$ $BASHPID $1
	}

for networkFile in ${networkFiles[@]}
do
	outDir=$(dirname ${networkFile})/GSEA
	mkdir -p ${outDir}
	echo ${outDir}
	# -- directory for output
	fileBaseName=$(basename ${networkFile} _sp.tsv)
	# -- base for output file names 

	outDirBase=${outDir}/${fileBaseName}_${geneSetAbbrev}
	echo ${tfTargetOutDir}
	# output directory for motifEnrichemnts_TFtargets.py

	for signOpt in ${signOpts[@]}
	do						
		if [ "${signOpt}" -eq "1" ]; then
			signInf="dir"
			echo "signed GSEA of TF targets"
		elif [ "${signOpt}" -eq "0" ]; then
			signInf="abs"
			echo "unsigned GSEA of TF targets"
		else
			echo "signOpt must be either 0 or 1 (for unsigned/signed GSEA)"
			exit 64
		fi

		log "GSEA of ${netBase} WRT ${geneSetAbbrev}"
		echo "python ${scriptHome}/tfTargetGSEA_sortPvals.py ${networkFile} ${geneSetFile} ${backgroundGeneFile} ${rawPcut} ${signOpt} ${interactionWeightCut} ${minSetSize} ${outDirBase}"
		python ${scriptHome}/tfTargetGSEA_sortPvals.py ${networkFile} ${geneSetFile} \
			${backgroundGeneFile} ${rawPcut} ${signOpt} ${interactionWeightCut} \
			${minSetSize} ${outDirBase}					

		gseaFileName=${fileBaseName}_${geneSetAbbrev}_Praw${rawPcut/\./0p}_${signInf}_wCut0p0_minSet${minSetSize}
		ls ${outDir}/${gseaFileName}.txt
	
		log "Convert to table form, based on FDR cutoff -- for unsigned analysis"
		echo "python ${scriptHome}/makeGSEAmatrices_sepPosNeg_fxn.py ${outDir} ${fdrCut} ${outDir} ${gseaFileName} ${geneSetAbbrev}"
		python ${scriptHome}/makeGSEAmatrices_sepPosNeg_fxn.py \
			${outDir} \
			${fdrCut} \
			${outDir} \
			${gseaFileName} \
			${geneSetAbbrev}
			
		log "Convert to table form, based on raw p-value cutoff -- for unsigned analysis"	
		gseaFileName=${fileBaseName}_${geneSetAbbrev}_Praw${rawPcutTable/\./0p}_${signInf}_wCut0p0_minSet${minSetSize}
		ls ${outDir}/${gseaFileName}.txt
		echo "python ${scriptHome}/makeGSEAmatrices_sepPosNeg_prawCut_fxn.py \
			${outDir} \
			${rawPcutTable} \
			${outDir} \
			${gseaFileName} \
			${geneSetAbbrev}"

		python ${scriptHome}/makeGSEAmatrices_sepPosNeg_prawCut_fxn.py \
			${outDir} \
			${rawPcutTable} \
			${outDir} \
			${gseaFileName} \
			${geneSetAbbrev}

	done	
done

echo "Finished!"
