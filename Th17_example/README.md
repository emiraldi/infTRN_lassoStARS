# scTRN

This folder contains workflows for TRN analyses:
* Gene-set enrichment analysis (GSEA) of a TF's positive and repressed target genes
	* (optional) filter TF-gene interactions according to partial correlation: [filter_TRNs_by_pcorr.sh](filter_TRNs_by_pcorr.sh)
	* runGSEA on each TF's positive and repressed target genes: [tfTargets_GSEA_loop.sh](tfTargets_GSEA_loop.sh)
	* visualize enrichments with a heatmap: [visGSEAenrich_heatmaps_comb.m](visGSEAenrich_heatmaps_comb.m)
* Identify "core" TF regulators for a subset of conditions or celltypes in the gene expression dataset (see: [Miraldi et al., 2019](https://genome.cshlp.org/content/early/2019/02/19/gr.238253.118))
	* (optional) filter TF-gene interactions according to partial correlation: [filter_TRNs_by_pcorr.sh](filter_TRNs_by_pcorr.sh)
	* run GSEA to test whether a TF's postive or repressed targets are enriched in condition or celltype gene signatures: [tfTargets_GSEA.sh](tfTargets_GSEA.sh)
	* visualize "Top N" core TFs for the condition or celltype in a bargraph: [th17_coreAnalysis.m](th17_coreAnalysis.m)
* Visualize TF degree, positive and negative edges, and overlap with TF-gene interactions in the prior: [viz_TF_degree.m](viz_TF_degree.m)
