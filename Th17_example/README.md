# Th17_example

This folder contains workflows for TRN analyses:
* Basic TRN inference pipeline: [filter_TRNs_by_pcorr.sh](filter_TRNs_by_pcorr.sh)
	* TRN construction
	* precision-recall analysis with a gold standard of TF-gene interactions
* TF-TF module analysis: Discovery of TFs that co-regulate gene pathways: [example_Th17_tfTfModules.m](example_Th17_tfTfModules.m)
* Out-of-sample gene expression prediction: [example_workflow_Th17_r2Pred.m](example_workflow_Th17_r2Pred.m)
	* includes calculation of R<sup>2</sup><sub>pred
* Modeling of time-series gene expression with linear differential equations: [example_workflow_Th17_timeLag.m](example_workflow_Th17_timeLag.m), as in [Bonneau et al. (2006) Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2006-7-5-r36)
* Combine TRNs using "max" or "mean" combine (e.g., models using TFA and TF mRNA): [combine_Th17_TRNs.m](combine_Th17_TRNs.m)
* Visualize inferred "core" networks [vis_th17_coreAnalysis.m](vis_th17_coreAnalysis.m)
	* Code to construct "core" networks is located [here](../scTRN)