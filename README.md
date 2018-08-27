# infTRN_lassoStARS

This repository contains a workflow for inference of transcriptional regulatory networks (TRNs) from gene expression data and prior information, as described in:

[Miraldi et al., Leveraging chromatin accessibility for transcriptional regulatory network inference in T Helper 17 Cells](https://www.biorxiv.org/content/early/2018/04/05/292987).

From gene expression data and tables of prior information, the [example Th17 workflow](Th17_example/example_workflow_Th17.m) can be used to infer a TRN using modified LASSO-StARS, and relies upon [GlmNet in Matlab](https://web.stanford.edu/~hastie/glmnet_matlab/index.html) to solve the LASSO. Workflow also includes TRN model evaluation based on precision-recall and ROC.

The resulting network can be visualized with TRN visualization software: [jp_gene_viz](https://github.com/simonsfoundation/jp_gene_viz).

Additional workflows are provided for:
* [Out-of-sample gene expression prediction](Th17_example/example_workflow_Th17_r2Pred.m), including calculation of R<sup>2</sup><sub>pred
* [Modeling of time-series gene expression with linear differential equations](Th17_example/example_workflow_Th17_timeLag.m), as in [Bonneau et al. (2006) Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2006-7-5-r36)
