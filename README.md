# infTRN_lassoStARS

This repository contains a workflow for inference of transcriptional regulatory networks (TRNs) from gene expression data and prior information, as described in:

[Miraldi et al., Leveraging chromatin accessibility for transcriptional regulatory network inference in T Helper 17 Cells](https://www.biorxiv.org/content/early/2018/04/05/292987).

From gene expression data and tables of prior information, the [example Th17 workflow](Th17_example/example_workflow_Th17.m) can be used to infer a TRN using modified LASSO-StARS, and relies upon [GlmNet in Matlab](https://web.stanford.edu/~hastie/glmnet_matlab/index.html) to solve the LASSO.

The resulting network can be visualized with TRN visualization software: [jp_gene_viz](https://github.com/simonsfoundation/jp_gene_viz).

The workflow also includes TRN model evaluation based on precision-recall and ROC.

