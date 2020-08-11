
We estimate TF-target gene interactions by associating putative TF binding sites within regions of open chromatin to nearby genes. This transcriptional regulatory network (TRN) serves as an input to the Inferelator.

Use the [construct_atac_prior.R](construct_atac_prior.R) script to construct a prior TRN from ATAC-seq peaks.

Prior construction requires the following:
1. bed file of accessible chromatin regions (ie, peaks)
2. A directory of TF motifs specified as position frequency matrices
3. bed file of genomic features (eg, TSS, gene body)

Follow the format of the example inputs in the [example_prior_inputs](example_prior_inputs) directory. See the [example_prior_outputs](example_prior_outputs) for example outputs.

R packages requirements:
* [motifmatchr](https://www.bioconductor.org/packages/release/bioc/html/motifmatchr.html) is a wrapper for the [MOODS](https://github.com/jhkorhonen/MOODS) motif scanning library.
* [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
* [TFBSTools](http://bioconductor.org/packages/release/bioc/html/TFBSTools.html)
* [BSgenome](https://www.bioconductor.org/packages/release/bioc/html/BSgenome.html) object for the necessary species. Here are a few commonly used genome builds:
    * [mm10](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Mmusculus.UCSC.mm10.html)
    * [hg19](http://www.bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html)
    * [hg38](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html)

Other requirements:
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
