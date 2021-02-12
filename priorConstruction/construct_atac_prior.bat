#BSUB -W 24:00
#BSUB -n 1
#BSUB -M 200000

module load bedtools/2.17.0
module load R/4.0.2

Rscript construct_atac_prior.R
