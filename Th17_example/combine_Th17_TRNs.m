%% combine_Th17_TRNs
% max combine TFA and TF mRNA Th17 TRNs (15 TFs per gene) to get a 15 TFs per
% gene model

clear all
close all
restoredefaultpath

currDir = '..';

addpath(fullfile(currDir,'infLassoStARS'))
addpath(fullfile(currDir,'glmnet'))
addpath(fullfile(currDir,'customMatlabFxns'))

%% parameters
combinedNetTsv = 'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb_sp.tsv';

combineOpt = 'max';

meanEdgesPerGene = 15;

nets2combine = {'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50.mat';
    'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_TFmRNA.mat'};

combineTRNs(combinedNetTsv,combineOpt,meanEdgesPerGene,nets2combine)
