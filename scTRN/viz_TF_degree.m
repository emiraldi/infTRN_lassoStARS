%% viz_TF_degree
% visualize TF degree distributions with "visTfDegree_priorOverlaps.m"
%% References: 
% (1) Miraldi et al. (2018) "Leveraging chromatin accessibility for 
% transcriptional regulatory network inference in T Helper 17 Cells"
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: Nov. 15, 2018

clear all
close all
restoredefaultpath

%% Inputs:

currDir = '..';

addpath(fullfile(currDir,'infLassoStARS'))
addpath(fullfile(currDir,'glmnet'))
addpath(fullfile(currDir,'customMatlabFxns'))

netDir = '/Users/way5hl/Desktop/CCHMC/research/projects/scTRN/hMP/analysis/GRN/outputs';
netStatsFolder = 'netStats';

%netFiles = {'TRN_scRNA_Bryson_hMP_10X_QC_scater_N500M1_Norm_tpm_Dist_CIDR_Sub_sct1k_Clust_NN_size20_adapt0_sigma1_Norm_tpm_pseudo1_log2_prior_ATAC_Bryson_MP_MACS_MOODS_p5_sGB10kb_bias50_sp_pCorrCut0p01.tsv','Impute_NN20_TFA';
%'TRN_scRNA_Bryson_hMP_10X_QC_scater_N500M1_Norm_tpm_Dist_CIDR_Sub_sct1k_Clust_NN_size20_adapt0_sigma1_Norm_tpm_pseudo1_log2_prior_ATAC_Bryson_MP_MACS_MOODS_p5_sGB10kb_bias50_TFmRNA_sp_pCorrCut0p01.tsv','Impute_NN20_TFmRNA';
%'TRN_scRNA_Bryson_hMP_10X_QC_scater_N500M1_Norm_tpm_pseudo1_log2_prior_ATAC_Bryson_MP_MACS_MOODS_p5_sGB10kb_bias50_sp_pCorrCut0p01.tsv','log2_TPM_TFA';
%'TRN_scRNA_Bryson_hMP_10X_QC_scater_N500M1_Norm_tpm_pseudo1_log2_prior_ATAC_Bryson_MP_MACS_MOODS_p5_sGB10kb_bias50_TFmRNA_sp_pCorrCut0p01.tsv','log2_TPM_TFmRNA';
%};
netFiles = {
'TRN_scRNA_Bryson_hMP_10X_QC_scater_N500M1_Norm_tpm_pseudo1_log2_prior_ATAC_Bryson_MP_MACS_MOODS_p5_sGB10kb_bias50_sp_pCorrCut0p01.tsv','log2_TPM_TFA';
'TRN_scRNA_Bryson_hMP_10X_QC_scater_N500M1_Norm_tpm_pseudo1_log2_prior_ATAC_Bryson_MP_MACS_MOODS_p5_sGB10kb_bias50_TFmRNA_sp_pCorrCut0p01.tsv','log2_TPM_TFmRNA';
};

topN = 20; % number of Top-Degree TFs to include in the zoomed-in degree bar graph
tfList = '/Users/way5hl/Desktop/CCHMC/research/projects/scTRN/hMP/inputs/TF_scRNA_Bryson_hMP_10X_QC_scater_N500M1_union_ExpSCT_TFASCT_50pctl.txt';
%tfList_figSuffix = '_litTFs'; % file suffix for degree bar graph from TFs provided in tfList
tfList_figSuffix = '_PotRegTFs'; % file suffix for degree bar graph from TFs provided in tfList


%% END Inputs

netStatsOut = fullfile(netDir,netStatsFolder);
mkdir(netStatsOut)

totNets = size(netFiles,2);

%%
for nind = 1:totNets
    netTsv = netFiles{nind,1};
    disp(netTsv)
    netFile = fullfile(netDir,netTsv);
    netNickName = netFiles{nind,2};
    statsOut = fullfile(netStatsOut,[netNickName '_stats']);
    fullBarOut = fullfile(netStatsOut,[netNickName '_full']);
    topN_barOut = fullfile(netStatsOut,[netNickName '_top' num2str(topN)]);
    tfListOut = fullfile(netStatsOut,[netNickName tfList_figSuffix]);
    
    visTfDegree_priorOverlaps(netFile,statsOut,fullBarOut,topN,topN_barOut,...
        tfList, tfListOut);    
end
