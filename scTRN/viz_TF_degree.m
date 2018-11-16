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

netDir = 'outputs';
netStatsFolder = 'netStats';

netFiles = {'prior_miraldi_Th17_48h_cut4_bp10000_sATAC_p1Em5_huA_bias50_TFmRNA_sp.tsv','scMeth_TFmRNA';
    'prior_miraldi_Th17_48h_cut4_bp10000_sATAC_p1Em5_huA_bias50_sp.tsv','scMeth_TFA'};

topN = 20; % number of Top-Degree TFs to include in the zoomed-in degree bar graph
tfList = 'inputs/geneLists/th17_literatureCore_TFs.txt';
tfList_figSuffix = '_litTFs'; % file suffix for degree bar graph from TFs provided in tfList


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
