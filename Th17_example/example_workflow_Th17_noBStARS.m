%% example_workflow_Th17_noBStARS
% Use mLASSO-StARS to build a TRN from gene expression and prior
% information in four steps. Please refer to each function's help
% annotations for descriptions of inputs, outputs and other information.
% This code does not run bStARS to define upper and lower bounds on lambda
% containing the target instability, but rather relies on a user-defined
% lambda range.
%% References: 
% (1) Miraldi et al. (2018) "Leveraging chromatin accessibility for 
% transcriptional regulatory network inference in T Helper 17 Cells"
% (2) Qian et al. (2013) "Glmnet for Matlab."
% http://www.stanford.edu/~hastie/glmnet_matlab/
% (3) Liu, Roeder, Wasserman (2010) "Stability Approach to Regularization 
%   Selection (StARS) for High Dimensional Graphical Models". Adv. Neural.
%   Inf. Proc.
% (4) Muller, Kurtz, Bonneau. "Generalized Stability Approach for Regularized
%   Graphical Models". 23 May 2016. arXiv.
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: March 30, 2018

clear all
close all

currDir = '..';

addpath(fullfile(currDir,'infLassoStARS'))
addpath(fullfile(currDir,'glmnet'))
addpath(fullfile(currDir,'customMatlabFxns'))

geneExprTFAdir = './outputs/processedGeneExpTFA';
mkdir(geneExprTFAdir)

%% 1. Import gene expression data, list of regulators, list of target genes
% into a Matlab .mat object
normGeneExprFile = './inputs/geneExpression/th17_RNAseq254_DESeq2_VSDcounts.txt';
targGeneFile = './inputs/targRegLists/targetGenes_names.txt';
potRegFile = './inputs/targRegLists/potRegs_names.txt';
tfaGeneFile = './inputs/targRegLists/genesForTFA.txt';
geneExprMat = fullfile(geneExprTFAdir,'geneExprGeneLists.mat');

disp('1. importGeneExpGeneLists.m')
importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,...
    tfaGeneFile,geneExprMat)

%% 2. Given a prior of TF-gene interactions, estimate transcription factor 
% activities (TFAs) using prior-based TFA and TF mRNA levels
priorName = 'ATAC_Th17';
priorFile = ['./inputs/priors/' priorName '.tsv']; % Th17 ATAC-seq prior
edgeSS = 50;
minTargets = 3;
[xx, priorName, ext] = fileparts(priorFile);
tfaMat = fullfile(geneExprTFAdir,[priorName '_ss' num2str(edgeSS) '.mat']);

disp('2. integratePrior_estTFA.m')
integratePrior_estTFA(geneExprMat,priorFile,edgeSS,...
     minTargets, tfaMat)

%% 3. Calculate network instabilities using a range of instabilities
tfaOpt = '';
lambdaBias = .5;
totSS = 50;
lambdaMin = .01;     % check estimateInstabilitiesTRN.m output figure to verify that target instability is reached
lambdaMax = 1;
logLambdaSteps = 25; % will have 1/logLambdaSteps per log10 step
subsampleFrac = .63;
leaveOutSampleList = '';
leaveOutInf = '';
netSummary = [priorName '_bias' strrep(num2str(lambdaBias),'.','p') tfaOpt];
instabilitiesDir = fullfile('./outputs',strrep(['instabilities_minL' num2str(lambdaMin) ...
    '_maxL' num2str(lambdaMax) '_SS' num2str(totSS) leaveOutInf],'.','p'));
mkdir(instabilitiesDir)
instabOutMat = fullfile(instabilitiesDir,netSummary);

disp('3. estimateInstabilitiesTRN.m')
estimateInstabilitiesTRN(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
    totSS,lambdaMin,lambdaMax,logLambdaSteps,subsampleFrac,instabOutMat,...
    leaveOutSampleList)

%% 4. For a given instability cutoff and model size, rank TF-gene
% interactions, calculate stabilities and network file for jp_gene_viz
% visualizations
priorMergedTfsFile = ['./inputs/priors/' priorName '_mergedTfs.txt'];
try % not all priors have merged TFs and merged TF files
    ls(priorMergedTfsFile) 
catch
    priorMergedTfsFile = '';
end
meanEdgesPerGene = 20;
targetInstability = .05;
networkDir = strrep(instabilitiesDir,'instabilities','networks');
instabSource = 'Network';
mkdir(networkDir);
networkSubDir = fullfile(networkDir,[instabSource ...
    strrep(num2str(targetInstability),'.','p') '_' ...
    num2str(meanEdgesPerGene) 'tfsPerGene']);
mkdir(networkSubDir)
trnOutMat = fullfile(networkSubDir,netSummary);
outNetFileSparse = fullfile(networkSubDir,[netSummary '_sp.tsv']);
networkHistDir = fullfile(networkSubDir,'Histograms');
mkdir(networkHistDir)
subsampHistPdf = fullfile(networkHistDir,[netSummary '_ssHist']);

disp('4. buildTRNs_mLassoStARS.m')
buildTRNs_mLassoStARS(instabOutMat,tfaMat,priorMergedTfsFile,...
    meanEdgesPerGene,targetInstability,instabSource,subsampHistPdf,trnOutMat,...
    outNetFileSparse)
