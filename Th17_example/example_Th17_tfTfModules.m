%% example_Th17_tfTfModules.m
% example workflow to find clusters of TFs ("TF-TF modules" with shared 
% gene regulatory programs
%% 1. Calculate normalized overlap in target genes between TFs
%% 2. Evaluate cluster solutions / visualize silhouette scores versus 
% number of clusters to help choose number of TF-TF clusters
%% 3. Cluster and visualize full matrix of z-scored, normalized overlaps,
% optionally overlay with gene expression and core TF annotations
%% 4. Cluster and visualize "Top N" TF-TF clusters, where clusters are
% empirically ranked according to target overlap between cluster TF members
% and size of the cluster (see Miraldi et al. 2019. Genome Research.)
%% 5. Output "Top N" TF-TF clusters as networks for visualization in 
% jp_gene_viz (see https://github.com/flatironinstitute/ILCnetworks and the
% notebook ILC_TRN_Notebooks/TfTfModules-Positive-c57.ipynb); several
% network construction rules are used (1. TFs and the union of their
% targets, 2. TFs and genes that are targets of >50% of the TFs, and 3. TFs
% and genes that are targets of at least two TFs.)
%% References: 
% (1) Miraldi et al. (2019) "Leveraging chromatin accessibility for 
% transcriptional regulatory network inference in T Helper 17 Cells"
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: July 5, 2019

clear all
close all
restoredefaultpath

matlabDir = '..';

addpath(fullfile(matlabDir,'infLassoStARS'))
addpath(fullfile(matlabDir,'customMatlabFxns'))

%% Inputs from other workflows:
%% Network
% Input network below combines TF mRNA and TFA models (via
%   combine_Th17_TRNs.m) and was filtered so that |partial corr| between TFs
%   and target genes is > .01 (via filter_Th17_TRNs_by_pcorr.sh)
% Network should be in tab-delimited "table" format for step 1; columns 
%   correspond to TF regulators, rows correspond to target genes,
%   and values correspond to signed interactions
inputNetwork = 'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb/ATAC_Th17_bias50_maxComb_cut01.tsv';
% "Sparse" network format is required for step 5 (network visualization of
%   TF-TF modules, col 1 = TF, col 2 = target gene, col 3 = signed weight
%   of interaction, additional columns are optional
inputNetworkSparse = 'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb/ATAC_Th17_bias50_maxComb_cut01_sp.tsv';
%% (Optional) a subset of averaged gene expression coditions, used to
% visual TF mRNA expression in step 3. See aveGeneExpMatrix_subset_4viz.m
% TO OMIT: set aveGeneExprMat = ''
aveGeneExprMat = './outputs/processedGeneExpTFA/geneExpHeatmapInputs/Th0_Th17_48hTh.mat'; 
%% (Optional) tables of TF target enrichments in genesets, e.g., as generated by
% tfTargets_GSEA.sh, table contains p-vals, rows = gene sets, columns are
% genesets. Good place to input "core" TF info. Given the large 
% number of gene sets from other analyses (involving Kegg, etc., as in 
% scTRN/tfTargets_GSEA_loop.sh) that info would be better integrated into 
% the analysis using scTRN/visGSEAenrich_heatmaps.m with TFs outputted from
% step 3
% TO OMIT: set annotations = ''
annotations = {... N x 3 cell, were each row corresponds to an enrichment analysis:
    ...col 1: sign of enriched regulatory edges (1 --> positive, -1 --> inhibition)
    ...col 2: nickname for enrichment anlaysis (goes in title & figure file)
	...col 3: tab-delimited table of p-values associated with enrichments (sets X TFs)
    ...col 4: cell with the desired ordered of set annotations -- if desired
    ...   order is unknown, this can be left as an empty string
    1,'posEdgeCore',...
    'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb/GSEA/ATAC_Th17_bias50_maxComb_cut01_Th17set_Praw0p1_dir_wCut0p0_minSet5/Th17set_praw10_up_adjp.txt',...
    {'Th17 upOnly'; 'Th17 downOnly'};
    ... Annotation 2: TFs core due to negative edges
    -1,'negEdgeCore',....
    'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb/GSEA/ATAC_Th17_bias50_maxComb_cut01_Th17set_Praw0p1_dir_wCut0p0_minSet5/Th17set_praw10_down_adjp.txt',...
    {'Th17 upOnly'; 'Th17 downOnly'};
    };

%% 1. Calculate normalized overlap in target genes between TFs
tfTargMin = 20;     % only consider TFs with at least this number of targets
targTfMin = 1;      % only consider targets with at least this number of TFs
fdrCut = .1;         % cutoff for TF pair inclusion
edgeOpt = 'comb';   % can be set to one of three options:
%   'pos' -- limits analysis to positive regulatory interactions
%   'neg' -- limits anlaysis to negative regulatory interactions
%   'comb' -- integrates both positive and negative regulatory
%       interactions. For example, if TF A and B both negatively regulate
%       Gene 1, that regulatory interaction contributes to overlap,
%       while, if Gene 2 is positively regulated by TF A and negatively
%       regulated by TF B, that regulatory interaction will not contribute
%       to target overlap between TFs (because sign is different).
[outDirBase,fileName,ext] = fileparts(inputNetwork);
outDir = fullfile(outDirBase,fileName,strjoin({'zOverlaps',...    
    [edgeOpt 'Edge'],...
    ['fdr' num2str(100*fdrCut)],...
    ['tfMin' num2str(tfTargMin)],...
    ['targMin' num2str(targTfMin)]},'_'));

disp(' 1. Calculate z-scored overlap TF regulatory interactions')
tfPairMat = calc_zscoredTfTfOverlaps(inputNetwork,tfTargMin,...
    targTfMin, fdrCut, edgeOpt, outDir);
% tfPairMat = 'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb/ATAC_Th17_bias50_maxComb_cut01/zOverlaps_combEdge_fdr10_tfMin20_targMin1/tfPair.mat';

%% 2. Evaluate cluster solutions / visualize silhouette scores versus 
% number of clusters to help choose number of TF-TF clusters
saveFig2 = 0;       % save figures? 1 --> yes, 0 --> no

disp('2. Evaluate cluster solutions using silhouette score.')
eval_clusterSolns_tfTfOverlap(tfPairMat, outDir, saveFig2)

%% 3. Cluster and visualize full matrix of z-scored, normalized overlaps,
% optionally overlay with gene expression and core TF annotations
datasetName = 'Th17';
desClusts = 50;     % desired number of clusters (clustering is 
%   hierarchical, so affects visualization only)
titleInf = [datasetName ', ' edgeOpt '-edge , FDR: ' num2str(100*fdrCut)...
    '%, min TF, target: ' num2str(tfTargMin) ', ' num2str(targTfMin)...
    ', clust = ' num2str(desClusts)]; % 
fullFigOutDir3 = fullfile(outDir,['fullFigs_clust' num2str(desClusts)]);
saveFig3 = 0;       % save figures? 1 --> yes, 0 --> no
axisFontSize3 = 2;   % heatmap fontsize
xSize = 7;
ySize = 12;
figureDimensions = [xSize ySize]; % 2D vector of (x-size, y-size) dimensions for
 %   saving pdfs of the figures, units are inches

disp('3. Cluster and visualize all TF-TF overlaps')
vis_allTfTfOverlaps(datasetName, tfPairMat, desClusts, titleInf,... 
    fullFigOutDir3, aveGeneExprMat, annotations, saveFig3, axisFontSize3,...
    figureDimensions)

%% 4. Cluster and visualize "Top N" TF-TF clusters, where clusters are
% empirically ranked according to target overlap between cluster TF members
% and size of the cluster (see Miraldi et al. 2019. Genome Research.)

saveFig4 = 1;       % save figures? 1 --> yes, 0 --> no
axisFontSize4 = 4;  % heatmap fontsize
topN = 15;          % number of top-ranking TF-TF clusters for heatmap viz
fullFigOutDir4 = fullfile(outDir,['Top' num2str(topN) '_Figs_clust' num2str(desClusts)]);

disp('4. Cluster and visualize "Top N" TF-TF overlaps')
vis_topN_TfTfOverlaps(datasetName, tfPairMat, desClusts, topN,...
    titleInf, fullFigOutDir4, aveGeneExprMat, annotations, saveFig4, ....
    axisFontSize4, figureDimensions)

%% 5. Output "Top N" TF-TF clusters as networks for visualization in 
% jp_gene_viz (see https://github.com/flatironinstitute/ILCnetworks and the
% notebook ILC_TRN_Notebooks/TfTfModules-Positive-c57.ipynb); several
% network construction rules are used (1. TFs and the union of their
% targets, 2. TFs and genes that are targets of >50% of the TFs, and 3. TFs
% and genes that are targets of at least two TFs.)

netOutDir = fullfile(outDir,['Top' num2str(topN) '_Networks_clust' num2str(desClusts)]);

disp('5. Generate subnetwork files for each of the "Top N" TF-TF clusters.')
output_topN_tfTfclusters_jp_gene_viz(inputNetworkSparse,...
    tfPairMat, topN, desClusts, netOutDir)