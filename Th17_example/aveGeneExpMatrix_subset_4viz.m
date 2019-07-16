%% aveGeneExpMatrix_subset_4viz
%% Subset and average particular conditions (or cells) from the gene
% expression data matrix for visualization. 
% * User specifies new names for conditions and experiments to average
% * User specifies the order of select conditions and identifies clusters
% of conditions (e.g., for visualization in heatmaps in vis_allTfTfOverlaps.m) 
%% OUTPUTS:
% 1. .txt file of gene expression values, with averaged conditions
%   ordered as specified. (e.g., can serve as input to jp_gene_viz)
% 2. .mat file of gene expresison values with cluster info (e.g., for 
%   visualization in heatmaps in vis_allTfTfOverlaps.m) 
%% NOTE: The output gene expression matrix is limited to TFA genes.

clear all
close all

%% INPUTS
matlabDir = '..';
addpath(fullfile(matlabDir,'customMatlabFxns'))

geneExprTFAdir = './outputs/processedGeneExpTFA';
outDir = fullfile(geneExprTFAdir,'geneExpHeatmapInputs');
geneExprMat = fullfile(geneExprTFAdir,'geneExprGeneLists.mat');

%% condGroups -- each row is 3D
% Col 0 = a field not used in this code
% Col 1 = name for group of conditions -- will form output file name
% Col 2 = table with three columns (1 --> add a group line on heatmap
%    2 --> a name for set of samples, 3--> cell containing individual
%    sample names (as found in conditionsc)
condGroups = {...
    {'Th0_Th17_48hTh';...
        {1,'media(1h)',{'media_only_nCD4_T_cell_wt_1h_SL2653'};
        1,'Th0(1h)',{'Th0_wt_1h_SL2654'};
        0,'Th0(3h)',{'Th0_wt_3h_SL2655'};
        0,'Th0(6h)',{'Th0_wt_6h_SL2656'};
        0,'Th0(16h)',{'Th0_wt_16h_SL2778';'Th0_wt_16h_SL8352';'Th0_wt_16h_SL8354'};
        0,'Th0(48h)',{'Th0_wt_48h_SL2673';'Th0_RORg_wt_pp_2_48h_SL1843';'Th0_RORg_wt_pp_3_48h_SL2679';'Th0_STAT3_wt_pCD4-Cre_7_48h_SL1847';'Th0_STAT3_wt_pCD4-Cre_8_48h_SL3540'};        
        1,'Th17(1h)',{'Th17_wt_1h_SL1851'};
        0,'Th17(3h)',{'Th17_wt_3h_SL1852'};
        0,'Th17(6h)',{'Th17_wt_6h_SL1853'};
        0,'Th17(9h)',{'Th17_wt_9h_SL1854'};
        0,'Th17(12h)',{'Th17_wt_12h_SL1855'};
        0,'Th17(16h)',{'Th17_wt_16h_SL1856';'Th17_RORg_wt_pp_1_16h_SL599';'Th17_wt_16h_SL8353';'Th17_wt_16h_SL8355'};
        0,'Th17(24h)',{'Th17_wt_24h_SL1857'};
        0,'Th17(48h)',{'Th17_wt_48h_SL1858';'Th17_RORg_wt_pp_2_48h_SL1844';'Th17_RORg_wt_pp_3_48h_SL2680';'Th17_STAT3_wt_pCD4-Cre_7_48h_SL1848';'Th17_STAT3_wt_pCD4-Cre_8_48h_SL3541'};%;'48hr_Th17_veh_1';'48hr_Th17_veh_2'}}};
        1,'Th1(48h)',{'Th1_IMDM_wt_48h_SL2683';'Th1_IMDM_wt_48h_SL3304'};
        1,'Th2(48h)',{'Th2_IMDM_wt_48h_SL2684';'Th2_IMDM_wt_48h_SL3306'};
        1,'Treg(48h)',{'iTreg_IMDM_wt_48h_SL2685';'iTreg_IMDM_wt_48h_SL3305'}}};
    };

%% END INPUTS

mkdir(outDir)
disp(outDir)

totGroups = length(condGroups);
load(geneExprMat)

for groupInd = 1:totGroups
    condGroup = condGroups{groupInd};
    condGroupName = condGroup{1}; condInf = condGroup{2};
    startSpotsConds = [condInf{:,1}];
    condPrintNames = {condInf{:,2}}';
    indConds = {condInf{:,3}}';

    outBase = fullfile(outDir,[condGroupName]);

    %% Limit gene expression matrix to specific conditions
    aveCounts = [];      % averaged conditions (corresponding to condPrintNames)
    totBigConds = length(condPrintNames);
    aveGeneExprVals = [];
    for cind = 1:totBigConds
        currSamps = indConds{cind}';
        totIndConds = length(currSamps);    
        sampInds = zeros(totIndConds,1);
        for icind = 1:totIndConds
            sampInds(icind) = find(ismember(conditionsc,currSamps{icind}));
        end
        currVals = tfaGeneMat(:,sampInds);
        aveCounts = [aveCounts median(currVals,2)];
    end

    %% generate output text file
    outExprFile = [outBase '.txt'];
    fout = fopen(outExprFile,'w');
    fprintf(fout,['\t' strjoin(condPrintNames,'\t') '\n']);
    for gene = 1:length(tfaGenes)
        fprintf(fout,[tfaGenes{gene} '\t' ...
            strjoin(cellstr(num2str(aveCounts(gene,:)')),'\t') '\n']);
    end
    disp([outExprFile ' generated.'])
    disp([num2str(length(tfaGenes)) ' TFA genes included in output matrix.'])
    
    %% generate output .mat file
    genesc = tfaGenes;
    totConds = length(condPrintNames);
    zAveCounts = zscore(aveCounts')';
    meanC = mean(2.^aveCounts,2);
    l2AveCounts = log2((2.^aveCounts)./repmat(meanC,1,totConds));


    matOut = [outBase '.mat'];
    save(matOut,...
        'zAveCounts',...
        'l2AveCounts',...
        'meanC',...
        'aveCounts',...
...        'conds2average',...
        'condPrintNames',...
        'startSpotsConds',...
        'genesc')
    disp([matOut ' created.'])

end
