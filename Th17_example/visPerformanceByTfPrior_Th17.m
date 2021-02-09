%% visPerformanceByTfPrior_Th17
% visualize per-TF AUPRs in a heatmap (e.g., as in Fig. 2C of Miraldi et 
% al. (2019) Genome Research). 
% Heatmap columns are TF AUPRs organized according to gold standard (gold
%   standards results are separated using solid lines)
% Rows can be grouped based on prior network (can be separated by solid lines)
%   then TF mRNA, then TFA (can be separated by dotted lines)
% Based on: aggregatePerformanceByTfPrior_Th17_greyChIP_lim_1803.m and
%   visPerformanceByTfPrior_Th17_greyChIP_1803.m 
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and 
%   Biomedical Informatics, Cincinnati Children's Hospital
%% Date: Feb. 8, 2021

clear all
close all
restoredefaultpath

matlabDir = '..';

addpath(fullfile(matlabDir,'infLassoStARS'))
addpath(fullfile(matlabDir,'glmnet'))
addpath(fullfile(matlabDir,'customMatlabFxns'))

%% INPUTS

%% Gold standard(s) info
% column 1: nickname for G.S. (to be used in outputs)
% column 2: location of precision-recall outputs for that gold standard
% (e.g., from calcPRinfTRNs.m)

gsInf = {'KO' 'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/PR_KO';
    'KC1p5' 'outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/PR_KC1p5'};

totGSs = length(gsInf);

%% Info about the networks

%% Network info
% column 1: label for the network, 
% column 2: file base for the network (e.g., as
% given by netSummary variable in example_workflow_Th17.m and likely used
% for downstream outputs like P-R)
% column 3: what kind of line (solid '-', dotted ':', or no line '') should
%   be used underneath the network in that row of the heatmap. 
networkInf = {'No Prior'	'ATAC_Th17_bias100_TFmRNA' '-';
    'ATAC Th17 b=.5, mRNA'	'ATAC_Th17_bias50_TFmRNA' '';
    'ATAC Th17 b=.25, mRNA'	'ATAC_Th17_bias25_TFmRNA' ':';
    'ATAC Th17 b=.5, TFA' 'ATAC_Th17_bias50' '';
    'ATAC Th17 b=.25, TFA'	'ATAC_Th17_bias25' '-';
    'ENCODE b=.5, mRNA'	'ENCODE_bias50_TFmRNA' '';
    'ENCODE b=.25, mRNA'	'ENCODE_bias25_TFmRNA' ':';
    'ENCODE b=.5, TFA'	'ENCODE_bias50' '';
    'ENCODE b=.25, TFA'	'ENCODE_bias25' '-'};
totNets = size(networkInf,1);

% specify output info
fontSize = 10;
outFileBase = ['outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/' 'relAUPRperTF_heatmap'];
% output file size dimensions (in inches)
xSize = 11;
ySize = 8.5;

%% END INPUTS

netNames = {networkInf{:,1}}';
lineLocs = find(ismember({networkInf{:,3}}','-'));  % locations of lines
dottedLineLocs = find(ismember({networkInf{:,3}}',':'));  % for, e.g., marking transition between no TFA --> TFA

% initialize AUPR objects
auprs = [];
auprsRand = [];
for gs = 1:totGSs
    auprs(gs).name = gsInf{gs,1};
    auprs(gs).outDir = gsInf{gs,2};
    auprs(gs).regs = '';
    auprs(gs).auprs = [];
    auprs(gs).auprsRands = [];
end

%% stack AUPRs per TF into rows of the AUPR matrix
for nind = 1:totNets
    currNet = networkInf{nind,2};
    
    % initialize AUPR objects
    if nind == 1
        for gs = 1:totGSs
            currPRresults = fullfile(auprs(gs).outDir,[currNet '.mat']);
            load(currPRresults)
            % get general info about the AUPRS
            auprs(gs).regs = gsInfs.regs;           
            auprs(gs).auprsRand = gsInfs.randAuprByTf;
            % add the AUPR for this method
        end
    end
    for gs = 1:totGSs
        currPRresults = fullfile(auprs(gs).outDir,[currNet '.mat']);
        load(currPRresults)
        auprs(gs).auprs = [auprs(gs).auprs; gsInfs.auprsByTf'];
        clear gsInfs 
    end
end

%% cluster TFs (within G.S.'s) and visualize relative AUPRs log2(AUPR / random AUPR)
relAuprMat = []; % table relative AUPRs
gsLineLocs = []; % location of solid line to separate G.S.'s
regs = ''; % for clustered TF columns per G.S.
for gs = 1:totGSs    
    auprMat = auprs(gs).auprs;
    regsTmp = auprs(gs).regs;        
    totRegsTmp = length(regsTmp);      
    relAuprMatTmp = log2(auprMat./repmat(auprs(gs).auprsRand',totNets,1));
            
    % cluster within gold standards
    pdis = pdist(relAuprMatTmp');
    link = linkage(pdis,'ward');
    [h t horderRegs] = dendrogram(link,0);
    regs = strvcat(regs,strvcat(regsTmp{horderRegs}));
    gsLineLocs = [gsLineLocs; length(regs)];
    relAuprMat = [relAuprMat relAuprMatTmp(:,horderRegs)];
end
totRegs = length(regs);

colorMax = max(abs(relAuprMat(:))); 
figure(1), clf
subplot(1,4,2:4)
imagesc(relAuprMat)
colormap redblue
set(gca,'XTick',1:totRegs,'XTickLabel',upper(regs),'XTickLabelRotation',90)
set(gca,'YTick',1:size(relAuprMat,1),'YTickLabel',netNames)
set(gca,'CLim',[-colorMax colorMax],'TickLength',[0 0],'FontSize',fontSize)
title(['log_2(AUPR_{model}/AUPR_{random})'])
colorbar
axis image
ax = axis();
hold on, % add lines separating models and gold standards
for clust = 1:length(lineLocs)-1
    plot([ax(1) ax(2)], [lineLocs(clust) lineLocs(clust)]+.5,'k',...
        'LineWidth',.5)
end
for clust = 1:length(dottedLineLocs)
    plot([ax(1) ax(2)], [dottedLineLocs(clust) dottedLineLocs(clust)]+.5,'k:',...
        'LineWidth',.5)
end
for gs = 1:length(gsLineLocs)
    plot(gsLineLocs(gs)*[1 1]+.5, ax(3:4),'k','LineWidth',.5)
end

figInf = outFileBase;
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
print('-painters','-dpdf','-r200',figInf)
saveas(gcf,figInf,'fig')
disp(figInf)

set(gca,'CLim',[-2 2])

figInf = [outFileBase '_clim2'];
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
print('-painters','-dpdf','-r200',figInf)
saveas(gcf,figInf,'fig')
disp(figInf)

