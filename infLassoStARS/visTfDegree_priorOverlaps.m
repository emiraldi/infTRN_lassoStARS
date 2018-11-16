function visTfDegree_priorOverlaps(netFile,statsOut,fullBarOut,topN,topN_barOut,...
  tfList, tfListOut)
%% function visTfDegree_priorOverlaps(netFile,statsOut,fullBarOut,topN,topN_barOut,
%   tfList, tfListOut)
%% GOAL: Visualize degree distribution of a TRN, as a horizontal bar graph
%   distinguishing (1) positive (red) and negative (blue) edges and (2)
%   prior-supported edges (pink, light blue). The resulting bar graph lacks
%   TF axis labels, but user can specify a value for topN >0 to generate a
%   zoomed-in view of the "Top N" highest degree TFs. The user can also
%   supply a text file list of TFs to visualize in a bar graph.
%% INPUTS:
%   netFile -- a sparse network file, e.g., expected output from 
%       buildTRNs_mLassoStARS.m, i.e., (top rows of example output):
%             TF	Target	SignedQuantile	NonzeroSubsamples(Stability)	pCorr	stroke	stroke-width	stroke-dasharray
%             Zhx1	D16Ertd472e	1	50.31	0.312	rgb(220,157,158)	2	None
%             Zkscan1	Il17f	-0.99738	50.28	-2.8E-01	rgb(157,181,227)	1.9994	2,2
%             Irf3	Il17f	-0.99477	50.26	-2.6E-01	rgb(160,182,227)	1.999	2,2
%          NOTE: last column is important: "2,2" denotes edge without prior
%           information, while "None" denotes prior information
%   statsOut -- a text file name for statistics about TRN (fraction of prior
%       edges retained in the final network, new edges learned, etc.)
%   fullBarOut -- specify a file name base for .fig + .pdf's generated of
%       the full degree distribution.
%       NOTE: provide empty string, if you don't want to save the figure
%   topN -- supply the number of top-degree TFs to visualize in a zoomed-in
%       version of the bar graph (includes labeling of TFs on Y axis).
%   topN_barOut -- specify a file name base for .fig + .pdf's of the zoomed
%       -in Top N bar graph.
%       NOTE: provide empty string, if you don't want to save the figure
%   tfList -- (optional) provide a text file of TFs you'd like to visualize 
%       w/ degree bar graph
%   tfListOut -- (optional) specify a file name base for .fig + .pdf's of 
%       the bar graph for selected TFs
%% OUTPUTS:
%   statsOut -- a text file of statistics about TRN (fraction of prior
%       edges retained in the final network, new edges learned, etc.)
%   fullBarOut -- .fig + .pdf's generated for the full degree distribution.
%   topN_barOut -- .fig + .pdf's of the Top N Highest Degree TFs
%   tfListOut -- .fig + .pdf's of degrees for selected TFs
%% Dependencies: customMatlabFxns
%% References: 
% Miraldi et al. (2018) bioRxiv "Leveraging chromatin accessibility for 
% transcriptional regulatory network inference in T Helper 17 Cells"
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: Nov. 15, 2018

%% begin DEBUG parameters:
% clear all
% close all
% restoredefaultpath
% currDir = '..';
% addpath(fullfile(currDir,'infLassoStARS'))
% addpath(fullfile(currDir,'glmnet'))
% addpath(fullfile(currDir,'customMatlabFxns'))
% currDir = '../scTRN';
% cd(currDir)
% netDir = 'outputs';
% netStatsFolder = 'netStats';
% netFiles = {'prior_miraldi_Th17_48h_cut4_bp10000_sATAC_p1Em5_huA_bias50_TFmRNA_sp.tsv','scMeth_TFmRNA';
%     'prior_miraldi_Th17_48h_cut4_bp10000_sATAC_p1Em5_huA_bias50_sp.tsv','scMeth_TFA'};
% topN = 20; % number of Top-Degree TFs to include in the zoomed-in degree bar graph
% tfList = 'inputs/geneLists/th17_literatureCore_TFs.txt';
% tfList_figSuffix = '_litTFs'; % file suffix for degree bar graph from TFs provided in tfList
% % END Inputs
% netStatsOut = fullfile(netDir,netStatsFolder);
% mkdir(netStatsOut)
% totNets = size(netFiles,2);
% nind = 1;
% netTsv = netFiles{nind,1};
% disp(netTsv)
% netFile = fullfile(netDir,netTsv);
% netNickName = netFiles{nind,2};
% statsOut = fullfile(netStatsOut,[netNickName '_stats.txt']);
% fullBarOut = fullfile(netStatsOut,[netNickName '_full']);
% topN_barOut = fullfile(netStatsOut,[netNickName '_top' num2str(topN)]);
% tfListOut = fullfile(netStatsOut,[netNickName tfList_figSuffix]);    
%% end DEBUG parameters

load cbcolors.mat
load linecolors.mat

%% Load inferred network
fid = fopen(netFile,'r');
% get first line and see what regulators we have
tline = fgetl(fid);    
totCols = length(cellstr(strsplit(tline,'\t')));
fclose(fid);
% get the rest of the data using textscan
fid = fopen(netFile,'r');

if totCols > 3
    C = textscan(fid,['%s%s%f' repmat('%s',1,totCols-3)],'Delimiter','\t','Headerlines',1);
else
    error('Network should be in sparse format (TF,gene,edge weight,...).')
end
fclose(fid);
infRegs = C{1};
infTargs = C{2};
infInts = [C{3}];
inPriorStr = C{8};
uniInfTargs = unique(infTargs);
uniInfRegs = unique(infRegs);
totTargs = length(uniInfTargs);
totRegs = length(uniInfRegs);

%% Determine number of + and - targets per TF, and whether interaction was in prior
negPrior = zeros(totRegs,1); 
negTotal = zeros(totRegs,1);
posPrior = zeros(totRegs,1); 
posTotal = zeros(totRegs,1);

inPrior = zeros(totRegs,1);
for tf = 1:totRegs
    currTf = uniInfRegs{tf};
    currInds = find(ismember(infRegs,currTf));
    currInts = infInts(currInds);
    currPriorInf = {inPriorStr{currInds}};
    currTargs = {infTargs{currInds}};
    totCurrTargs = length(currTargs);
    currPos = {currTargs{find(currInts>0)}};
    currNeg = {currTargs{find(currInts<0)}};
    % intersect with prior
    currTargsPrior = {currTargs{find(ismember(currPriorInf,'None'))}};
    if length(find(ismember(currPriorInf,'None'))) > 0
        inPrior(tf) = 1;
    end
    posTotal(tf) = length(currPos); negTotal(tf) = length(currNeg);
    posPrior(tf) = length(intersect(currPos,currTargsPrior));
    negPrior(tf) = length(intersect(currNeg,currTargsPrior));
end

% Print some stats
stats = '';
sum(negTotal);
ratioPosToNeg = sum(posTotal)/sum(negTotal);
totNet = sum(posTotal) + sum(negTotal);
stats = [stats num2str(totNet) ' Edges Total.\n'];
stats = [stats num2str(sum(posTotal)) ' + Edges\n' num2str(sum(negTotal)) ' - Edges.\n'];
stats = [stats roundstring2(ratioPosToNeg) ' is the ratio + to - Edges.\n'];
totPrior = sum(posPrior) + sum(negPrior);
stats = [stats num2str(totRegs) ' total TFs in the TRN.\n' num2str(sum(inPrior)) ' TFs were in the Prior.\n'];
stats = [stats num2str(totPrior) ' edges in the Prior.\n' num2str(round(totPrior*100/totNet)) ' percent of total edges are in the prior.\n'];
totNewWPrior = sum(posTotal(find(inPrior))) + sum(negTotal(find(inPrior))) - totPrior;
fractionNewEdgesPriorTFs = totNewWPrior / (totNewWPrior+totPrior);
stats = [stats num2str(totNewWPrior) ' new edges for TFs in the Prior.\n' num2str(round(100*fractionNewEdgesPriorTFs)) ' percent of edges for Prior TFs are new.\n'];
totNewNonPriorTfs = totNet - totPrior - totNewWPrior;
fractionTFsNotInPrior = totNewNonPriorTfs/ totNet;
stats = [stats num2str(totNewNonPriorTfs) ' new edges for TFs NOT in the Prior.\n' num2str(round(100*fractionTFsNotInPrior)) ' percent of total edges connect to TFs NOT in the Prior.\n'];
if statsOut
    fprintf(stats); % to the screen
    fid = fopen(statsOut,'w');
    fprintf(fid, stats);
    fclose(fid);
    disp('  ')
    disp([statsOut ' generated.'])
end


%% generate a figure
colors = [cbcolors(1,:); % red
cbcolors(9,:); % pink    
linecolors(1,:); % blue      
linecolors(11,:)]; % light blue
    
[tfsOrd,inds] = sort(sum([negTotal posTotal],2),'ascend');

barData = [(posTotal) posPrior -(negTotal)  -negPrior];
legendInf = {'+, not in Prior','+, w/Prior Support','-, not in Prior','-, w/Prior Support'};
barDataSorted = barData(inds,:);
figure
% subplot(4,1,1:3)
currFont = 14;
hold on
for bd = 1:4
    barh(barDataSorted(:,bd),'FaceColor',colors(bd,:))
end
axis tight
grid on, box on
set(gca,'FontSize',currFont)
legend(legendInf,'Location','SouthOutSide','FontSize',currFont+2)
xlabel('# of Negative and Positive Targets per TF','FontSize',currFont+2)
ylabel('TFs ordered by degree','FontSize',currFont+2)

currFig = fullBarOut;
if currFig
    disp('Improve figure dimensions manually, then hit enter.')
    pause
    save2pdf(currFig,gcf,200)
    saveas(gcf,currFig,'fig')
    disp(currFig)
end

%% zoom in on figure and add TF names
if topN
    barDataSorted = barData(inds(end-topN+1:end),:);
    figure(2), clf
    % subplot(4,1,1:3)
    currFont = 14;
    hold on
    for bd = 1:4
        barh(barDataSorted(:,bd),'FaceColor',colors(bd,:))
    end
    axis tight
    grid on, box on
    set(gca,'FontSize',currFont)
    legend(legendInf,'Location','SouthOutSide','FontSize',currFont+2)
    xlabel('# of Negative and Positive Targets per TF','FontSize',currFont+2)
    set(gca,'YTick',1:topN,'YTickLabel',upper({uniInfRegs{inds(end-topN+1:end)}}))

    currFig = topN_barOut;
    if currFig
        disp('Improve figure dimensions manually, then hit enter.')
        pause
        save2pdf(currFig,gcf,200)
        saveas(gcf,currFig,'fig')
        disp(currFig)
    end
end
    
%% Plot a subset of TFs in specified order:
if tfList
    fid = fopen(tfList,'r');
    C = textscan(fid,'%s','Headerlines',0);
    fclose(fid);
    tfOrd = flipud(C{1}); 
    
    totTfs = length(tfOrd);    
    tfInds = zeros(totTfs,1);
    for tt = 1:totTfs
        tfInds(tt) = find(ismember(uniInfRegs,tfOrd{tt}));
    end
    
    barDataSorted = barData(tfInds,:);
    figure(3), clf
    subplot(1,2,2)
    currFont = 14;
    hold on
    for bd = 1:4
        barh(barDataSorted(:,bd),'FaceColor',colors(bd,:))
    end
    
    axis tight
    ax = axis();
    axis([ax(1:2) 0 totTfs+1])    
    grid on, box on
    grid minor
    set(gca,'FontSize',currFont)
%     legend(legendInf,'Location','SouthEast','FontSize',currFont+2)
    xlabel('# of - and + Targets per TF','FontSize',currFont+2)
    set(gca,'YTick',1:totTfs,'YTickLabel',upper(tfOrd))   
    currFig = tfListOut;
    if currFig
        disp('Improve figure dimensions manually, then hit enter.')
        pause
        save2pdf(currFig,gcf,200)
        saveas(gcf,currFig,'fig')
        disp(currFig)
    end    
end
   
    
