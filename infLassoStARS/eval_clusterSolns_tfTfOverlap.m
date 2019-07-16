function eval_clusterSolns_tfTfOverlap(tfPairMat, outDir, saveFig,...
    solnsOfInt)
%% GOAL: Evaluate cluster solutions / visualize silhouette scores versus 
% number of clusters to help choose number of TF-TF clusters
%% INPUTS:
% tfPairMat -- .mat file of TF-TF normalized overlaps and other statistics, 
%   as calculated by calc_zscoredTfTfOverlaps.m
% outDir -- directory for output figure folders (if saved)
% saveFig -- save figures? 1 --> yes, 0 --> no
% solnsOfInt -- row vector of cluster numbers to consider (OPTIONAL);
%   default is [2, 3,... .6*"total sig TFs"]

% parameters unlikely to change
xSize = 7;  % output .pdf figure x-dimension (inches)
ySize = 6;  % output .pdf figure y-dimension (inches)

%% load overlap analysis
load(tfPairMat)
zMatDist = tfPairAnal.zMatDist;

%% set cluster solution range, if needed
if nargin < 4
    disp('Using default cluster range: [2, 3,... .6*"total sig TFs"]')
    sigTfs = tfPairAnal.sigTfs;
    totTfs = length(sigTfs);
    solnsOfInt = 2:round(.6*totTfs); % vector of solutions considered
end
totSolns = length(solnsOfInt);
figBase = ['silhouette_' num2str(solnsOfInt(1)) '_to_' num2str(solnsOfInt(end))];

%% get clustering solns and silhouette scores   
pdis = tril(zMatDist,-1);
pdis = squareform(pdis);
link = linkage(pdis,'ward');    
solnMat = zeros(totTfs,totSolns);
silScoreMat = zeros(totTfs,totSolns);
aveSilScore = zeros(totSolns,1);
stdSilScore = zeros(totSolns,1);
for cind = 1:totSolns
    desClusts = solnsOfInt(cind); % desired clusters        
    currSoln = cluster(link,'maxclust',desClusts);
    [silScores,aveSilScore(cind),stdSilScore(cind)] = evalSilouetteDist(zMatDist,currSoln);
    solnMat(:,cind) = currSoln;
    silScoreMat(:,cind) = silScores;
%     disp(['Clust = ' num2str(desClusts) ', mean_{Sil} = ' num2str(mean(silScores))])
end

% boxplot
figure(1), clf
boxplot(silScoreMat,'Labels',num2str(solnsOfInt'))
xlabel('Number of Clusters','FontSize',12)
ylabel('Silhouette Scores','FontSIze',12)
set(gca,'FontSize',12)
plotRange = 5:5:solnsOfInt(end);
set(gca,'XTick',plotRange,'XTickLabel',num2str(plotRange'),'XTickLabelRotation',90)
grid on
title(['Silhouette Score'],'FontSize',14)

if saveFig
    subfigname = fullfile(outDir,[figBase '_boxplot']);
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]); 
    print('-painters','-dpdf','-r100',[subfigname '.pdf'])
    saveas(gcf,subfigname,'fig')
    disp(subfigname)
end

% means w/ standard deviation error bars
figure(2), clf
errorbar(solnsOfInt,aveSilScore,stdSilScore,'bo-','LineWidth',2)
axis tight
xlabel('Number of Clusters','FontSize',12)
ylabel('Silhouette Scores','FontSIze',12)
set(gca,'FontSize',12)
grid on
title(['Silhouette Score'],'FontSize',14)

if saveFig
    subfigname = fullfile(outDir, [figBase '_mean']);
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]); 
    print('-painters','-dpdf','-r100',[subfigname '.pdf'])
    saveas(gcf,subfigname,'fig')
    disp(subfigname)
end