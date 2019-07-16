function output_topN_tfTfclusters_jp_gene_viz(inputNetworkSparse,...
    tfPairMat, topN, desClusts, netOutDir)
%% output_topN_tfTfclusters_jp_gene_viz(inputNetworkSparse,...
%     tfPairMat, topN, desClusts, netOutDir)
%% GOAL: Get sparse networks for the "Top N" most significant TF-TF modules. 
% Three networks will be constructed for each TF-TF module, with the 
% following rules for including target genes:
%   1. Union of all TF module members' target genes
%   2. target Gene must be in >1 TF module members' target genes
%   3. target Gene must be in > 50% of module members' target genes
%% INPUTS:
% inputNetworkSparse -- network in "sparse" format, where col 1 = TF, 
%   col 2 = target gene, col 3 = signed weight of interaction (additional 
%   columns are optional)
% tfPairMat -- .mat file of TF-TF normalized overlaps and other statistics, 
%   as calculated by calc_zscoredTfTfOverlaps.m
% topN -- number of "most significant" clusters to visualize
% desClusts -- desired number of clusters / clustering solution
% netOutDir -- output directory for TF-TF module subnetworks
%% OUTPUTS:
% sparse-format networks for the "Top N" most significant TF-TF modules
%   (see GOAL) for details

%% Load inferred network w/ prior information
mkdir(netOutDir)
disp(netOutDir)

fid = fopen(inputNetworkSparse,'r');
% get first line and see what regulators we have
header = [fgetl(fid) '\n'];    
totCols = length(cellstr(strsplit(header,'\t')));
fclose(fid);
% get the rest of the data using textscan -- note treating everything as
% strings
fid = fopen(inputNetworkSparse,'r');
if totCols == 3
    infNet = textscan(fid,['%s%s%s'],'Delimiter','\t','Headerlines',1);
elseif totCols > 3
    infNet = textscan(fid,['%s%s%s' repmat('%s',1,totCols-3)],'Delimiter','\t','Headerlines',1);
else
    error('Prior matrix should be in sparse format (TF,gene,edge weight).')
end
fclose(fid);
infRegs = infNet{1};
infTargs = infNet{2};

load(tfPairMat)
%% get clustering solutions     
zMatDist = tfPairAnal.zMatDist;
zMat = tfPairAnal.zMat;
sigTfs = tfPairAnal.sigTfs;
totTfs = length(sigTfs);
pdis = tril(zMatDist,-1);
pdis = squareform(pdis);
link = linkage(pdis,'ward');    
currSoln = cluster(link,'maxclust',desClusts);

solnNstats = [];
solnNstats.clustInds = cell(desClusts,1);
solnNstats.clustMemNames = cell(desClusts,1);
solnNstats.clustMedzOv = zeros(desClusts,1);
solnNstats.clustMeanzOv = zeros(desClusts,1);
solnNstats.clustMadPraw = zeros(desClusts,1);
solnNstats.soln = currSoln;

zMatDistribution = squareform(zMat);    
totPairs = length(zMatDistribution);
clustPvals = zeros(desClusts,1);
clustSizes = zeros(desClusts,1);
for clust = 1:desClusts
    clustInds = find(currSoln == clust);
    clustSize = length(clustInds);       
    clustSizes(clust) = clustSize;
    solnNstats.clustInds{clust} = clustInds;
    solnNstats.clustMemNames{clust} = strjoin(cellstr(strvcat(sigTfs{clustInds})),'_');
    currZmat = squareform(zMat(clustInds,clustInds));
    currPairs = length(currZmat);
    currPs = zeros(currPairs,1);
    for pind = 1:currPairs %      currMed = median(currZmat(pind,setdiff(1:clustSize,pind))); % calculate median
        currPs(pind) = length(find(zMatDistribution >= currZmat(pind)))/totPairs;
    end        
    % combine the p-values
    weights = (clustSize/currPairs)*ones(currPairs,1); % these weights with weighted z-method will count each TF vs. TF pairs
    clustPvals(clust) = combineP_stouffersZMethod(currPs,weights);        
end   
solnNstats.clustSizes = clustSizes;

[psOrd, inds] = sort(clustPvals,'ascend');

% generate subnetworks for each of the top N clusters
clusts2keep = inds(1:topN);
for kc = 1:topN
    tfInds = find(currSoln==clusts2keep(kc));
    tfNames = {sigTfs{tfInds}};
    currInts = tfPairAnal.sigRegInts(:,tfInds);
    tfsPerTarg = sum(sign(abs(currInts)),2);
    uniTargs = unique(tfPairAnal.sigTargs); 
    
    unionTargs = find(tfsPerTarg);
    p50Targs = find(tfsPerTarg > length(tfInds)/2);
    min2Targs = find(tfsPerTarg > 1);

    currFileBase = strjoin(tfNames,'_');
    foutUnion = fopen(fullfile(netOutDir,[currFileBase '_unionTargs_sp.tsv']),'w'); % include target in network if at least 1 TF regulates it
    fout50per = fopen(fullfile(netOutDir,[currFileBase '_p50Targs_sp.tsv']),'w');    % > 50% of TFs must regulate it
    foutMin2 = fopen(fullfile(netOutDir,[currFileBase '_min2Targs_sp.tsv']),'w');   % at least 2 TFs must regulate it
    fprintf(foutUnion,header);
    fprintf(foutMin2,header);
    fprintf(fout50per,header);
    for ti = 1:length(unionTargs)
        currTargInd  = unionTargs(ti);
        targName = tfPairAnal.sigTargs{currTargInd};
        textOut = '';
        targInts = currInts(currTargInd,:);
        currTfInds = find(targInts);
        totTargTfs = length(currTfInds);
        for tf = 1:totTargTfs
            tfInd = currTfInds(tf);
            tfName = tfNames{tfInd};
            infNetInd = intersect(find(ismember(infRegs,tfName)),...
                find(ismember(infTargs,targName)));
            outList = cell(totCols,1);
            for icount = 1:totCols
                outList{icount} = infNet{icount}{infNetInd};
            end
            textOut = [textOut strjoin(outList,'\t') '\n' ];
        end
        fprintf(foutUnion,textOut); % output to union network
        if totTargTfs > 1 % output to Min2 network
            fprintf(foutMin2,textOut);
        end
        if totTargTfs > clustSizes(clusts2keep(kc))/2 % output to > 50% network
            fprintf(fout50per,textOut);
        end            
    end
    fclose(fout50per); fclose(foutMin2); fclose(foutUnion);
    unix(['wc -l ' fullfile(netOutDir, [currFileBase '*'])]);        
end