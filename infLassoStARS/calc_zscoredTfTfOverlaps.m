function outMat = calc_zscoredTfTfOverlaps(inputNetwork,tfTargMin,...
    targTfMin, fdrCut, edgeOpt, outDir)
%% outMat = calc_zscoredTfTfOverlaps(inputNetwork,tfTargMin,...
%%    targTfMin, fdrCut, edgeOpt, outDir)
% Taking inspiration from the context likelihood of relatedness (CLR) 
% (normalized mutual information, proposed by Faith et al. PLoS Comp. Bio.
% 2007), this script caculates the z-scored overlap between TFs' regulatory 
% interactions, resulting in a "normalized overlap" between TFs:
%               Z(i,j)  = sqrt(Z(i)^2 + Z(j)^2), for Z(i), Z(j) > 0
%                       = -sqrt(Z(i)^2 + Z(j)^2), for Z(i), Z(j) < 0 **
%                       = 0, otherwise
% where Z(i) is the z-score of overlap for TF i with j, based on the 
% distribution of TF i's overlap with other TFs, as described in
%% Reference: Miraldi et al. (2019) Genome Research. [Equation 7]
%% ** with the modification indicated by "**" above
%% INPUTS:
% inputNetwork -- a tab-delimited "table" formatted network file, where
%   columns correspond to TF regulators, rows correspond to target genes,
%   and values correspond to signed interactions
% tfTargMin -- limit analysis to TFs with at least this number of targets
%   (applied first)
% targTfMin -- limit analysis to targets with at least this number of TFs
%   (applied second)
% fdrCut -- cutoff for TF pair inclusion (applied third)
% edgeOpt -- can be set to one of three options:
%   'pos' -- limits analysis to positive regulatory interactions
%   'neg' -- limits anlaysis to negative regulatory interactions
%   'comb' -- integrates both positive and negative regulatory
%       interactions. For example, if TF A and B both negatively regulate
%       Gene 1, that regulatory interaction contributes to overlap,
%       while, if Gene 2 is positively regulated by TF A and negatively
%       regulated by TF B, that regulatory interaction will not contribute
%       to target overlap between TFs (because sign is different).
% outDir -- output directory for results
%% OUTPUTS:
% outMat -- .mat output containing raw overlaps, p-vals of overlap,
%   z-scored overlap matrix, etc.

%% debugging:
% matlabDir = '..';
% addpath(fullfile(matlabDir,'infLassoStARS'))
% addpath(fullfile(matlabDir,'customMatlabFxns'))
% 
% inputNetwork = '../Th17_example/outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/ATAC_Th17_bias50_maxComb/ATAC_Th17_bias50_maxComb_cut01.tsv';
% tfTargMin = 20;     % only consider TFs with at least this number of targets
% targTfMin = 1;      % only consider targets with at least this number of TFs
% fdrCut = 1;         % cutoff for TF pair inclusion
% edgeOpt = 'comb';   % can be set to one of three options:
% [outDirBase,fileName,ext] = fileparts(inputNetwork);
% outDir = fullfile(outDirBase,fileName,strjoin({'zOverlaps',...    
%     [edgeOpt 'Edge'],...
%     ['fdr' num2str(100*fdrCut)],...
%     ['tfMin' num2str(tfTargMin)],...
%     ['targMin' num2str(targTfMin)]},'_'));

disp(outDir)
mkdir(outDir);

% load network file
fid = fopen(inputNetwork,'r');
tline=fgetl(fid);
tfs = cellstr(strvcat(strsplit(tline,'\t')));
totTfsRaw = length(tfs);  % determine how many columns of TFs file has
fclose(fid);
fid = fopen(inputNetwork,'r'); % get the rest of the data matrix
C = textscan(fid,['%s' repmat('%f',1,totTfsRaw)],'Delimiter','\t','HeaderLines',1);
targs = C{1};
totTargsRaw = length(targs);
regInts = [C{2:end}];
fclose(fid);

% note we're using genes in network as background:
intsTot = length(targs);

tfPairAnal = [];
tfPairAnal.edgeType = edgeOpt;

currRegInts = zeros(totTargsRaw,totTfsRaw);
if length(intersect({edgeOpt},{'pos'}))
    disp('Positive Edge Mode')
    for tfi = 1:totTfsRaw
        currTargs = regInts(:,tfi);
        keepInds = find(currTargs>0);
        currRegInts(keepInds,tfi) = currTargs(keepInds);
    end
elseif length(intersect({edgeOpt},{'neg'}))
    disp('Negative Edge Mode')
    for tfi = 1:totTfsRaw
        currTargs = regInts(:,tfi);
        keepInds = find(currTargs<0);
        currRegInts(keepInds,tfi) = currTargs(keepInds);
    end
elseif length(intersect({edgeOpt},{'comb'}))
    disp('Combined, signed edge mode')
    disp('FDR estimate and statistics are conservative in this')
    disp('mode, as positive and negative regulatory edges to a')
    disp('particular gene are mutually exclusive (not independent).')        
    % because FDR filter below uses total genes as background, while
    % here we track gene and sign, thus overlap is the result of two
    % processes (hypergeometric PDF to select gene and a second random
    % process to select sign, which we don't account for, thus p-values 
    % overestimate probability of overlap due to random chance)
    currPosInts = currRegInts;
    currNegInts = currRegInts;
    for tfi = 1:totTfsRaw
        currTargs = regInts(:,tfi);
        keepInds = find(currTargs<0); % negative
        currNegInts(keepInds,tfi) = currTargs(keepInds);
        keepInds = find(currTargs>0); % positive
        currNegInts(keepInds,tfi) = currTargs(keepInds);
    end                    
    currRegInts = [currPosInts; currNegInts];
    intsTot = 2*intsTot; % double because we could + and - edges to gene separately
    targs = cellstr(strvcat(strvcat(targs),strvcat(targs))); % will double 
    % targets as well (needed for output_topN_tfTfclusters_jp_gene_viz.m)
else
    error('For "edgeOpt" please choose from "pos", "neg", or "comb".')
    
end

% filter the interaction matrix      
targsPerTf = sum(abs(sign(currRegInts)));   % first TFs
keepTfsCut = find(targsPerTf>=tfTargMin);
currTfs = cellstr(strvcat(tfs{keepTfsCut}));
totTfs = length(keepTfsCut);        
currRegInts = currRegInts(:,keepTfsCut);
tfsPerTarg = sum(abs(sign(currRegInts)),2); % then targets 
keepTargs = find(tfsPerTarg>=targTfMin);
currTargs = cellstr(strvcat(targs{keepTargs}));
currRegInts = currRegInts(keepTargs,:);

% calculate the overlaps and signficance for each TF pair
currTotPairs = totTfs*(totTfs-1)/2;
currNames = cell(currTotPairs,1);
olapPvals = ones(currTotPairs,1);
perOverlaps = zeros(currTotPairs,1);
perMaxOverlaps = zeros(currTotPairs,1);
perMinOverlaps = zeros(currTotPairs,1);
overlaps = zeros(currTotPairs,1);
count = 0;
for tf1 = 1:totTfs-1
    tfTargs1 = find(currRegInts(:,tf1));
    totTf1 = length(tfTargs1);
    for tf2 = tf1+1:totTfs
        count = count + 1;
        currNames{count} = strjoin(sort({tfs{tf1},tfs{tf2}}),'--');        
        tfTargs2 = find(currRegInts(:,tf2));
        totTf2 = length(tfTargs2);
        overlaps(count) = length(intersect(tfTargs1,tfTargs2));
        if overlaps(count) > 0
            perOverlaps(count) = overlaps(count)/mean(totTf1,totTf2);
            perMaxOverlaps(count) = overlaps(count)/min(totTf1,totTf2); % upper limit on percent
            perMinOverlaps(count) = overlaps(count)/max(totTf1,totTf2); % lower limit on percent
            currOv = max(0,overlaps(count)-1);
            olapPvals(count) = hygecdf(max(currOv-1,0),intsTot,totTf1,totTf2,'upper'); % raw p-value for overlap                
        end
    end
end            
% calculate adjusted p-values
[adjPvals, boolSig] = bh_adjust_pval(olapPvals,fdrCut);

% store initial results, all overlaps
tfPairAnal.currRegInts = currRegInts;
tfPairAnal.targsPerTf = targsPerTf;
tfPairAnal.currTotPairs = currTotPairs;
tfPairAnal.currNames = currNames;
tfPairAnal.olapPvals = olapPvals;
tfPairAnal.perOverlaps = perOverlaps;
tfPairAnal.perMaxOverlaps = perMaxOverlaps;
tfPairAnal.perMinOverlaps = perMinOverlaps;
tfPairAnal.overlaps = overlaps;

sqPerOverlaps = squareform(perOverlaps);
sqOverlaps = squareform(overlaps);
sqSigs = squareform(-log10(adjPvals));
sqIsSig = squareform(boolSig);
sqRawSigs = squareform(log10(olapPvals));

% limit to signficantly overlapping TFs
keepTfs = find(sum(sqIsSig));       
sigTfs = cellstr(strvcat(currTfs{keepTfs}));
sigPerOverlaps = sqPerOverlaps(keepTfs,keepTfs);
sigOverlaps = sqOverlaps(keepTfs,keepTfs);
sigSigs = sqSigs(keepTfs,keepTfs);
sigIsSig = sqIsSig(keepTfs,keepTfs);   
sigRawSigs = sqRawSigs(keepTfs,keepTfs);
% limit to targets with regulatory edges
sigRegIntsTmp = currRegInts(:,keepTfs);
keepTargs = find(sum(abs(sigRegIntsTmp),2)>0);
sigRegInts = currRegInts(keepTargs,keepTfs);    

tfPairAnal.sigTfs = sigTfs;
tfPairAnal.sigTargs = cellstr(strvcat(currTargs{keepTargs}));
tfPairAnal.sigPerOverlaps = sigPerOverlaps;
tfPairAnal.sigOverlaps = sigOverlaps;
tfPairAnal.sigSigs = sigSigs;
tfPairAnal.sigIsSig = sigIsSig;
tfPairAnal.sigRawSigs = sigRawSigs;
tfPairAnal.sigRegInts = sigRegInts;    

%% Get normalized overlap / pairwise "geometric mean" of z-score
zMat = pairwiseZnormSigned(sigOverlaps);
% define a distance from the similarity matrix zMat
maxZmat = max(abs(zMat(:))); % ensures distances are all > 0
zMatDist = maxZmat - zMat;

tfPairAnal.zMatDist = zMatDist;
tfPairAnal.zMat = zMat;

outMat = fullfile(outDir,'tfPair.mat');
save(outMat,'tfPairAnal')
disp([outMat '.mat created.'])