function [hygePs, overlaps, netSizes]...
        = calcRankListHygePs(inputValues,setList)
%% [aucpr, aroc,precisions, recalls, fprs, uniInVals, f1scores]...
%       = aupr_step_outVals(inputValues,setList)
%% GOALS: 
%       Calculate enrichment p-values (hypergeometric CDF) for a range of
%       model sizes
%% INPUTS:
%       inputValues - [N x 1] vector of real-valued confidences, unsorted,
%           higher value indicates higher confidence
%       setList - [N x 1] vector denoting whether entities belonging to a
%           particular set (e.g., disease-associated genes OR a "gold 
%           standard" set). 
%% OUTPUTS:
%       aucpr - area under the precision recall-curve
%       aroc - are under ROC
%       precisions - [M + 1 X 1] vector of precisions as you walk down the
%           ranked list.  Precision = fraction of retrieved instances that 
%           are relevant. where M is the number of levels in the input
%           value
%           (i.e., | genes @ top of list & in set| / |genes @ top of list|
%       recalls - [M + 1 x 1] vector of recalls as you walk down the ranked
%           list. Recall = fraction of relevant instances retrieved.
%           (i.e., |genes in set & top of list| / |genes in set|)
%       fprs - [M +1 x 1] vector of false positive rates (FPRs)
%           FPR = FPs / All Negatives, FP = false postive, TN = true negative
%       uniInVals -- [M x 1] outputs a vector unique inputValues
%       f1scores -- [M x 1] outputs a vector of f1scores: harmonic mean of
%           precision and recall
%% Reference:
% Miraldi et al. "Leveraging chromatin accessibility data for 
%   transcriptional regulatory network inference in T Helper 17 Cells"
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital

%% make sure that inputValues and setList have the same dimensions
if length(inputValues) ~= length(setList)
    error('inputValues and setLists vectors must have the same length')
end

%% find out how many levels of inputValues there are 
uniInVals = sort(unique(inputValues),'descend');
totLevels = length(uniInVals);
inSet = find(setList);
totInSet = length(inSet);
allEdges = length(inputValues);

hygePs = zeros(totLevels,1);
overlaps = zeros(totLevels,1); 
netSizes = zeros(totLevels,1);
% get p-r curve and ROC curve values
for lev = 1:totLevels
    predInds = find(inputValues>=uniInVals(lev));  % find all prediction at this confidence level
    overlap = length(intersect(predInds,inSet)); % true positives
    totPreds = length(predInds);
    hygePs(lev) = hygecdf(max(0,overlap-1),allEdges,totInSet,totPreds,'upper');
    netSizes(lev) = totPreds;
    overlaps(lev) = overlap;
end

