function [aucpr, aroc,...
    precisions, recalls, fprs, uniInVals, f1scores]...
        = aupr_step_outVals(inputValues,setList)
%% [aucpr, aroc,precisions, recalls, fprs, uniInVals, f1scores]...
%       = aupr_step_outVals(inputValues,setList)
%% GOALS: 
%       Calculate precision-recall, FPR, TPR for P-R and ROC curves, AUPR,
%       AROC, and f1-scores
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

precisions = zeros(totLevels,1);
recalls = zeros(totLevels,1); 
fprs = zeros(totLevels,1);
tns = length(setList) - totInSet;
% get p-r curve and ROC curve values
for lev = 1:totLevels
    predInds = find(inputValues>=uniInVals(lev));  % find all prediction at this confidence level
    hits = length(intersect(predInds,inSet)); % true positives
    totPreds = length(predInds);
    recalls(lev) = hits/totInSet; 
    precisions(lev) = hits/totPreds;
    fps = totPreds-hits;    % false positives
    fprs(lev) = fps/tns;
end
nonzeroPrecisions = find(precisions);
if length(nonzeroPrecisions)
    recalls = recalls(nonzeroPrecisions);
    precisions = precisions(nonzeroPrecisions);
    fprs = fprs(nonzeroPrecisions);
end
f1scores = 2*precisions.*recalls./(precisions + recalls);
recalls = [0; recalls];
precisions = [precisions(1); precisions];
fprs = [0; fprs];
% calculate aupr and aroc
% p-r
heights = (precisions(2:end)+precisions(1:end-1))/2;
widths = recalls(2:end) - recalls(1:end-1);
aucpr = heights'*widths;
% roc
heights = (recalls(2:end) + recalls(1:end-1))/2;
widths = fprs(2:end) - fprs(1:end-1);
aroc = heights'*widths;
