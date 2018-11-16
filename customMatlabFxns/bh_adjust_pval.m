function [adjpvals,bool_minFDRs] = bh_adjust_pval(pvals,FDRcutoffs)
%% function [adjpvals,bool_minFDRs] = bh_adjust_pval(pvals,FDRcutoffs)
% Generates adjusted p-values, i.e., the FDR at which each test would be
% significant, according to Benjamini Hochenburg correction
% Inputs:
%   pvals - [H X 1] vector of pvalues
%   FDRcutoffs - (optional) [T X 1] vector of FDR cutoffs to test
% Outputs:
%   adjpvals - [H X 1] vector of pvalues adjusted for multiple hypothesis testing
%   bool_minFDRs - (if FDRcutoffs is specified) [H X T] boolean matrix,
%       where entry (h,t) = 1 <==> the adjusted p-value (h,t) is
%       significant at the specified FDR and not significant for a smaller
%       FDR specified in FDRcutoffs.  (i.e., if adjpvals(2,1) = .01, and
%       FDRcutoffs = [.05 .1], then bool_minFDRs(2,1) = 1 and
%       bool_minFDRs(2,2) = 0)

%% debug
% pvals = [0 .05 .25 .97 .99 .98]';

npostests = length(pvals);
tested = find(isfinite(pvals));     % hypotheses tested
ntests = length(tested);            % # "           "

adjpvals = NaN*ones(npostests,1);   % set adjusted p-values to NaN

% rank the hypotheses tested according to p-value
[pordered indordered] = sort(pvals(tested),'ascend');
bhdenom = [ntests:-1:1]';           % denominator for each ranked p-value 
%   for which FDR at a certain level would be compared to

% final adjusted p-values
ps = pordered.*bhdenom;
onez = find(ps >= 1);   % find first "p-value" >= 1
if length(onez) > 1
    ps(onez(1):end) = 1;    % set that p-value and all subsequent p-values to 1
end
adjpvals(tested(indordered)) = ps;  

used = [];
if nargin > 1
    ncutoffs = length(FDRcutoffs);
    bool_minFDRs = zeros(npostests,ncutoffs);
    [cutoffs_ordered,inds_sorted] = sort(FDRcutoffs,'ascend');
    for pcut = 1:ncutoffs
        currcut = cutoffs_ordered(pcut);
        sigs = setdiff(find(adjpvals<currcut),used);
        bool_minFDRs(sigs,inds_sorted(pcut)) = 1;
        used = [used; sigs];  % don't count something as significant twice
    end
end
    
