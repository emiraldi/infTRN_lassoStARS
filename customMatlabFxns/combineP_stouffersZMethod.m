function pcomb = combineP_stouffersZMethod(pvals,weights)
%% pcomb = combineP_stouffersZMethod(pvals,weights)
% combine P-values using Stouffer's Method AKA Weighted Z method, which
% is symmetric so does not suffer from asymmetry bias like Fisher's Method
% REFERENCE: Whitlock. J. Evol. Bio. 2005.
% METHOD: convert the "k" input p-values to z-scores (k>=1), using the normal
% distribution, calculate Zs = weighted average of z-scores.  
% ASSUMES: one-side p-values, divide answer by two for two-sided test
% INPUTS:
%   pvals -- k X 1 vector of p-values
%   weights -- [OPTIONAL] k X 1 vector of weights (could correspond to
%       study size), default sets all weights to one
% OUTPUT: 
%   pcomb -- a combined p-value

% debugging:
% pvals = [.001, .01, .1]';
% weights = ones(size(pvals));

k = length(pvals);
if nargin < 2 % default: treat all methods equally
    weights = ones(k,1);
end

Zs = norminv(pvals,0,1);
Zave = (weights'*Zs / sqrt(sum(weights)));
pcomb = normcdf(Zave,0,1);
