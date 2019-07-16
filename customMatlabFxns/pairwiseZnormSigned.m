function zMat = pairwiseZnormSigned(interactions)
%% function zMat = pairwiseZnormSigne(interactions)
% Taking inspiration from the context likelihood of relatedness (CLR) 
% (normalized mutual information, proposed by Faith et al. PLoS Comp. Bio.
% 2007), this script caculates the z-scored overlap between TFs' regulatory 
% interactions, resulting in a "normalized overlap" between TFs:
%               Z(i,j)  = sqrt(Z(i)^2 + Z(j)^2), for Z(i), Z(j) > 0
%                       = -sqrt(Z(i)^2 + Z(j)^2), for Z(i), Z(j) < 0 **
%                       = 0, otherwise
% where Z(i) is the z-score of overlap for TF i with j, based on the 
% distribution of TF i's overlap with other TFs, as described in:
%% Reference: Miraldi et al. (2019) Genome Research. [Equation 7]
%% ** with the modification indicated by "**" above
%% INPUTS: 
% interactions -- TFs X targets matrix, where nonzero value
%   indicates a regulatory interaction
%% OUTPUTs: 
% zOvMat -- TF X TF matrix of z-scored overlaps
%% Description:
% For each TF i, calculate overlaps with each TF j (j != i) per row to
% yield a distribution of overlaps for TF i  and z-scores (note that
% the self-overlap of TF i (100%) is exluded from the calculation). Store
% those z-scores in row i of TF x TF matrix Z
% To arrive at a symmetric matrix of overlap "similarities" Zmat, for each 
% TF pair i, j, set Zmat(i,j) and Zmat(j,i) equal to the geometric mean of
% Z(i,j) and Z(j,i), when sign is shared. Specifically, 
%               Zmat(i,j)   = sqrt(Z(i,j)^2 + Z(j,i)^2), for Z(i,j), Z(j,i) > 0
%                           = -sqrt(Z(i,j)^2 + Z(j,i)^2), for Z(i,j), Z(j,i) < 0
%                           = 0, otherwise
% where Z(i,j) is the z-score of TF i with j (based on i's overlap)

[rows,cols] = size(interactions);

% 1. get a zscore for each row:
zscores = zeros(rows);
for rind = 1:rows
    currInds = setdiff(1:rows,rind); % don't include self-self values
    zscores(rind,currInds) = zscore(interactions(rind,currInds));
end

% 2. calculate pairwise normalized overlap similarity
zMat = zeros(rows);
for ro = 1:rows-1
    for co = ro+1:rows
        if and(zscores(ro,co)>0,zscores(co,ro)>0)
            zMat(ro,co) = sqrt(zscores(ro,co)^2 + zscores(co,ro)^2);
            zMat(co,ro) = zMat(ro,co);
        elseif and(zscores(ro,co)<0,zscores(co,ro)<0) 
            zMat(ro,co) = -sqrt(zscores(ro,co)^2 + zscores(co,ro)^2);
            zMat(co,ro) = zMat(ro,co);
        end
    end
end

