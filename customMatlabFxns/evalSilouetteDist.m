function [silScores, meanSilScore, stderrSilScore] = evalSilouetteDist(pairDists,clustSoln)
%% function [silScores, meanSilScore, stderrSilScore] = evalSilouetteDist(pairDists,clustSoln)
%% Given a symmetric matrix of pairwise distances and clustering solution, 
% this function estimates silouette scores for all clustered objects
%% Silhouette Value (definition from Mathworks website)
% The silhouette value for each point is a measure of how similar that point 
% is to points in its own cluster, when compared to points in other clusters. 
% The silhouette value for the ith point, Si, is defined as:
% Si = (bi-ai)/ max(ai,bi)
% where ai is the average distance from the ith point to the other points 
% in the same cluster as i, and bi is the minimum average distance from the 
% ith point to points in a different cluster, minimized over clusters.
% The silhouette value ranges from -1 to +1. A high silhouette value 
% indicates that i is well-matched to its own cluster, and poorly-matched 
% to neighboring clusters. If most points have a high silhouette value, 
% then the clustering solution is appropriate. If many points have a low 
% or negative silhouette value, then the clustering solution may have 
% either too many or too few clusters. The silhouette clustering evaluation 
% criterion can be used with any distance metric.
%% INPUTS: 
%   pairDists -- symmetrix matrix of distances between clustered
%       objects (dimmensions: obs X obs)
%   clustSoln -- vector of indicating which object belongs to which cluster
%       (obs X 1), 
%% OUTPUTS: 
%   silScores -- vector of each object's silhouette score (obs X 1)
%   meanSilScore -- mean silhouette score
%   stderrSilScore -- standard error of silhouette score

% clustSoln = currSoln;
% pairDists = zMatDist;

obs = length(clustSoln);
uniClusts = unique(clustSoln);
totClusts = length(uniClusts);

silScores = zeros(obs,1);%NaN*ones(obs,1);

% make sure that distances along the diagonal are zero
pairDists = pairDists - eye(obs).*pairDists;

for cind = 1:totClusts
    currClust = uniClusts(cind);
    clustInds = find(clustSoln==currClust);    
    clustSize = length(clustInds);
    if clustSize > 1
        othInds = setdiff(1:obs,clustInds);
        % note we calculate the average dividing by (clustSize - 1), because we
        % don't want to count self-distance between the same object, and the 
        aveWithinDist = sum(pairDists(clustInds,clustInds),2)/(clustSize-1);  
        othClustDists = zeros(clustSize,totClusts-1);
        othClusts = setdiff(1:totClusts,currClust);
        for oi = 1:totClusts-1
            currOthClust = find(clustSoln==othClusts(oi));
            othClustDists(:,oi) = mean(pairDists(clustInds,currOthClust),2);
        end
        minOutsideDist = min(othClustDists,[],2);

        silScores(clustInds) = (minOutsideDist - aveWithinDist) ./ ...
            max(minOutsideDist,aveWithinDist);
    end
end

realSilScores = silScores(find(isfinite(silScores)));
meanSilScore = mean(realSilScores);
stderrSilScore = std(realSilScores)/sqrt(length(realSilScores)-1);
