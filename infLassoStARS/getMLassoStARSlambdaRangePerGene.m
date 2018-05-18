function [minLambdas, maxLambdas, maxedOut, notSmallEnough,...
    minLambdaNet,maxLambdaNet,maxOutNet,minOutNet,...
    netInstabilitiesLb,netInstabilitiesUb,...
    instabilitiesLb,instabilitiesUb] = ...
  getMLassoStARSlambdaRangePerGene(predictorMat,...
    responseMat,priorWeightsMat,lambdaRange,...
    targetMinInstability,targetMaxInstability,subsampleSize, totSS)
%% function [minLambdas, maxLambdas, maxedOut, notSmallEnough,...
%     minLambdaNet,maxLambdaNet,maxOutNet,minOutNet] = ...
%   getMLassoStARSlambdaRangePerGene(predictorMat,...
%     responseMat,priorWeightsMat,lambdaMax,lambdaMin,logLambdaStep,...
%     targetMinInstability,subsampleSize, totSS)
%% GOAL: provide both per-gene and network-wide edge average 
% instabilities for a given range of lambda penalties for LASSO-StARS
% using bStARS average instability bounds to speed calculation:
% -- bootstrap = 2 instability calculation provides a lower bound lambda
% -- analytical upper bound derived from Poisson Binomial Distribution
% Reference:
% Muller, Kurtz, Bonneau. "Generalized Stability Approach for Regularized
%   Graphical Models". 23 May 2016. arXiv.
% Average instability is defined as in Liu, Roeder, Wasserman,
%   "Stability Approach to Regularization Selection (StARS) for High
%   Dimensional Graphical Models". 16 June 2010. arXiv.  Specifically
%   Ave. Instability (Liue et al.) = 1/2 * Ave. Instability (Muller et al.)
% We find a range of lambda penalties for each response, separately, using 
% LASSO regression (variables defined below in inputs):
% find B* that min B ( ||response - B * predictorMat|| +...
%   sum_ij | priorWeightsMat_ij * B_ij | )
% This version solves the LASSO optimization problem with Glmnet, Reference:
% Glmnet for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, 
% R. and Simon, N. -- http://www.stanford.edu/~hastie/glmnet_matlab/
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Main Reference: Miraldi et al. "Leveraging ATAC-seq data for 
%   transcriptional regulatory network inference in T Helper 17 Cells"
%% INPUTS:
%   predictorMat -- predictors X samples matrix
%   responseMat -- response X samples matrix
%   priorWeightsMat -- response X predictors matrix, showing what lambda
%       weights will be multiplied by -- note: can contain infinite weights
%       if TF-gene interaction should be filtered (e.g., based on CLR)
%   lambdaMax -- max lambda value to consider, [0,1]
%   lambdaMin -- min lambda value to consider, [0,1]
%   logLambdaStep -- number of lambda values to sample for each order
%       of magnitude (e.g., logLambdaStep = n, we'd get every 1/n log10
%       step) 10 has worked well in practice
%   targetMinInstability -- scalar (0,.5) for target minimum instability
%   targetMaxInstability -- scalar (0,.5) for target maximum instability
%   subsampleSize -- { 1,2,... samples-1 }, number of samples per subsample
%       The references above recomend subsampleSize =
%       floor(10*sqrt(samples)), but this is only feasible for >= 100
%       samples
%   totSS -- number of subsamples, 2 might be sufficient for network-wide
%       lambda bounds, more might be required for reproducible
%       gene-specific lambda bounds
%   NOTE: predictorMat and responseMat will be z-scored predictor and
%       response-wise, predictorMat z-scoring will be done for each gene
%       (as only a subset of predictors are expected per gene)
%% OUTPUTS:
%   minLambdas -- [response X 1] corresponds to the lambda at which 
%        targetMaxInstability value is reached for each gene
%   maxLambdas -- [response X 1] corresponding to max lambda for 
%       targetMinInstability (supplied as input) for ALL responses (so some 
%       response might reach the desired target instability for a smaller 
%       lambda value)
%   maxedOut -- vector of response indices for which lambda upper bound was
%       the lambda range max, suggesting that function should be rerun with
%       a larger "lambdaMax" parameter to ensure that targetMinInstability is
%       reached
%   notSmallEnough -- vector of response indices for which the lambda lower
%       bound might have been too small (as maximum instability was reach at
%       the boundary, lambdaMin).  Try running the function with a smaller
%       "lambdaMin" so that output is guaranteed to
%       included the lambda corresponding to maximum instability
%   minLambdaNet -- scalar, network average corresponding to maximum
%       instability lambda value
%   maxLambdaNet -- scalar, network average instability corresponding to
%       targetMinInstability
%   maxOutNet -- 1 --> maxLambda corresponds to input "lambdaMax",
%       suggesting that test range should be extended higher, maxLambdaNet might
%       be an underestimate, 0 --> otherwise
%   minOutNet -- 1 --> minLambda correspond to input "lambdaMin",
%       suggesting that the minLambdaNet might be an overestimate and test
%       range should be extended lower
%   <<Figure>> Will also output a figure(100) with LASSO solution paths
%       every 500 gene models to ward against wonkiness
%   <<Figure>> 5 subplots: (1) instability cutoff minimum lambdas per gene, 
%       (2) max instability lambdas per gene, (3) Instabilities Upper Bound
%       per gene, (4) Instabilities Lower Bound per gene, (5) Network-wide
%       Instabilities U.B. and L.B., with instability cutoff lambda and max
%       lambda marked.

%% debugging inputs:
% predictorMat = predictorMat(:,trainInds);
% responseMat = responseMat(:,trainInds);
% save('StARSfood.mat',...
% 'predictorMat',...
% 'responseMat',...
% 'priorWeightsMat',...
% 'lambdaMax',...
% 'lambdaMinMaxRatio','logLambdaStep','targetMinInstability',...
% 'subsampleSize',...
% 'totSS')
% load StARSfood.mat
% lambdaMax = 10;
% lambdaMinMaxRatio = 1E-3;
% figure(1), clf
% load StARSfoodVector.mat
% addpath('~/erm/MATLAB/glmnet')
% addpath('~/erm/MATLAB/emily_functions')    
% totSS = 4;

%% END Debugging Inputs

totLambdas = length(lambdaRange);
options = glmnetSet;
% For glmnet, we have to provide lambda's in decreasing order and we'll
% have to flip back the output from glmnet. We don't want to use default
% glmnet lambda, as we have to keep the same lambda for each gene-model
% regression
options.lambda = fliplr(lambdaRange);

%% calculate instabilities with two subsamples
fontSize = 12;
[totResponses,totSamps] = size(responseMat);
totPreds = size(predictorMat,1);

% ssMatrix = zeros(totLambdas,totResponses,totPreds);
instabilitiesLb = zeros(totResponses,totLambdas);
instabilitiesUb = zeros(totResponses,totLambdas);
minLambdas = zeros(totResponses,1);
maxLambdas = zeros(totResponses,1);

% use soft-thresholding to solve the LASSO problem

% get subsamp indices
subsamps = zeros(totSS,subsampleSize);
for ss = 1:totSS
    subsamps(ss,:) = randperm(totSamps,subsampleSize);
end

% get (finite) predictor indices for each response
responsePredInds = cell(totResponses,1);
for res = 1:totResponses
    currWeights = priorWeightsMat(res,:);
    responsePredInds{res} = find(isfinite(currWeights));
end    

tic
% network-level instabilities, will be calculated as a weighted average of
% gene instabilities, so that each edge has equal weight
netInstabilitiesUb = zeros(totLambdas,1); % 
netInstabilitiesLb = zeros(totLambdas,1);
totEdges = 0;   % denominator for network Instabilities
% now we're building each model genewise
% totResponses = 50;
for res = 1:totResponses % can be a parfor loop
    % limit to predictors with finite edges
    predInds = responsePredInds{res};
    currPredNum = length(predInds);
    currOptions = options;
    currOptions.penalty_factor = priorWeightsMat(res,predInds)';
    totEdges = totEdges + currPredNum;
    ssVals = zeros(totLambdas,currPredNum);
    for ss = 1:totSS
        subsamp = subsamps(ss,:);
        currPreds = zscore(predictorMat(predInds,subsamp)'); % nobs X nvars
        currResponses = zscore(responseMat(res,subsamp)'); % nobs X 1
        lsoln = glmnet(currPreds,currResponses,'gaussian',currOptions);        
        % lsoln.beta == predictors X lambda range, coefficient matrix
        currBetas = fliplr(lsoln.beta); % flip so that the lambdas are increasing
            % abs(sign()) as we only want to track nonzero edge occurrences
        ssVals = ssVals + abs(sign(currBetas))'; 
    end
    % calculate instabilities for the gene
    theta2 = (1/totSS)*ssVals; % empirical edge probability
    instabilitiesLb(res,:) = 2*(1/currPredNum)*sum(theta2.*(1-theta2),2); % bStARS lower bound
    netInstabilitiesLb = netInstabilitiesLb + currPredNum*instabilitiesLb(res,:)'; % weighted sum
    theta2mean = sum(theta2,2)./currPredNum; 
    instabilitiesUb(res,:) = 2*theta2mean.*(1-theta2mean); % bStARS upper bound
    netInstabilitiesUb = netInstabilitiesUb + currPredNum*instabilitiesUb(res,:)'; % weighted sum
    if rem(res,500) == 0
        % make sure nothing wonky is happening with lasso solution paths
        disp([num2str(res) ' models built.'])
%         figure(100), clf
%         glmnetPlot(lsoln)      
%         shg
    end
end
for res = 1:totResponses
    % take the supremum, find max Lambda, and set all smaller lambdas equal
    % to that value 
    [maxLb, maxLbInd] = max(instabilitiesLb(res,:));
    instabilitiesLb(res,1:maxLbInd(end)) = maxLb;
    [maxUb, maxUbInd] = max(instabilitiesUb(res,:));
    instabilitiesUb(res,1:maxUbInd(end)) = maxUb;
    % find the minimum lambda for the gene, based on maximum for upper bound
    % we are less interested in high instability lambdas, so okay to use
    % upper bound
    [xx, maxInstInd] = min(abs(instabilitiesLb(res,:)-targetMaxInstability));
    minLambdas(res) = lambdaRange(maxInstInd(end)); % to the right
    % find the lambda nearest the min instability worth considering, use
    % upper bound as that will be sure to find an lambda >= target instability lambda 
    [xx, minInstInd] = min(abs(instabilitiesUb(res,:)-targetMinInstability));
    maxLambdas(res) = lambdaRange(minInstInd(end));  % to the right
    % note for typical bStARS, where you know what instability cutoff you
    % want you'd use the upperbound to find the min lambda and the lb to
    % find the max lambda    
end
toc

% get network-level instabilities
netInstabilitiesUb = netInstabilitiesUb/totEdges;
netInstabilitiesLb = netInstabilitiesLb/totEdges;
[maxLb, maxLbInd] = max(netInstabilitiesLb);
netInstabilitiesLb(1:maxLbInd(end)) = maxLb; % take supremum for lambdas smaller than instability max
[maxUb, maxUbInd] = max(netInstabilitiesUb);
netInstabilitiesUb(1:maxUbInd(end)) = maxUb; 
[xx, maxInstInd] = min(abs(netInstabilitiesLb - targetMaxInstability));
minLambdaNet = lambdaRange(maxInstInd(end));
[xx, minInstInd] = min(abs(netInstabilitiesUb - targetMinInstability));
maxLambdaNet = lambdaRange(minInstInd(end));

% find out if maxLambda or minLambda were at the lambda range extremes
maxedOut = find(maxLambdas == lambdaRange(end));
notSmallEnough = find(minLambdas == lambdaRange(1));

if maxedOut
    disp(['Warning: ' num2str(length(maxedOut)) ' gene(s) might not have hit target min instability.'])
    disp(['If using gene-level instabilities, re-run with larger "lambdaMax" parameter. Examine output figures.'])
end
if notSmallEnough
    disp(['Warning: ' num2str(length(notSmallEnough)) ' gene(s) might not have hit max instability.'])
    disp(['If usinging gene-level instabilities, re-run with a smaller "lambdaMin" parameter.  Examine output figures.'])
end

maxOutNet = (maxLambdaNet == lambdaRange(end));
minOutNet = (minLambdaNet == lambdaRange(1));
if maxOutNet
    disp(['Warning: Network lambda range might not have hit target instability.'])
    disp(['If using network-level instabilities, re-run with larger "lambdaMax" parameter. Examine output figures.'])
end
if minOutNet
    disp(['Warning: Network lambda range might not have hit max instability.'])
    disp(['If using network-level instabilities, re-run with a smaller "lambdaMin" parameter.  Examine output figures.'])
end

% plot lambda mins
subplot(5,1,1)
hist(minLambdas,lambdaRange)
grid on
xlabel('Minimum Lambda','FontSize',fontSize)
ylabel('Genes','FontSize',fontSize)
set(gca,'xscale','log','FontSize',fontSize)
axis tight
ax = axis();
hold on
% plot network min lambda
plot(minLambdaNet*[1 1], ax(3:4),'r','LineWidth',2)

% plot lambda maxes
subplot(5,1,2)
hist(maxLambdas,lambdaRange)
grid on
xlabel('Maximum Lambda','FontSize',fontSize)
ylabel('Genes','FontSize',fontSize)
set(gca,'xscale','log','FontSize',fontSize)
axis tight
ax = axis();
hold on
% plot network min lambda
plot(maxLambdaNet*[1 1], ax(3:4),'r','LineWidth',2)

subplot(5,1,4)
plot(repmat(lambdaRange,totResponses,1)',instabilitiesLb')
set(gca,'xscale','log','FontSize',fontSize)
grid on
xlabel('Lambda Range','FontSize',fontSize)
ylabel('Instabilities (L.B.)','FontSize',fontSize)
axis([ax(1:2) 0 .5])
subplot(5,1,3)
plot(repmat(lambdaRange,totResponses,1)',instabilitiesUb')
set(gca,'xscale','log','FontSize',fontSize)
axis([ax(1:2) 0 .5])
grid on
xlabel('Lambda Range','FontSize',fontSize)
ylabel('Instabilities (U.B.)','FontSize',fontSize)

subplot(5,1,5)
plot(lambdaRange,netInstabilitiesUb,'bo-','LineWidth',2)
hold on
plot(lambdaRange,netInstabilitiesLb,'ro-','LineWidth',2)
legend({'U.B.','L.B.'},'Location','East')
grid on, grid minor
set(gca,'xscale','log','FontSize',fontSize)
plot(minLambdaNet*[1 1], [0 .5],'k','LineWidth',1)
plot(maxLambdaNet*[1 1], [0 .5],'k','LineWidth',1)
axis([ax(1:2) 0 .5])
xlabel('Lambda Range','FontSize',fontSize)
ylabel('Instabilities','FontSize',fontSize)
