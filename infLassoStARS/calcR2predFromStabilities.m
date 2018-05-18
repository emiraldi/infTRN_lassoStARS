function calcR2predFromStabilities(instabMat,stabilitiesMat,r2OutMat,modSizes)
% calcR2predFromStabilities(instabMat,stabilitiesMat,r2OutMat,modSizes)
%% BASED on stabilities from getMLassoStARSnetworksQuantilePCorr_instabCut
% Does not work for kfoldCv > 1, only intended for single LO analyses
%% GOAL: Calculate R^2_pred, R^2_fit, SSE_pred, SSE_fit for each 
%   instability cutoff, also calculate other statistics, such as average
%   number of TFs per target gene (SSE = Sum of Square Errors)
%% INPUTS:
% instabMat -- .mat file of instabilities and related data from
%   estimateInstabilitiesMLassoStars
% stabilitiesMat -- .mat output from
%   getMLassoStARSnetworksQuantilePCorr_instabCut, where edges are ranked
%   based on an instability cutoff (gene-wise or network-wise)
% r2OutMat -- a path and name for R^2_pred output from this file
% modSizes -- a vector of target average model sizes (TFs / gene), note
%   this function will automatically remove model sizes that would exceed
%   the total number of nominally stable edges in the stabilitiesMat
%% OUTPUTS:
% Files output to this directory:
% .mat: ${priorName}_instRang_63biasX_SSX_CVK_ZOOM_${noTfaOpt}_r2.mat
%   -- contains SSE's, R^2's, # parameters per model,...
% .fig+.pdf: ${priorName}_instRang_63biasX_SSX_CVK_ZOOM_${noTfaOpt}_r2_v_preds
%   -- plot of R2's vs. instability cutoff
%   -- plot of R2's vs. model size
%   Both can be used to help determine an appropriate stability cutoff
%% Debugging:
% clear all, close all
% % save('r2predMats.mat','instabMat','stabilitiesMat','r2OutMat')
% load('~/erm/MATLAB/infLassoStARS/r2predMats.mat')
% modSizes = [.25 .5 .75 1:25 30:5:40]; % desired model sizes

%% load gene expression

load([instabMat '.mat'])
load([stabilitiesMat '.mat'])

[xx, instabTitleInf, ext] = fileparts(instabMat);
instabTitleInf = strrep(strrep(strrep(instabTitleInf,'instRange_63',''),'ZOOM',''),'_',' ');


totPredTargs = length(targGenes);
totConds = length(conditionsc);
totPreds = length(allPredictors);
% modSizes = [.25 .5 .75 1:25]; % desired model sizes
totModSizes = length(modSizes);
instabRangeNet = zeros(totModSizes,1);

[stabOrd] = sort(abs(allStabsMergedTFs(find(allStabsMergedTFs))),'descend'); % make sure cutoffs are in order
totInts = length(stabOrd);

subsampleCutoffs = modSizes; % get actual model cutoffs
rmInds = [];
for ms = 1:totModSizes
    currNetSize = round(totPredTargs*modSizes(ms));
    if currNetSize <= totInts
        subsampleCutoffs(ms) = stabOrd(currNetSize);
        currTheta = round(subsampleCutoffs(ms))/totSS;
        instabRangeNet(ms) = 2*currTheta*(1-currTheta); % get instability at the cutoff
    else
        rmInds = [rmInds; ms];
        disp([num2str(modSizes(ms)) ' TFs/gene cutoff yields ' num2str(currNetSize) ' network size > ' num2str(totInts) ' TF-gene interactions in network'])        
    end    
end
% in the event that a target TF/gene model size is larger than # of
% interactions
% keepModSizes(rmInds) = [];%find(instabRangeNet);
modSizes(rmInds) = [];% = modSizes(keepModSizes);
totModSizes = length(modSizes);
instabRangeNet(rmInds) = [];% = instabRangeNet(keepModSizes);
subsampleCutoffs(rmInds) = [];% = subsampleCutoffs(keepModSizes);

% output matrices
kfoldCv = 1; % this is legacy -- but would make implementing for k-fold CV easier, so leaving it in
SSE_predMat = zeros(totModSizes,totConds,totPredTargs);
SSE_predCvBreak = zeros(totModSizes,kfoldCv);
SSE_predMeanCvBreak = zeros(kfoldCv,1);
SSE_fitMat = zeros(totModSizes,totConds,totPredTargs);
SSE_fitMatMean = zeros(totConds,totPredTargs);
SSE_predMatMean = zeros(totConds,totPredTargs);
r2_pred_byTarg = zeros(totModSizes,totPredTargs);
r2_fit_byTargs = zeros(totModSizes,totPredTargs);
r2_pred_byCond = zeros(totModSizes,totConds);
r2_fit_byCond = zeros(totModSizes,totConds);
SSE_pred = zeros(totModSizes,totConds);
SSE_fit_mean = zeros(totModSizes,totConds);
r2_pred = zeros(totModSizes,1);
r2_fit = zeros(totModSizes,1);
r2_predNz = zeros(totModSizes,1);
r2_fitNz = zeros(totModSizes,1);
SSE_predMatNz = zeros(totModSizes,1);
SSE_predMatMeanNz = zeros(totModSizes,1);
SSE_fitMatNz = zeros(totModSizes,1);
SSE_fitMatMeanNz = zeros(totModSizes,1);
SSE_predMatNzCV = zeros(kfoldCv,totModSizes);
SSE_predMatMeanNzCV = zeros(kfoldCv,totModSizes);
SSE_fitMatNzCV = zeros(kfoldCv,totModSizes);
SSE_fitMatMeanNzCV = zeros(kfoldCv,totModSizes);
SSE_pred_kfoldInstabTargs = zeros(kfoldCv,totModSizes,totPredTargs);
SSE_predMean_kfoldInstabTargs = zeros(kfoldCv,totModSizes,totPredTargs);
r2_pred_kfoldInstabTargs = zeros(kfoldCv,totModSizes,totPredTargs);

totalModelEdges = zeros(kfoldCv,totModSizes);
priorOverlaps = zeros(kfoldCv,totModSizes);
totNzModels = zeros(kfoldCv,totModSizes);
r2_fitNzCV = zeros(kfoldCv,totModSizes);
r2_predNzCV = zeros(kfoldCv,totModSizes);
numParamsMat = zeros(totModSizes,kfoldCv,totPredTargs);
numParams = zeros(totModSizes,totPredTargs);

kcount = 0;

for kind = 1:max(kfoldCv,1)
    %% load CV samples
%     outBase = fullfile(outFolder,cvExtraFolder,priorName);
%     figInf = [outBase '_instRang_' num2str(round(100*subsampleFrac)) ...
%        'bias' num2str(round(100*lambdaBias)) '_SS' num2str(totalSS) outCV '_ZOOM' noTfaOpt];
%     load(figInf)
%     instabOut = fullfile(LOfolder,...
%         [priorName '_instRang_' num2str(round(100*subsampleFrac)) ...
%        'bias' num2str(round(100*lambdaBias)) '_SS' num2str(totalSS) outCV noTfaOpt]);
%     disp(instabTitleInf)
    
    % get the train and test data matrices for the CV subset
    % as well as the train and test TFAs
    totTrainConds = length(trainInds);
    testInds = setdiff(1:totConds,trainInds);
    totTestConds = length(testInds);
    targGenesTrain = responseMat(:,trainInds);       
    targGenesTest = responseMat(:,testInds);
    tfaTrain = predictorMat(:,trainInds);
    tfaTest = predictorMat(:,testInds);
    % mean-center training and test target genes, according to training mean
    meanTargGene = mean(targGenesTrain,2);        
    cTargGenesTrain = targGenesTrain - repmat(meanTargGene,1,totTrainConds);
    cTargGenesTest = targGenesTest - repmat(meanTargGene,1,totTestConds);
    % z-score TFA, according to training mean and std
    meanTfa = mean(tfaTrain,2);
    stdTfa = std(tfaTrain')';
    zTfaTrain = (tfaTrain - repmat(meanTfa,1,totTrainConds))./repmat(stdTfa,1,totTrainConds);
    zTfaTest = (tfaTest - repmat(meanTfa,1,totTestConds))./repmat(stdTfa,1,totTestConds);        

    % Most basic model is that the mean gene expression from the training
    % data is the best predictor of test condition gene expression
    % Because we already mean-centered the test data according to the
    % training data, we can represent the base model as a matrix of
    % zeros 
    baseModelTest = cTargGenesTest;
    baseModelTrain = cTargGenesTrain;
    % SSE of train data mean for test data (we'll calculate R^2
    % = 1 - [ sum(SSE_predMat) / sum(SSEpredMatMean) ]
    SSE_predMatMean(testInds,:) = baseModelTest'.^2;
    SSE_predMeanCvBreak(kind) = sum(sum(baseModelTest'.^2));
    % SSE of train data mean for train data
    SSE_fitMatMean(trainInds,:) = SSE_fitMatMean(trainInds,:) + ...
        baseModelTrain'.^2;        
        
    currSS = zeros(totPreds,1);
    for targ = 1:totPredTargs
        % start with least stringent to most stringent instability cutoff
        currTargName = targGenes{targ};
        currTargValsTrain = cTargGenesTrain(targ,:)';
        currTargValsTest = cTargGenesTest(targ,:)';
        
        currSubsamples = allStabsMergedTFs(targ,:);
        
        for iind = 1:totModSizes
            currCut = subsampleCutoffs(iind);
            
            intInds = find(currSubsamples>=currCut);
            % get target expression levels
            totParams = length(intInds);
            numParamsMat(iind,kind,targ) = totParams;
            numParams(iind,targ) = numParams(iind,targ) + totParams;
            if totParams > 0 % make sure there's at least one feature for regression
                currPredValsTrain = zTfaTrain(intInds,:)';
                currPredValsTest = zTfaTest(intInds,:)';
                % fit a linear model with least squares (use pseudo
                % inverse, just in case we end up in the under-determined
                % regime)
                coefs = currPredValsTrain \ currTargValsTrain;
                trainPreds = currPredValsTrain * coefs;
                testPreds = currPredValsTest * coefs;
                % get SSE for non-zero models:
                SSE_predMatNz(iind) = SSE_predMatNz(iind) + ...
                    sum((testPreds'-currTargValsTest').^2);
                SSE_predMatMeanNz(iind) = SSE_predMatMeanNz(iind) + ...
                    sum(SSE_predMatMean(testInds,targ));
                SSE_fitMatNz(iind) = SSE_fitMatNz(iind) + ...
                    sum((trainPreds'-currTargValsTrain').^2);
                SSE_fitMatMeanNz(iind) = SSE_fitMatMeanNz(iind) + ...
                    sum(SSE_fitMatMean(trainInds,targ));
                SSE_predMatNzCV(kind,iind) = sum((testPreds'-currTargValsTest').^2) + ...
                    SSE_predMatNzCV(kind,iind);
                SSE_predMatMeanNzCV(kind,iind) = sum(SSE_predMatMean(testInds,targ)) + ...
                    SSE_predMatMeanNzCV(kind,iind);
                SSE_fitMatNzCV(kind,iind) = sum((trainPreds'-currTargValsTrain').^2) + ...
                    SSE_fitMatNzCV(kind,iind);
                SSE_fitMatMeanNzCV(kind,iind) = sum(SSE_fitMatMean(trainInds,targ)) + ...
                    SSE_fitMatMeanNzCV(kind,iind);
                SSE_pred_kfoldInstabTargs(kind,iind,targ) = sum((testPreds'-currTargValsTest').^2);
                SSE_predMean_kfoldInstabTargs(kind,iind,targ) = sum(SSE_predMatMean(testInds,targ));
            else
                trainPreds = zeros(totTrainConds,1);
                testPreds = zeros(totTestConds,1);
            end
            % SSE between model vs. train data mean
            SSE_predMat(iind,testInds,targ) = (testPreds'-currTargValsTest').^2;
            % note there will be "kfoldCV-1" estimates of model fit,
            % we'll divide by this denominator later.  For now, sum:
            SSE_fitMat(iind,trainInds,targ) = SSE_fitMat(iind,trainInds,targ) + ...
                (trainPreds'-currTargValsTrain').^2;
            SSE_predCvBreak(iind,kind) = sum((testPreds'-currTargValsTest').^2) +...
                SSE_predCvBreak(iind,kind);            
        end        
    end    
    kcount = kcount + 1;          
    disp([' CV = ' num2str(kcount)])
    clear ssMatrix % I had some issues where ssMatrix was not saved (because it was too big)
    % and the ssMatrix from the past iteration was used if the new .mat
    % file didn't have an ssMatrix
    % clearing ssMatrix here will ensure that next iteration will be based on its own ssMatrix 
    % and cause an error if ssMatrix not found
end
numParams = numParams/kcount; % take the average
if kcount % there were Inferelator results
    %% get output base name
%     mkdir(figOutBase)
    paramNumsMean = zeros(totModSizes,totPredTargs);
    paramNumsStd = zeros(totModSizes, totPredTargs);
    paramNumsMax = zeros(totModSizes, totPredTargs);
    paramNumsMin = zeros(totModSizes, totPredTargs);  

    % R^2 for non-zero models
    r2_predNz = 1 - SSE_predMatNz./SSE_predMatMeanNz;
    r2_fitNz = 1 - SSE_fitMatNz./SSE_fitMatMeanNz;
    r2_predNzCV = 1 - SSE_predMatNzCV./SSE_predMatMeanNzCV;
    r2_fitNzCV = 1 - SSE_fitMatNzCV./SSE_fitMatMeanNzCV;
    r2_pred_kfoldInstabTargs = 1 - SSE_pred_kfoldInstabTargs./ ...
        SSE_predMean_kfoldInstabTargs;
        
    for bind = 1:totModSizes
        % calculate R^2_pred and R^2_fit
        % calculate SSE according to condition (sum over targets)
        currSSE_preds = zeros(totConds,totPredTargs);
        currSSE_preds(:,:) = SSE_predMat(bind,:,:);
        r2_pred(bind) = 1 - (sum(currSSE_preds(:))/sum(SSE_predMatMean(:)));  
        r2_pred_byTarg(bind,:) = 1 - (sum(currSSE_preds)./sum(SSE_predMatMean))';
        r2_pred_byCond(bind,:) = 1 - (sum(currSSE_preds')./sum(SSE_predMatMean'))';
        SSE_pred_byCond(bind,:) = sum(currSSE_preds');
        ratioSSE_preds = (log2(currSSE_preds./SSE_predMatMean)); % conds X targs                
        currSSE_fit = zeros(totConds,totPredTargs);
        currSSE_fit(:,:) = SSE_fitMat(bind,:,:);
        r2_fit(bind) = 1 - (sum(currSSE_fit(:))/sum(SSE_fitMatMean(:)));
        r2_fit_byTargs(bind,:) = 1 - (sum(currSSE_fit)./sum(SSE_fitMatMean))';
        r2_fit_byCond(bind,:) = 1 - (sum(currSSE_fit')./sum(SSE_fitMatMean'))';
        SSE_fit_byCond(bind,:) = sum(currSSE_fit');  
        ratioSSE_fit = (log2(currSSE_fit./SSE_fitMatMean)); % conds X targs
        % figure out how many parameters per gene model
        currParamNums = zeros(kfoldCv,totPredTargs);
        currParamNums(:,:) = numParamsMat(bind,:,:);
        paramNumsStd(bind,:) = std(currParamNums);
        paramNumsMean(bind,:) = mean(currParamNums);
        paramNumsMin(bind,:) = min(currParamNums,[],1);
        paramNumsMax(bind,:) = max(currParamNums,[],1);
        % plot heatmaps of SSEs (genes X conditions) for each bootstrap cutoff
        % merge SSE_fit and SSE_pred, cluster by genes and conditions
        % 1. cluster conditions
%         condClustIn = [ratioSSE_preds ratioSSE_fit]; % cluster both pred + fit to get conditions order
%         pdis = pdist(condClustIn,'euclidean');
%         link = linkage(pdis,'average');
        
    end
    % plot R^2's    
    fontSize = 12;
    lineWidth = 2;
    xSize = 8;
    ySize = 11;
    figure
    subplot(2,1,1)
    plot(100*instabRangeNet,r2_pred,'r-+','LineWidth',2)
    hold on
    plot(100*instabRangeNet,r2_fit,'b-+','LineWidth',2)
    plot(100*instabRangeNet,r2_predNz,'g-o','LineWidth',2)
    plot(100*instabRangeNet,r2_fitNz,'k-o','LineWidth',2)
    xlabel('Instability cutoff','FontSize',fontSize)
    ylabel('R^2 (1 - SSE_{model}/SSE_{mean})','FontSize',fontSize)
    set(gca,'FontSize',fontSize)
    grid on, grid minor
    legend({'Predicted','Fit','nzPredicted','nzFit'},'Location','EastOutside')
     title(instabTitleInf,'FontSize',fontSize) 
    axis tight
    ax = axis();
    axis([ax(1:2) 0 1])
    
    % Plot Average Model Sizes
    subplot(2,1,2)
    
    totEdges = sum(numParams,2)/totPredTargs;
    plot(totEdges,r2_pred,'k-','LineWidth',2)
    hold on, grid on, set(gca,'FontSize',fontSize,'FontWeight','Bold')
    plot(median(numParams,2),median(r2_pred_byTarg,2),'b-','LineWidth',lineWidth)
    hold on
    plot(mean(numParams,2),mean(r2_pred_byTarg,2),'r-','LineWidth',lineWidth)
    legend({'Overall R^2_{pred}','Gene Median R^2_{pred}','Gene Average R^2_{pred}'},'Location','EastOutside')
    axis([0 60 0 1])

%     set(gca,'XTick',1:totModSizes,'XTickLabel',num2str((round(100*instabRangeNet'))),...
%         'FontSize',fontSize-2,'XTickLabelRotation',90)
    title(instabTitleInf,'FontSize',fontSize)
%     xlabel('Instability cutoff','FontSize',fontSize)
    grid on, grid minor
%     axis tight
%     ax = axis();
%     axis([ax(1) ax(2) 0 1])
%     set(gca,'yscale','log')
    xlabel('# of Predictors per Gene','FontSize',fontSize)
    ylabel('R^2 (1 - SSE_{model}/SSE_{mean})','FontSize',fontSize)
    currFig = [r2OutMat,'_v_preds' ];
    saveas(gcf,currFig,'fig')
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
    print('-painters','-dpdf','-r150',[currFig '.pdf'])
    disp(currFig)

    %% save variables of interest
    disp('Save variables of interest!')
    figOutBase = [r2OutMat];
    save([figOutBase '.mat'],...
        'SSE_predCvBreak',... = zeros(totModSizes,kfoldCv);
        'SSE_predMeanCvBreak',... = zeros(kfoldCv,1);
        'SSE_pred_kfoldInstabTargs',... = zeros(kfoldCv,totModSizes,totTargs);
        'SSE_predMean_kfoldInstabTargs',... = zeros(kfoldCv,totModSizes,totTargs);
        'r2_pred_kfoldInstabTargs',... = zeros(kfoldCv,totModSizes,totTargs);
        'instabRangeNet',...
        'modSizes',...
        'kcount',...
        'kfoldCv',...
        ...'priorName',... currInf{3};
        'SSE_predMatNzCV',...
        'SSE_predMatMeanNzCV',...
        'SSE_fitMatNzCV',...
        'SSE_fitMatMeanNzCV',...   
        'r2_fitNzCV',... = zeros(kfoldCv,totModSizes);
        'r2_predNzCV',... = zeros(kfoldCv,totModSizes);
        'SSE_predMat',... zeros(totModSizes,totConds,totPredTargs);
        'SSE_fitMat',... zeros(totModSizes,totConds,totPredTargs);
        'SSE_fitMatMean',... zeros(totConds,totPredTargs);
        'SSE_predMatMean',... zeros(totConds,totPredTargs);
        'r2_pred_byTarg',... zeros(totModSizes,totPredTargs);
        'r2_fit_byTargs',... zeros(totModSizes,totPredTargs);
        'r2_pred_byCond',... zeros(totModSizes,totConds);
        'r2_fit_byCond',... zeros(totModSizes,totConds);
        'SSE_pred',... zeros(totModSizes,totConds);
        'SSE_fit_mean',... zeros(totModSizes,totConds);
        'r2_pred',... zeros(totModSizes,1);
        'r2_fit',... zeros(totModSizes,1);
        'r2_predNz',... = zeros(totModSizes,1);
        'r2_fitNz',... = zeros(totModSizes,1);
        'SSE_pred_byCond',... zeros(totModSizes,totConds);    
        'SSE_fit_byCond',... zeros(totModSizes,totConds);
        'SSE_predMatNz',... = zeros(totModSizes,1);
        'SSE_predMatMeanNz',... = zeros(totModSizes,1);
        'SSE_fitMatNz',... = zeros(totModSizes,1);
        'SSE_fitMatMeanNz',... = zeros(totModSizes,1);
        'numParams',... zeros(totModSizes,kfoldCv*totPredTargs);
        'totEdges',...
        'numParamsMat',... zeros(totModSizes,kfoldCv,totPredTargs);
        'SSE_fitMat',... SSE_fitMat/(kfoldCv-1);    % now we have average SSE for fit
        'SSE_fitMatMean',... SSE_fitMatMean/(kfoldCv-1);    
        'paramNumsMean',... zeros(totModSizes,totPredTargs);
        'paramNumsStd',... zeros(totModSizes, totPredTargs);
        'paramNumsMax',... zeros(totModSizes, totPredTargs);
        'paramNumsMin')%,... zeros(totModSizes, totPredTargs);    
    disp([figOutBase '.mat created.'])
   
end                

