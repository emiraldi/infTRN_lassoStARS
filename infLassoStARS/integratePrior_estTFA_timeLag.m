function integratePrior_estTFA_timeLag(geneExprMat,priorFile,edgeSS,...
     minTargets, tfaOutMat, timeLagFile, timeLag)
%% integratePrior_estTFA_timeLag(geneExprMat,priorFile,edgeSS,...
%     minTargets, tfaOutMat, timeLagFile, timeLag)
%% GOALS
% Integrate prior information with gene expression data for
%   1. TRN inference based on TF mRNA
%   2. TRN inference based on prior-based TFA with the following equation:
%                         X = P * A,    [Equation 1]
% where P  [genes X TFs] is the prior matrix of known TF-gene interactions, 
% X  [genes X samples] is the expression matrix for genes in the prior
% A  [TFs X samples] contains the unknown protein activities for TFs in the
% prior. (Ortiz et al. (2015) Molec. Sys. Bio.)
% Here, we introduce a time-lag, as proposed by Bonneau et al. 2006. Genome
%   Biology. For Sample Q taken at time point at U time units proceeded
%   temporally by Sample P taken at time point U-V time units:
%       TFA(Q) = TFA(Q) - ((TFA(Q) - TFA(P)) / V ) * (timelag)
%   where timelage represents (1) delay between TF mRNA becoming active
%   protein TF and (2), for prior-based TFA, the biological interpretation
%   of the time lag is unclear (but has been done in at least one other
%   study (Tchourine, Vogel and Bonneua (2018) Cell))
%% INPUTS:
% geneExprMat -- .mat file of gene expression and gene lists for TFA 
%   estimation (e.g., as generated by importGeneExprGeneLists.m)
% priorFile -- a prior in table form, rows = genes, columns = regulators
% edgeSS -- number of prior edge subsamples, NOTE: if edgeSS = 0, no 
%   subsampling will occur and all edges will be used to calculate TFA
% minTargets -- minimum number of targets per TF to stay in prior
% tfaOutMat -- filename for TFA output .mat
% timeLagFile -- four-column text file (w/ header), each row corresponds to 
%   a time point that needs to have time-lagged TFA. col 1 = Sample Name Q, 
%   col 2 = time point Sample Q (number), col 3 = Sample Name P, which
%   immediately proceeds Sample Q temporally, col 4 = time point Sample P.
%   NOTE: The timeLagFile can be a subset of total samples, as not every
%       sample need be part of a time-series with a proceeding time point
% timeLag -- time lag, or expected delay between TF gene expression and 
%       protein TFA, should be in same units as timeLagFile (hr, min...)
%% OUTPUTS:
% tfaOutMat -- whos
%   medTfas                 341x254              692912  double              
%   noPriorRegs             369x1                 45554  cell                
%   noPriorRegsMat          369x254              749808  double              
%   pRegs                   341x1                 41542  cell                
%   pRegsNoTfa              346x1                 41996  cell                
%   pTargs                14485x1               1792680  cell                
%   pTargsNoTfa           14513x1               1796176  cell                
%   priorMatrix           14485x341            39515080  double              
%   priorMatrixNoTfa      14513x346            40171984  double 
%% NOTE (!!): If there are TFs with identical targets in the prior, it is
% assumed that there is also a second priorFile (specifically, 
% mergedFile = strrep(priorFile,'.tsv','_merged.tsv') and associated file
% mergedTFs = strrep(priorFile,'.tsv','_mergedTFs.txt') to decode gene
% names of merged TFs, as can be generated by mergeDegeneratePriorTFs.py 
% in the priorParsingFxns folder associated with this repo)
%% Reference:
% Miraldi et al. "Leveraging chromatin accessibility for 
%   transcriptional regulatory network inference in T Helper 17 Cells"
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital

%% debugging
% load('~/erm/infTRN_lassoStARS/infLassoStARS/timeLag.mat')


%% load input gene expression
load(geneExprMat)
[totGenes, totConds] = size(targGeneMat);
  
%% Case 1: TF mRNA, use the unmerged prior, limit prior to potential 
% regulator list and target genes
fid = fopen(priorFile,'r');
% get first line and see what regulators we have    
tline = fgetl(fid);    
pRegsTmp = strsplit(tline,'\t');
pRegsTmp = cellstr(strvcat(pRegsTmp)); % get rid of first \t if it exists
totPRegs = length(pRegsTmp);        
fclose(fid);
% get the rest of the data using textscan
fid = fopen(priorFile,'r');
C = textscan(fid,['%s' repmat('%f',1,totPRegs)],'Delimiter','\t','Headerlines',1);
fclose(fid);
pTargsTmp = C{1};
pIntsTmp = [C{2:end}];    
        
%% Limit prior regulators and gene targets to those included in input data (.mat)
[pRegsNoTfa, pRegIndsTmp, potRegInds] = intersect(pRegsTmp,potRegs);
[pTargsNoTfa, pTargIndsTmp, expTargInds] = intersect(pTargsTmp,tfaGenes);
priorMatrixNoTfa = pIntsTmp(pTargIndsTmp,pRegIndsTmp); 

% output mRNA estimates for regulators without TFA
[noPriorRegs, expInds] = setdiff(potRegs_mRNA,pRegsNoTfa);
noPriorRegsMat = potRegMat_mRNA(expInds,:);

%% Case 2: prior-based TFA is used
% Check whether there were degenerate TFs and outputs from e.g.,
% mergeDegenreatePriorTFs.py exist
mergedTFsExist = 1;
try
    mergedFile = strrep(priorFile,'.tsv','_merged.tsv');
    mergedTFs = strrep(priorFile,'.tsv','_mergedTFs.txt');
    ls(mergedFile)
    ls(mergedTFs)
catch
    mergedTFsExist = 0;
end

if mergedTFsExist
    mtIn = fopen(mergedTFs,'r');
    C = textscan(mtIn,'%s%s','Delimiter','\t','HeaderLines',0);
    mergedTFs = C{1};
    individualTFs = C{2};
    totMergedSets = length(individualTFs);
    keepMergedTFs = []; % keep track of merged TF names to keep
    for mind = 1:totMergedSets
        currSet = strsplit(individualTFs{mind},', ');
        usedTfs = intersect(currSet,potRegs);
        if length(usedTfs)
            keepMergedTFs = [keepMergedTFs; mind];
        end
    end
    if length(keepMergedTFs) % add merged potential regulators to our list       
        potRegs = union(potRegs,{mergedTFs{keepMergedTFs}});
        priorFile = mergedFile;
        % regulator list and target genes
        fid = fopen(priorFile,'r');
        % get first line and see what regulators we have    
        tline = fgetl(fid);    
        pRegsTmp = strsplit(tline,'\t');
        pRegsTmp = cellstr(strvcat(pRegsTmp)); % get rid of first \t if it exists
        totPRegs = length(pRegsTmp);        
        fclose(fid);
        % get the rest of the data using textscan
        fid = fopen(priorFile,'r');
        C = textscan(fid,['%s' repmat('%f',1,totPRegs)],'Delimiter','\t','Headerlines',1);
        fclose(fid);
        pTargsTmp = C{1};
        pIntsTmp = [C{2:end}];   
    end        
end
%% now go through and enforce minimum number of targets
%% Limit prior regulators and gene targets to those included in input data (.mat)
[pRegs, pRegIndsTmp, potRegInds] = intersect(pRegsTmp,potRegs);
[pTargs, pTargIndsTmp, expTargInds] = intersect(pTargsTmp,tfaGenes);
pInts = pIntsTmp(pTargIndsTmp,pRegIndsTmp);    
intsPerTf = sum(abs(sign(pInts)));
%% remove predictors that have fewer than the specified minimum number of targets
minTargets = max(minTargets,0);
keepRegs = find(intsPerTf>minTargets);
pInts = pInts(:,keepRegs);
pRegs = cellstr(strvcat(pRegs{keepRegs}));  
% possible that there now are some target genes with no regulators
intsPerTarg = sum(abs(sign(pInts)),2);
keepTargs = find(intsPerTarg > 0);
pInts = pInts(keepTargs,:);
pTargs = cellstr(strvcat(pTargs{keepTargs}));      
% make sure that target indices in prior and expression matrix
% coordinate
[pTargs2, pTargIndsTmp, expTargInds] = intersect(pTargs,tfaGenes);
priorMatrix = pInts(pTargIndsTmp,:);
% get target gene matrix
targExp = tfaGeneMat(expTargInds,:);

%% calculate TFA, ( X = P * A, A = P \ X)
[totTargs, totPreds] = size(priorMatrix);
if edgeSS > 0
    tfas = zeros(edgeSS,totPreds,totConds);
    for ss = 1:edgeSS
        sPrior = zeros(totTargs,totPreds);
        for col = 1:totPreds
            currTargs = priorMatrix(:,col);
            targInds = find(currTargs);
            totCurrTargs = length(targInds);
            ssample = targInds(randperm(totCurrTargs,ceil(.63*totCurrTargs)));
            sPrior(ssample,col) = priorMatrix(ssample,col);
        end        
        tfas(ss,:,:) = sPrior \ targExp;    
    end
    medTfas = zeros(totPreds,totConds);
    medTfas(:,:) = median(tfas,1);
    disp(['Median from ' num2str(edgeSS) ' subsamples used for prior-based TFA.'])
else
    medTfas = priorMatrix \ targExp;
    disp(['No subsampling for prior-based TFA estimate.'])
end

%% load time-lag info
fid = fopen(timeLagFile,'r');
C = textscan(fid,'%s%f%s%f','Delimiter','\t','Headerlines',1);
fclose(fid);

timePointConds = C{1};
timePointProceedConds = C{3};
distanceBetweenPoints = C{2}-C{4}; % current - proceeding time point
totTimes = length(timePointConds);

noPriorRegsMatTmp = noPriorRegsMat;
medTfasTmp = medTfas;
potRegMat_mRNA_tmp = potRegMat_mRNA;

proceedInds = zeros(totTimes,1);
currInds = zeros(totTimes,1);

for tp = 1:totTimes
    currPoint = timePointConds{tp};
    currInd = find(ismember(conditionsc,currPoint));
    currInds(tp) = currInd;
    proceedPoint = timePointProceedConds{tp};
    proceedInd = find(ismember(conditionsc,proceedPoint));
    proceedInds(tp) = proceedInd;
    % TF mRNA in prior
    delta_TFmRNA = noPriorRegsMatTmp(:,currInd) - noPriorRegsMatTmp(:,proceedInd);
    slopes_TFmRNA = delta_TFmRNA/distanceBetweenPoints(tp);    
    noPriorRegsMat(:,currInd) = noPriorRegsMat(:,currInd) - slopes_TFmRNA*timeLag;
    % TFA in prior
    delta_TFA = medTfasTmp(:,currInd) - medTfasTmp(:,proceedInd);
    slopes_TFA = delta_TFA/distanceBetweenPoints(tp);
    medTfas(:,currInd) = medTfas(:,currInd) - slopes_TFA*timeLag;    
    % TF mRNA
    delta_TFmRNA = potRegMat_mRNA_tmp(:,currInd) - potRegMat_mRNA_tmp(:,proceedInd);
    slopes_TFmRNA = delta_TFmRNA/distanceBetweenPoints(tp);    
    potRegMat_mRNA(:,currInd) = potRegMat_mRNA(:,currInd) - slopes_TFmRNA*timeLag;
    
end

% vis TFA
% tfName = 'Rorc';
% rorInd = find(ismember(pRegs,tfName));
% figure, plot([C{4} C{2}]',[medTfasTmp(rorInd,proceedInds)' medTfasTmp(rorInd,currInds)']','o-')
% hold on
% plot([C{2}-.5],[medTfas(rorInd,currInds)'],'*')
% xlabel('Time','FontSize',14)
% ylabel(tfName,'FontSize',14)
% grid on
% grid minor
% set(gca,'FontSize',12)
% title('time-lagged TFA, .5 hours','FontSize',14)
% 
% % vis TF mRNA
% rorInd = find(ismember(potRegs_mRNA,tfName));
% figure, plot([C{4} C{2}]',[potRegMat_mRNA_tmp(rorInd,proceedInds)' potRegMat_mRNA_tmp(rorInd,currInds)']','o-')
% hold on
% plot([C{2}-.5],[potRegMat_mRNA(rorInd,currInds)'],'*')
% xlabel('Time','FontSize',14)
% ylabel(tfName,'FontSize',14)
% grid on
% grid minor
% set(gca,'FontSize',12)
% title('time-lagged TF mRNA, .5 hours','FontSize',14)

% disp('Note: We overwrite')
% disp('potRegMat_mRNA in the geneExprMat, so the tfaMat must be loaded')
% disp('beforehand in e.g., as in estimateInstabilitiesTRNbStARS.')

% output the .mat file
save(tfaOutMat,...
    'pTargsNoTfa',...
    'priorMatrixNoTfa',...
    'pRegsNoTfa',...
    'medTfas',...
    'pRegs',...
    'pTargs',...
    'priorMatrix',...
    'noPriorRegs',...
    'potRegMat_mRNA',...
    'noPriorRegsMat')
disp(tfaOutMat)
 