function calcPRinfTRNs(infTrnFile,gsFile,rankColTrn,...
    targGeneFile,gsRegsFile,outFileBase,figOutBase)
%% calcPRinfTRNs(infTrnFile,gsFile,rankColTrn,...
%    targGeneFile,gsRegsFile,outFileBase,saveFigs)
%% GOAL: calculates network precision-recall (and TPR-FPR) for a gold standard
% P-R / ROC is provided for the full G.S. network, as well as each TF in
% the G.S. individually
%% Reference:
% Miraldi et al. "Leveraging chromatin accessibility data for 
%   transcriptional regulatory network inference in T Helper 17 Cells"
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% INPUTS:
% infTrnFile -- sparse network file (col 1 = TF, col 2 = target gene, col
%   "rankColTrn" = edge confidence.  NOTE: Code ranks edges based on
%   ABSOLUTE VALUE of edge confidence
% gsFile -- gold standard in sparse network file form (col 1 = TF, col 2 = 
%   target gene, col 3 = edge confidence). NOTE: Any nonzero TF-gene
%   interaction will be included in the gold standard.
% rankColTrn -- column number corresponding to the values in the infTrnFile
%   that will be used to rank interactions
% targGeneFile -- the set of target genes to be considered for P-R analysis
%   (e.g., could be the intersection of genes (1) for which TRNs were built 
%   and (2) that were evaluated in constructing the gold standard)
%   OPTION: Providing '' (empty string argument) will result in inclusion of
%   all target genes in the gsFile for P-R and ROC
% gsRegsFile -- list of regulators to be considered for P-R analysis
%   OPTION: Providing '' (empty string argument) will result in inclusion of
%   all regulators in the gsFile for P-R and ROC
% outFileBase -- output file name used for ${outFileBase}.mat as well as
%   figure output ${outFileBase}.pdf ${outFileBase}.fig
% figOutBase -- name for figure files, 
%   OPTION: Provide '' (empty string argument) to avoid saving figure
%% OUTPUTS:
% outFileBase -- contains precision, recall, TPR, FPR for the network
%   overall as well as each gene individually
% (optional) figOutBase.fig / .pdf of P-R and ROC for full network

%% load genes considered
if targGeneFile
    geneIn = fopen(targGeneFile,'r');
    C = textscan(geneIn,'%s');
    fclose(geneIn);
    potTargGenes = C{1};
    totTargGenes = length(potTargGenes);
end
%% load regulators considered
if gsRegsFile
    geneIn = fopen(gsRegsFile,'r');
    C = textscan(geneIn,'%s');
    fclose(geneIn);
    potRegs = C{1};
    totGsRegs = length(potRegs);
end
%% Load the gold standard
gsInfs = [];    % file base, text, regs, reg-targs, random P-R
gsInfs.regs = {};      % keep track of the gs TFs across all sets, and target genes, because LASSO networks are big
fid = fopen(gsFile,'r');
% get first line and see what regulators we have
tline = fgetl(fid);    
totCols = length(cellstr(strsplit(tline,'\t')));
fclose(fid);
% get the rest of the data using textscan
fid = fopen(gsFile,'r');
if totCols == 3
    C = textscan(fid,['%s%s%f'],'Delimiter','\t','Headerlines',1);
elseif totCols > 3
    C = textscan(fid,['%s%s%f' repmat('%s',1,totCols-3)],'Delimiter','\t','Headerlines',1);
else
    error('Prior matrix should be in sparse format (TF,gene,edge weight).')
end
fclose(fid);
gsRegsTmp = C{1};
gsTargsTmp = C{2};
gsWeights = C{3};
keepWeights = find(gsWeights); % indices of edges with nonzero weight
% limit to TF-gene interactions considered by the model
if gsRegsFile
    gsRegInds = find(ismember(gsRegsTmp,potRegs));
else
    gsRegInds = 1:length(gsRegsTmp); % consider all TFs
end
if targGeneFile
    gsTargInds = find(ismember(gsTargsTmp,potTargGenes));
else
    gsTargInds = 1:length(gsRegsTmp); % consider all target genes
end
keepInds = intersect(intersect(gsRegInds,gsTargInds),keepWeights);
gsRegs = cellstr(strvcat(gsRegsTmp{keepInds}));
gsTargs = cellstr(strvcat(gsTargsTmp{keepInds}));
totGsInts = length(keepInds);
uniGsRegs = unique(gsRegs);
totGsRegs = length(uniGsRegs);
gsInfs.regs = uniGsRegs;
if not(length(targGeneFile))  % all targets in the prior are used if no target gene list was supplied
    potTargGenes = unique(gsTargs);
    totTargGenes = length(potTargGenes);
end
gsInfs.totPotInts = totTargGenes*totGsRegs;
randPR = totGsInts/gsInfs.totPotInts;
[dd, currFileBase,ext] = fileparts(gsFile);
gsInfs.fileBase = gsFile;
gsInfs.figText = strrep(strrep(currFileBase,'_sp',''),'_',' '); % use file name 
        
ints = [gsRegs gsTargs];
totInts = length(ints);
gsInfs.edges = cell(totInts,1);
for ii = 1:totInts
    gsInfs.edges{ii} = strjoin({ints{ii,:}},',');
end    
gsInfs.randPR = randPR;
gsInfs.targs = unique(gsTargs);
    
gsInfs.edgesByTf = cell(totGsRegs,1);
gsInfs.randAuprByTf = zeros(totGsRegs,1);
for gind = 1:totGsRegs
    currInds = find(ismember(gsRegs,uniGsRegs(gind)));
    gsInfs.edgesByTf{gind} = {gsInfs.edges{currInds}}';
    gsInfs.randAuprByTf(gind) = length(gsInfs.edgesByTf{gind})/totTargGenes;
end        
% we will fill in these values and then save a version of gsInfs for
% each network
gsInfs.precisions = {};
gsInfs.recalls = {};
gsInfs.fprs = {};
gsInfs.arocs = NaN;
gsInfs.auprs = NaN;
gsInfs.auprsByTf = zeros(totGsRegs,1);
gsInfs.arocsByTf = zeros(totGsRegs,1);
gsInfs.precisionsByTf = cell(totGsRegs,1);
gsInfs.recallsByTf = cell(totGsRegs,1);
gsInfs.fprsByTf = cell(totGsRegs,1);
    
fprintf([num2str(totGsInts) '\t' gsInfs.figText '\n']);

%% Calculate ROC, P-R performance
%% Need to get: regulators, edges, edge weights
% check number of columns
disp(infTrnFile)
fid = fopen(infTrnFile,'r');            
tline = fgetl(fid);
tline = fgetl(fid); % use 2nd line, to make sure that it works on R output (often missing first tab in a table)
fclose(fid);
totInfOuts = length(strsplit(tline,'\t')');
% import the data
fid = fopen(infTrnFile,'r');
C = textscan(fid,[repmat('%s',1,totInfOuts)],'Delimiter','\t','Headerlines',1);
fclose(fid);    
regsRaw = C{1};
targs = C{2};

rankings = abs(str2double([C{rankColTrn}]));  % assume that distance from zero is proportional to confidence
% limit network edges to specified TFs and target genes
% limit to TF-gene interactions considered by the model
if gsRegsFile
    gsRegInds = find(ismember(regsRaw,potRegs));
else
    gsRegInds = find(ismember(regsRaw,gsRegs)); % consider all TFs
end
if targGeneFile
    gsTargInds = find(ismember(targs,potTargGenes));
else
    gsTargInds = find(ismember(targs,gsTargs)); % consider all target genes
end
keepEdges = intersect(gsRegInds,gsTargInds);
regs = {regsRaw{keepEdges}}';
targs = {targs{keepEdges}}';
rankings = rankings(keepEdges);
stabRange = unique(rankings);
ints = [regs targs];
totTrnInts = length(ints);
infEdges = cell(totTrnInts,1);
for ii = 1:totTrnInts
    infEdges{ii} = strjoin({ints{ii,:}},',');
end   

gsRegs = gsInfs.regs;
gsEdges = gsInfs.edges;
totGsRegs = length(gsRegs);

% limit analysis to regs specific to current gold standard
setdiffEdges = setdiff(gsEdges,infEdges);
totGsEdges = length(setdiffEdges);
unpredictedEdges = gsInfs.totPotInts - totTrnInts - totGsEdges;

inputValues = [rankings; zeros(totGsEdges+unpredictedEdges,1)];                        
setList = [ismember(infEdges,gsEdges);
    ones(totGsEdges,1);zeros(unpredictedEdges,1)];
disp(['Total Interactions w/ GS TFs (' num2str(totGsRegs) '): ' num2str(totTrnInts)])                
        
if totTrnInts > 0 % there was at least one G.S. TF in the network
    [aupr,...
        aroc,...
        precisionV, recallV, fprV, stepVals, f1scores]...
            = aupr_step_outVals(inputValues,...
                setList);   
    gsInfs.precisions = precisionV;
    gsInfs.recalls = recallV;
    gsInfs.fprs = fprV;
    gsInfs.auprs = aupr;
    gsInfs.arocs = aroc;                       
    gsInfs.f1scores = f1scores;
    gsInfs.stepVals = stepVals;                
else
    disp(['No TFs from the G.S. ' gsInfs.figText ' found.'])             
end

[dd, currFileBase,ext] = fileparts(infTrnFile);
infTitleBase = strrep(strrep(currFileBase,'_sp',''),'_',' '); % use file name 

% P-R
subplot(2,1,1)
plot([0 1],gsInfs.randPR*[1 1],':','LineWidth',2,'Color',.75*[1 1 1]) % random
hold on
plot(gsInfs.recalls,gsInfs.precisions,'b-','LineWidth',2) % method performance
grid on, grid minor, box on
axis([0 1 0 1])
xlabel('Recall','FontSize',12)
ylabel('Precision','Fontsize',12)
title([infTitleBase],'FontSize',14)    

% ROC
subplot(2,1,2)
plot([0 1],[0 1],':','LineWidth',2,'Color',.75*[1 1 1]) % random
hold on
plot(gsInfs.fprs,gsInfs.recalls,'b-','LineWidth',2) % method performance
grid on, grid minor, box on
axis([0 1 0 1])
xlabel('FPR','FontSize',12)
ylabel('TPR','Fontsize',12)
title([infTitleBase],'FontSize',14)  

if figOutBase    
    saveas(gcf,figOutBase,'fig')
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [4.25 8]);
    print('-painters','-dpdf','-r150',[figOutBase '.pdf'])
    disp([figOutBase '.fig + .pdf created.'])
end

%% Breakdown By TF
gsRegs = gsInfs.regs;
totGsRegs = length(gsRegs);

for gind = 1:totGsRegs
    currReg = gsRegs{gind};
    currGsEdges = gsInfs.edgesByTf{gind};                        
    infInds = find(ismember(regs,currReg));
    infInputEdges = {infEdges{infInds}}';
    currRankings = rankings(infInds);
    setdiffEdges = setdiff(currGsEdges,infInputEdges);
    totGsOnlyEdges = length(setdiffEdges);
    unpredictedEdges = totTargGenes - length(infInputEdges) - totGsOnlyEdges;
    inputValues = [currRankings; zeros(totGsOnlyEdges+unpredictedEdges,1)];      
    predMade = ismember(infInputEdges,currGsEdges);
    if length(predMade) > 1                
        unpred = [ones(totGsOnlyEdges,1);zeros(unpredictedEdges,1)];
        perm4rest = randperm(length(unpred));
        setList = [predMade;
            unpred(perm4rest)];
        % get predictions and gold standard ready for aupr curve
        [aupr, aroc, precisionV, recallV, fprV, uniInVals, f1scores]...
            = aupr_step_outVals(inputValues,setList);
    else % there were no predictions, so random performance
        aupr = gsInfs.randAuprByTf(gind);
        aroc = .5;
        precisionV = [aupr aupr];                
        recallV = [0 1];
        fprV = [0 1];
%                 disp(['Total overlap = ' num2str(length(ismember(infEdges,currGsEdges)))])
    end                
    gsInfs.precisionsByTf{gind} = precisionV;
    gsInfs.recallsByTf{gind} = recallV;
    gsInfs.fprsByTf{gind} = fprV;
    gsInfs.auprsByTf(gind) = aupr;
    gsInfs.arocsByTf(gind) = aroc;   
    disp([gsInfs.figText ' ' currReg ', AUPR = ' num2str(aupr) ', AROC = ' num2str(aroc)])
end

%% save results
save([outFileBase '.mat'],...
    'stabRange',...
    'gsInfs',...
    'regs',...
    'rankings',...
    'infEdges',...
    'infTitleBase')                
disp(outFileBase)

disp('Finished')