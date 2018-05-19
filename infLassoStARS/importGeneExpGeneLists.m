function importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,...
    tfaGeneFile,outMat)
%% function importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,...
%    tfaGeneFile,outMat)
%% GOALS:
% 1. get target gene matrix for all relevant conditions
% 2. get mRNA levels of predictors in a matrix
% 3. get the gene expression matrix that can be used for TFA based on prior
%       knowledge of TF-target relationships
% 4. Put these items in a .mat object to be used for mLASSO-StARS
%% INPUTS:
% normGeneExprFile -- tab-delimited table of gene expression data (normalized
%   across samples but not necessarily across genes), columns correspond to
%   samples, rows correspond to genes
% targGeneFile -- text file with list of target genes for gene models
% potRegFile -- text file with list of potential regulators to be
%   considered as predictors in the gene models
% tfaGeneFile -- text file list of genes that should be included in TFA 
%   estimation using TFA = pseudoinverse(P) * X,... TFA will not be
%   calculated in this script
%   NOTE: Set to '', to use all genes, which is the standard default.
% outMat -- output .mat file name
%% OUTPUTS:
% outMat -- .mat file containing, gene expression matrices, row and column
%   labels, names of potential regulators, mRNA expression if available for
%   potential regulators...
%% Reference:
% Miraldi et al. "Leveraging chromatin accessibility for 
%   transcriptional regulatory network inference in T Helper 17 Cells"
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital

eps = 1E-10; % target genes whose standard deviation across all samples is 
% less than eps will be removed from the target gene matrix

%% input gene expression data
currFile = normGeneExprFile;
ls(currFile)
fid = fopen(currFile);
tline=fgetl(fid); % line 1
tline2 = fgets(fid); % line 2
totSamps = length(cellstr(strvcat(strsplit(tline2,'\t')))') - 1; % number of samples
conditionsc = cellstr(strvcat(strsplit(tline,'\t')))';
conditionsc = conditionsc(end-totSamps+1:end);
fclose(fid);
% get input data
fid = fopen(currFile);
C = textscan(fid,['%s' repmat('%f',1,totSamps)],'Delimiter','\t','Headerlines',1);
fclose(fid);
genesc = C{1};
ncounts = [C{2:end}];
[vars obs] = size(ncounts);
 
%% load target genes, predictors, nominally expressed genes and get matrices for each
%% target genes
fid = fopen(targGeneFile,'r');
C = textscan(fid,'%s','Headerlines',0);
fclose(fid);
targGenesTmp = C{1};
% get target gene matrix
[xx, indsMat] = intersect(genesc,targGenesTmp);
targGenes = cellstr(strvcat(genesc{indsMat}));
targGeneMat = ncounts(indsMat,:);
% make sure that genes actual show variation
stds = mad(targGeneMat')';
Zstd = find(stds<eps);
remove = strvcat(targGenes{Zstd});
if remove
    disp('Target gene without variation, removed from analysis:')
    disp(remove)
end
keep = setdiff(1:length(indsMat),Zstd);
targGenes = {targGenes{keep}};
targGeneMat = targGeneMat(keep,:);
fprintf([num2str(length(targGenes)) ' target genes total.\n\n'])
% see if any target genes failed to map
if length(targGenesTmp) > length(targGenes)
    [missing, missInds] = setdiff((targGenesTmp),(targGenes));
    disp(['The following ' num2str(length(missing)) ' target genes were not found:'])
    fprintf([strjoin({targGenesTmp{missInds}},'\n') '\n\n']);
end

%% import potential regulators and find mRNA expression levels, if available
fid = fopen(potRegFile,'r');
C = textscan(fid,'%s','Headerlines',0);
fclose(fid);
potRegs = C{1};
% get mRNA levels of potential regulators matrix
[xx, indsMat] = intersect((genesc),(potRegs));
potRegs_mRNA = cellstr(strvcat(genesc{indsMat}));
potRegMat_mRNA = ncounts(indsMat,:);
fprintf([num2str(length(potRegs_mRNA)) ' potential regulators with expression data.\n\n'])
% see if any regulators lack gene expression data
if length(potRegs) > length(potRegs_mRNA)
    [missing, missInds] = setdiff((potRegs),(potRegs_mRNA));
    disp(['The following ' num2str(length(missing)) ' regulators have no expression data:'])
    fprintf([strjoin({potRegs{missInds}},'\n') '\n\n']);
end

%% genes to be included in TFA estimation (e.g., nominally expressed)
if tfaGeneFile
    fid = fopen(tfaGeneFile,'r');
    C = textscan(fid,'%s','Headerlines',0);
    fclose(fid);
    tfaGenesTmp = C{1};
else
    disp('No TFA gene file found, all genes will be used to estimate TFA.')
    tfaGenesTmp = genesc;
end
% get TFA gene matrix
[xx, indsMat] = intersect((genesc),(tfaGenesTmp));
tfaGenes = cellstr(strvcat(genesc{indsMat}));
tfaGeneMat = ncounts(indsMat,:);

%% SAVE target, TFA target mRNA levels, etc.
save(outMat,...
    'conditionsc',...
    'genesc',...
    'potRegMat_mRNA',...
    'potRegs',...
    'potRegs_mRNA',...
    'targGeneMat',...
    'targGenes',...
    'tfaGeneMat',...
    'tfaGenes')

disp(['Saved: ' outMat])