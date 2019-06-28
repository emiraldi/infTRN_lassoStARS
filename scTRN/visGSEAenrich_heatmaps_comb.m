%% visGSEAenrich_heatmaps_comb
% visualize enrichment of TFs' target genes in gene sets
% clusters sets and TFs based on signficance of enrichments
% assumes that GSEA has already been performed using tfTarget_GSEA_loop.sh
%   or tfTargets_GSEA.sh
% the resulting heatmap colors enrichment of TF's positive target genes red
%   and suppressed targets blue
%% Author: Emily Miraldi; Date: June 28, 2019

clear all
close all

addpath(fullfile('..','customMatlabFxns'));

%% Parameters
topNtfSets = Inf;   % keep the top N most significant sets per TF (set to Inf to include all sets per TFs)
FDR_cutoff = .01;    % cutoff for inclusion of enriched sets in heatmap
filterTfs = 1;      % filter TFs according to FDR_cutoff?, 1 --> on; 0 --> off
padjSat = .001;     % heatmap will be saturated for p-values more significant than this

% plotting
clusterTfs = 1;         % cluster TFs?, 1 --> on; 0 --> off
liveFigFormat = 1;      % enables user to interactively set figure dimensions
    % 1 --> on (requires user interaction), 0 --> off
fontSize = 7;           % font size for axes labels (TF and set names)
clim = -log10(padjSat); % saturation scale for the heatmap (-log10 p-values)

%% Gene sets database info
geneSetDBs = {'PathwayCommons';'Kegg'; 'MAPP'}; % 'Signatures_MSigDB'};
    % gene set databases must exactly match setNames in
    % tfTarget_GSEA_loop.sh
geneSetNickNames = {'P','K','Ma'};%,'Ms'};
    % OPTIONAL, user-defined gene set database nicknames (e.g., can specify
    % "K" for "Kegg database", so that "Cytokine Signaling" appears as 
    % "Cytokine Signaling (K)" in the heatmaps' gene-set labels; order must
    % match geneSetDBs
nickNamesOn = 1;    % 1 --> use nicknames, 0 --> don't use nicknames
totGeneSets = size(geneSetDBs,1);

%% Networks
% info for each GSEA network analysis is represented as a 4-column cell:
% 1. GSEA output folder format (with GENESETplaceholder for geneset database), 
%       (points to results from tfTarget_GSEA_loop.sh or tfTargets_GSEA.sh
% 2. optional list of TFs to limit analysis to (use empty string if all TFs 
%       are to be included in the visualization)
% 3. string to add to output file name
% 4. output directory
networkInfs = {... 
    {... % GSEA directory format "GENESETplaceholder" for geneset database name
    ['outputs/GRN_10X_pCorrCut0p01_combine_max/GSEA/GRN_10X_pCorrCut0p01_combine_max_GENESETplaceholder_Praw0p1_dir_wCut0p0_minSet5']; 
    '';         % .txt file list to limit TFs in heatmap, can be left as an empty string
    '_allTFs';  
    'outputs/GRN_10X_pCorrCut0p01_combine_max/GSEA/heatmaps'};
    };
totNetworks = size(networkInfs,1);

%% END parameters

setInf = [ 'combSets_fdr' num2str(100*FDR_cutoff) '_top' num2str(topNtfSets) 'tfSets']; % for output
titleInf = strrep([setInf],'_',' ');
disp(titleInf)       

for nind = 1:totNetworks
    tfsOfIntFile = networkInfs{nind}{2};
    outBit = networkInfs{nind}{3};
    outDir = networkInfs{nind}{4};
    mkdir(outDir)
    disp(outDir)

    %% set up p-val structure
    pvalStruct = [];
    pvalStruct(1).type = 'up';
    pvalStruct(2).type = 'down';
    pvalStruct(1).factor = 1;
    pvalStruct(2).factor = -1;
    totTypes = length(pvalStruct);

    %% load adjusted pvals per gene set database
    for gs = 1:totGeneSets    
        geneSetDB=geneSetDBs{gs};        
        outFileBase = strrep(networkInfs{nind}{1},'GENESETplaceholder',geneSetDB);
        setInfIn = ['' geneSetDB '_praw10'];   % for input set in table form
    
        %% get up and down-regulatory interaction p-vals
        for tind = 1:totTypes
            pvalFile = fullfile(outFileBase,[setInfIn '_' pvalStruct(tind).type '_adjp.txt']);
            fid = fopen(pvalFile);  % get first line, TFs
            tline=fgetl(fid);
            tfTmps = cellstr(strvcat(strsplit(tline,'\t')));
            fclose(fid);          
            totTfs = length(tfTmps);
            fid = fopen(pvalFile);  % get pvals + sets
            C = textscan(fid,['%s' repmat('%f',1,totTfs)],'Delimiter','\t','Headerlines',1);
            setsTmp = strrep(C{1},'_',' ');
            padjsTmp = [C{2:end}];
            totSets = length(setsTmp);
                        
            if nickNamesOn % optionally add nicknames to sets
                geneSetNickName = geneSetNickNames{gs};
                setNamesN = setsTmp;
                setsIn = length(setsTmp);
                for si = 1:setsIn
                    setsTmp{si} = [setNamesN{si} '(' geneSetNickName ')'];
                end
            end

            if gs == 1    % initialize for first gene set
                transVals = pvalStruct(tind).factor * -log10(max(padjsTmp,1E-10)); % multiply by factor
                pvalStruct(tind).pvals = transVals;
                pvalStruct(tind).sets = strvcat(setsTmp);
                pvalStruct(tind).tfs = tfTmps;
            else                          
                pvalStruct(tind).sets = strvcat(pvalStruct(tind).sets,strvcat(setsTmp)); % append sets
                transVals = pvalStruct(tind).factor * -log10(max(padjsTmp,1E-10));
                [newTfs,newTfInds] = setdiff(tfTmps,pvalStruct(tind).tfs);
                totNew = length(newTfInds);
                [oldTfs,oldTfInds] = setdiff(tfTmps,newTfs);
                oldTfNewLoc = find(ismember(pvalStruct(tind).tfs,oldTfs));
                pvalsForNewSetsOldTfs = zeros(totSets,size(pvalStruct(tind).pvals,2)); 
                pvalsForNewSetsOldTfs(:,oldTfNewLoc) = transVals(:,oldTfInds);
                pvalStruct(tind).tfs = cellstr(strvcat(strvcat(pvalStruct(tind).tfs),strvcat(newTfs))); % append new TFs
                pvalStruct(tind).pvals = [pvalStruct(tind).pvals zeros(size(pvalStruct(tind).pvals,1),totNew); % append p-values
                        pvalsForNewSetsOldTfs transVals(:,newTfInds)];
            end
        end
    end
    
    % convert sets to strings    
    % find indices corresponding to the end of + and - enrichments
    typeEndLoc = zeros(totTypes,1); 
    totSets = 0;
    for tind = 1:totTypes            
        pvalStruct(tind).sets = cellstr(pvalStruct(tind).sets);
        totSets = totSets + length(pvalStruct(tind).sets);        
        typeEndLoc(tind) = totSets;
    end
    
    %% optionally limit to TFs
    if length(tfsOfIntFile)
        fid = fopen(tfsOfIntFile,'r');
        C = textscan(fid,'%s','HeaderLines',0);
        fclose(fid);
        allTfs = C{1};
    else
        allTfs = union(pvalStruct(:).tfs);
    end
    totTfs = length(allTfs);

    %% order TFs between + and - enrichments
    allEnrichesTmp = zeros(totSets,totTfs);
    [xx,subInds,finalInds] = intersect(pvalStruct(1).tfs,allTfs);
    allEnrichesTmp(1:typeEndLoc(1),finalInds) = pvalStruct(1).pvals(:,subInds);
    [xx,subInds,finalInds] = intersect(pvalStruct(2).tfs,allTfs);
    allEnrichesTmp(typeEndLoc(1)+1:end,finalInds) = pvalStruct(2).pvals(:,subInds);
    setNamesTmp = cellstr(strvcat(strvcat(pvalStruct(1).sets),strvcat(...
        pvalStruct(2).sets)));

    %% threshhold sets based on FDR-cutoff
    setMaxes = max(abs(allEnrichesTmp),[],2);
    keepSets = find(setMaxes>-log10(FDR_cutoff));
    allEnriches = allEnrichesTmp(keepSets,:);
    setNames = {setNamesTmp{keepSets}}';
    
    %% if no TF file was provided, optionally filter TFs
    if filterTfs
        setMaxes = max(abs(allEnriches),[],1);
        keepTfs = find(setMaxes>-log10(FDR_cutoff));
        allEnriches = allEnriches(:,keepTfs);
        allTfs = {allTfs{keepTfs}}';
        totTfs = length(allTfs);
    end

    if and(topNtfSets,isfinite(topNtfSets)) % go tf by tf and find top N enriched sets
        topTfSets = [];
        for tind = 1:totTfs
            currTFvec = abs(allEnriches(:,tind));
            nonZero = find(currTFvec);
            [vals, inds] = sort(currTFvec,'descend');
            topTfSets = union(topTfSets,intersect(inds(1:topNtfSets),nonZero));
            disp(allTfs{tind})
            disp(strvcat(setNames{intersect(inds(1:topNtfSets),nonZero)}))
        end
        allEnriches = allEnriches(topTfSets,:);
        setNames = {setNames{topTfSets}}';
    end
   
    %% cluster all gene sets together
    binEnrich = sign(abs(allEnriches));
    
    % optionally cluster TFs
    if clusterTfs
        pdis = pdist(binEnrich');
        link = linkage(pdis,'ward');
        [h,t,horderTfs] = dendrogram(link,0);
    else
        horderTfs = 1:totTfs
    end

    % cluster sets, separately cluster + and - edge enrichments
    % 1. negative
    setSigns = sum(allEnriches,2);
    negInds = find(setSigns<0);
    currPs = binEnrich(negInds,:);
    if size(currPs,1) > 1
        pdis = pdist(currPs,'correlation');
        link = linkage(pdis,'ward');
        [h t horderNegs] = dendrogram(link,0);
        horderNegs = reshape(negInds(horderNegs),length(negInds),1);
    elseif size(currPs,1) == 1
        horderNegs = negInds;
        setNames{negInds};
    else
        horderNegs = [];
    end
    
    % 2. positive
    posInds = find(setSigns>0);
    currPs = binEnrich(posInds,:);
    if size(currPs,1) > 1
        pdis = pdist(currPs,'correlation');
        link = linkage(pdis,'ward');
        [h t horderPos] = dendrogram(link,0);
        horderPos = posInds(horderPos);
        horderPos = reshape(horderPos,length(posInds),1);
    elseif size(currPs,1) == 1
        horderPos = posInds;
    else 
        horderPos = [];
    end
    
    % 3. combine
    horderSets = [horderPos;horderNegs];
    typeEndLoc = length(horderPos);
    allEnriches = allEnriches(horderSets,horderTfs);
    setNames = {setNames{horderSets}}';

    %% generate a heatmap
    figure,
    subplot(5,5,[4:5,8:10, 13: 15, 19:20])
    imagesc(allEnriches)
    colormap redblue
    hold on
    ax = axis();
    grid on

    set(gca,'XTick',1:totTfs,'XTickLabel',strvcat(strrep(({allTfs{horderTfs}}),'_',' ')),...
        'FontSize',fontSize,'FontWeight','Bold')
    set(gca,'YTick',1:length(setNames),'YTickLabel',setNames,...
        'XTickLabelRotation',90)
    axis image
    set(gca,'CLim',[log10(padjSat) -log10(padjSat)],'TickLength',[0 0])
        colorbar('YTick',[log10(padjSat) -1 0 1 -log10(padjSat)],'YTickLabel',...
            {['Repressive, P_{adj} <= ' roundstring3(padjSat)];'Repressive, P_{adj} = .1';...
            'P_{adj} = 1';'Activating, P_{adj} = .1';['Activating, P_{adj} <= ' roundstring3(padjSat)]},...
            'FontWeight','Bold')

    title(titleInf,'FontSize',fontSize + 2)
    currFig = fullfile(outDir, [setInf outBit]);
    saveas(gcf,currFig,'fig')
    if liveFigFormat
        disp('Hit enter once satisfied with figure dimensions.')
        pause
    end
    save2pdf(currFig,gcf,150)
    disp(currFig)

end
