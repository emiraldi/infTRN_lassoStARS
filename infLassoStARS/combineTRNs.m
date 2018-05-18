function combineTRNs(combinedNetTsv,combineOpt,meanEdgesPerGene,...
    nets2combine)
%% combineTRNs(combinedNetTsv,combineOpt,meanEdgesPerGene,...
%   nets2combine)
%% GOAL: Algebraically combine networks, current combination options 
% include "max" (maximum) and "mean". Briefly, max will leverage 
% non-redundant strengths of individual models, while mean highlights edges
% supported by both models. Useful background is contained in:
%% Kittler et al. (1998) "On combining classifiers". IEEE.
%% NOTE: File is limited to combining networks from .mat TRN outputs from
%   buildTRNs_mLassoStARS.m, so that resulting sparse network file has
%   annotations for visualization in jp_gene_viz (e.g., edge width, style,
%   color
%% INPUTS:
% combinedNetTsv -- an output name for the resulting sparse network file
%   format, as described below
% combineOpt -- 'max' or 'mean' are the only two options
% meanEdgesPerGene -- an upper bound on size of resulting combined network,
%   corresponding to the average number of TFs per gene
% nets2combine -- a column cell (N X 1) of N .mat TRN outputs from 
%   buildTRNs_mLassoStARS.m to be combined
%% OUTPUTS:
% combinedNetTsv -- 3-column network file format for visualization in 
%       jp_gene_viz, limit models to size "meanEdgesPerGene"
%        0.  Edge confidence is a quantile, where total edges is set to
%              meanEdgesPerGene * total Gene Models
%        1.  Edge thickness (in output sparse network) is max or average
%               stability of combined TRNs
%        2.  Edges signs are calculated based on max or average partial 
%               correlation of combined TRNs
%% Additional Reference:
% Miraldi et al. "Leveraging chromatin accessibility data for 
%   transcriptional regulatory network inference in T Helper 17 Cells"
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: April 22, 2018

%% Network Edge colors, for jp_gene_viz
medBlue = [0,85,255];
medRed = [228,26,28];
lightGrey = [217,217,217];

%% Analysis Parameters
totNets = size(nets2combine,1);

%% 1. Consider all possible edges when combining networks
allEdges = {};
allRegs = {};
allTargs = {};
infEdgeCollections = cell(totNets,1);
for nind = 1:totNets
    load(nets2combine{nind})
    totInts = length(regs);
    infEdges = cell(totInts,1);
    for ii = 1:totInts
        infEdges{ii} = [regs{ii} ',' targs{ii}];
    end    
    allEdges = union(allEdges,infEdges);
    allRegs = union(allRegs,regs);
    allTargs = union(allTargs,targs);
    infEdgeCollections{nind} = infEdges;
end
totInfInts = length(allEdges);
totTargs = length(allTargs);
aveQuantVec0 = zeros(totInfInts,1);
aveStabVec0 = zeros(totInfInts,1);
aveCoefVec0 = zeros(totInfInts,1);
aveInPriorVec0 = zeros(totInfInts,1);

if find(ismember({combineOpt},'max'))
    disp('max-combine')
    for nind = 1:totNets
        load(nets2combine{nind})
        infEdges = infEdgeCollections{nind};
        [vals, currIndsAll, currIndsInd] = intersect(allEdges,infEdges);
        aveStabVec0(currIndsAll) = max(aveStabVec0(currIndsAll),rankings(currIndsInd)); % take the max
        % take most extreme partial correlation with caution (to preserve sign change)
        oldNewCoefs = [aveCoefVec0(currIndsAll) coefVec(currIndsInd)];
        [maxes, inds] = max(abs(oldNewCoefs),[],2);
        oldInds = find(inds == 1);
        newInds = find(inds == 2);
        aveCoefVec0(currIndsAll(oldInds)) = oldNewCoefs(oldInds,1);
        aveCoefVec0(currIndsAll(newInds)) = oldNewCoefs(newInds,2); % update new values
        aveInPriorVec0(currIndsAll) = aveInPriorVec0(currIndsAll) + abs(inPriorVec(currIndsInd)); % running sum
        uniStabs = sort(setdiff(unique(rankings),0),'descend');
        totStabs = length(uniStabs);
        totVals = 0;
        for sind = 1:totStabs % from highest stability to lowest
            rankInds = find(rankings(currIndsInd) == uniStabs(sind));
            totVals = totVals + length(rankInds);
            aveQuantVec0(currIndsAll(rankInds)) = max(aveQuantVec0(currIndsAll(rankInds)),1-totVals/totInfInts); % quantile in terms of totInfInts
        end
    end
elseif find(ismember({combineOpt},'mean'))
    disp('mean-combine')
    for nind = 1:totNets
        load(nets2combine{nind})
        infEdges = infEdgeCollections{nind};
        [vals, currIndsAll, currIndsInd] = intersect(allEdges,infEdges);
        aveStabVec0(currIndsAll) = aveStabVec0(currIndsAll) + rankings(currIndsInd); % sum now, average later
        aveCoefVec0(currIndsAll) = aveCoefVec0(currIndsAll) + coefVec(currIndsInd);
        aveInPriorVec0(currIndsAll) = aveInPriorVec0(currIndsAll) + abs(inPriorVec(currIndsInd));
        uniStabs = sort(setdiff(unique(rankings),0),'descend');
        totStabs = length(uniStabs);
        totVals = 0;
        for sind = 1:totStabs % from highest stability to lowest
            rankInds = find(rankings(currIndsInd) == uniStabs(sind));
            totVals = totVals + length(rankInds);
            aveQuantVec0(currIndsAll(rankInds)) = aveQuantVec0(currIndsAll(rankInds)) + 1-totVals/totInfInts;
        end
    end
    aveCoefVec0 = aveCoefVec0 / totNets;
    aveStabVec0 = aveStabVec0 / totNets;
else
    error('combineOpt must be either "mean" or "max".')
end

%% 2. Take top meanEdgesPerGene and recalculate quantiles
[quantsOrd, quantOrdInds] = sort(aveQuantVec0,'descend');
totQuantEdges = totTargs * meanEdgesPerGene;
combQuantiles = zeros(totQuantEdges,1);

if totInfInts > totQuantEdges
    ranks4quant = quantsOrd(1:totQuantEdges); % note there might be stability
    % ties at the end of the ranks4quant matrix
    disp(['Total networks edges (' num2str(totInfInts) ') > meanEdgesPerGene (' num2str(meanEdgesPerGene) ', ' num2str(totQuantEdges) ').']) 
else
    ranks4quant = zeros(totQuantEdges,1);
    ranks4quant(1:totInfInts) = quantsOrd;
    disp(['Total networks edges (' num2str(totInfInts) ') < meanEdgesPerGene (' num2str(meanEdgesPerGene) ', ' num2str(totQuantEdges) ').']) 
end
uniRanks = sort(setdiff(unique(ranks4quant),[0]),'descend');
totRanks = length(uniRanks);
totVals = 0;
for rind = 1:totRanks
    rankInds = find(ranks4quant==uniRanks(rind));
    totVals = totVals + length(rankInds);
    combQuantiles(rankInds) = 1 - totVals/totQuantEdges;
end
totQuantEdges = length(find(combQuantiles));
aveStabVec = aveStabVec0(quantOrdInds);
aveCoefVec = aveCoefVec0(quantOrdInds);
aveInPriorVec = aveInPriorVec0(quantOrdInds);
allEdgesOrd = {allEdges{quantOrdInds}}';

%% 3. Print out interactions
%% get list of edges to go with rankings and output a sparse network file
fout = fopen(combinedNetTsv,'w');
fprintf(fout,'TF\tTarget\tSignedQuantile\tcombStability\tcombPCorr\tstroke\tstroke-width\tstroke-dasharray\n');
% For network visualization in jp_gene_viz (on Github):
% Stroke -- denotes color in R,G,B format? -- current "stroke" only
%   accepts color names, so have capital Stroke until R,G,B is
%   recognized by "stroke"
% stroke-width will go from [1,2] and will be proportional to stability
% stroke-dash-array will incorporate prior information:
% None --> will be solid and is for Prior-supported edges
% 2,2 --> will be for edges Not supported by the prior
minRank = min(aveStabVec); maxRank = max(aveStabVec); rrange = maxRank-minRank;
for ii = 1:totQuantEdges
    intsBit = strsplit(allEdgesOrd{ii},',');
    reg = intsBit{1}; targ = intsBit{2};
    strokeWidth = 1 + (aveStabVec(ii)-minRank)/rrange;  % 
    currPrho = aveCoefVec(ii);
    if currPrho > 0
        strokeVals = cellstr(num2str(round([currPrho*medRed + (1-currPrho)*lightGrey]')));        
    elseif currPrho < 0             
        strokeVals = cellstr(num2str(round([-currPrho*medBlue + (1+currPrho)*lightGrey]')));
    else % currPrho == 0
        strokeVals = cellstr(num2str([lightGrey]'));
        currPrho = 1; % set currPrho = 1, so that edge will appear
    end
    stroke = ['rgb(' strjoin(strokeVals,',') ')'];
    if aveInPriorVec(ii) % solid line for thing in the prior (jp_gene_viz)
        fprintf(fout,[reg '\t' targ '\t' num2str(sign(currPrho)*combQuantiles(ii)) ...
            '\t' num2str(aveStabVec(ii))...
            '\t' strrep(roundstring3(aveCoefVec(ii)),'.0E+00','') '\t' stroke ...
           '\t' num2str(strokeWidth) '\tNone\n']);
    else % dotted line for things not in the prior
        fprintf(fout,[reg '\t' targ '\t' num2str(sign(currPrho)*combQuantiles(ii)) ...
            '\t' num2str(aveStabVec(ii))...
           '\t' strrep(roundstring3(aveCoefVec(ii)),'.0E+00','') '\t' stroke ...
           '\t' num2str(strokeWidth) '\t2,2\n']);
    end
end
fclose(fout);
disp(combinedNetTsv)
end
