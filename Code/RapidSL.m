function [LethalSets, LPSolved] = RapidSL...
    (model, maxCardinality, cutOff, eliList, Mode)
% RapidSL constructs the first level of the DFS algorithim for
% Rapid-SL and reports all SLs with a predefined maximum cardinality in the
% corresponding model. To perform parallel computation, the first level of
% the DFS algorithm is treated separately. This function calls the
% "RapidSLRecall" function. 
%
% USAGE:
%
%    [LethalSets, LPSolved] = RapidSL...
%    (model, maxCardinality, cutOff, eliList, Mode)
%
% INPUT:
%    model:           COBRA model structure including reaction names.
%    maxCardinality:  Maximum desired cardinality of a synthetic lethal
%                     set
%
% OPTIONAL INPUTS:
%    cutOff:          Threshold of lethality (Default = 0.01 * Maximum
%                     growth rate of the wild-type strain). Set the value
%                     to [], if it is not specified.
%    eliList:         List of reactions or genes that should be excluded
%                     from the analysis. (Default = [])
%    Mode:            The mode of the lethality analysis: 'Rxn' for
%                     reaction-based analysis and 'Gene' for gene-based
%                     analysis. (Default = 'Rxn')
%
% OUTPUTS:
%    LethalSets:      List of synthetic lethal targets - cell.
%    LPSolved:        Number of linear programing problems solved.
%
% .. Author:
%       - Mehdi Dehghan Manshadi 06/2018

if nargin < 4
    eliListIDs = [];
elseif nargin >= 4
    eliListIDs = findRxnIDs(model, eliList)';
end

if ~exist('Mode', 'var') || isempty(Mode)
    Mode = 'Rxn';    
end
formatSpec = '\nStarting the lethality analysis... \n\n';
fprintf(formatSpec)
%% First Step: Finding the Seed_Space
[cplxModel, grRateWT] = buildCplexModel(model);
[cplxModelL1Norm, Flux] = buildCplexModelMinNorm(model, grRateWT);

if nargin == 2 || isempty(cutOff)
    cutOff = 0.01*grRateWT;
end

Jnz = find(~eq(Flux, 0))';
if strcmp(Mode, 'Rxn')
    Jnz = setdiff(Jnz, eliListIDs);
    targetList = Jnz;
    Gnz = [];
    formatSpec = 'Number of flux-carrying reactions in the WT strain: %.0f... \n\n';
    fprintf(formatSpec, length(Jnz))
elseif strcmp(Mode, 'Gene')
    [~, Gnz] = find(model.rxnGeneMat(Jnz, :));
    Gnz = setdiff(Gnz, eliListIDs);
    Gnz = unique(Gnz');
    targetList = Gnz;
    Jnz = [];
    formatSpec = 'Number of genes that code flux-carrying reactions in the WT strain: %.0f... \n\n';
    fprintf(formatSpec, length(Gnz))
end
%% Second Step: Searching Within the Seed Space
formatSpec = 'Identification of the root nodes... \n';
fprintf(formatSpec)
[LethalSetsIdx, NonLethalSets, GrowthRates, LPSolved1] = SearchWithinSeedSpace(model, maxCardinality, targetList, cutOff, Mode, cplxModel, grRateWT);
%% Repeating Step 1&2 using the DFS algorithm
NoCases = zeros(1, maxCardinality - 1);
for i = 1 : maxCardinality - 1
    if ~isempty(NonLethalSets{i})
        NoCases(i) = length(NonLethalSets{i}(:, 1));
    else
        NoCases(i) = 0;
    end
end

totalNoCases = cumsum(NoCases);
LethalSetsSubBranch = cell(sum(NoCases), 1);
LPSolved2 = zeros(sum(NoCases), 1);
formatSpec = 'Root nodes were identified. The number of root nodes: %.0f... \n\n';
fprintf(formatSpec, sum(NoCases))
formatSpec = 'Brancing is started... \n';
fprintf(formatSpec)
parfor p = 1 : sum(NoCases)
    s = find(p <= totalNoCases, 1);
    if s == 1
        i = p;
    else
        i = p - totalNoCases(s - 1);
    end
    
    if strcmp(Mode, 'Rxn')
        constrainedRxns = NonLethalSets{s}(i, :);
        Set = constrainedRxns;
        eliListRecall = unique([eliListIDs, Jnz]);
    elseif strcmp(Mode, 'Gene')
        constrainedGenes = NonLethalSets{s}(i, :);
        Set = constrainedGenes;
        constrainedRxns = evaluateRules(model, constrainedGenes);
        eliListRecall = unique([eliListIDs, Gnz]);
    end
    
    gr = GrowthRates{s}(i);
    strain = model;
    strain.lb(constrainedRxns) = 0;
    strain.ub(constrainedRxns) = 0;
    [NewLethalSets, LPs] = RapidSLRecall(strain, maxCardinality - s, cutOff, eliListRecall, Mode, gr, cplxModel, cplxModelL1Norm);
    LPSolved2(p) = LPs;
    Lethal = cell(maxCardinality, 1);
    for q = s + 1 : maxCardinality
        Nrows = size(NewLethalSets{q - s});
        Lethal_q = [repmat(Set, [Nrows(1), 1]), NewLethalSets{q - s}];
        Lethal(q) = {Lethal_q};
    end
    LethalSetsSubBranch(p) = {Lethal};
end
formatSpec = 'Brancing is completed... \n';
fprintf(formatSpec)
%% Classifying the results
for p = 1 : sum(NoCases)
    for q = 2 : maxCardinality
        Lethal_q = LethalSetsIdx{q};
        Lethal_p = LethalSetsSubBranch{p}{q};
        LethalSetsIdx(q) = {[Lethal_p; Lethal_q]};
    end
end

for r = 3 : maxCardinality
    Secondary = LethalSetsIdx{r};
    for q = 2 : r - 1
        Primary = LethalSetsIdx{q};
        for p = 1:length(Primary(:, 1))
            exclude = sum(ismember(Secondary, Primary(p, :)), 2) >= q;
            Secondary(exclude,:) = [];
        end
    end
    LethalSetsIdx{r} = Secondary;
end
LPSolved = {(LPSolved1 + 2); LPSolved2};
LethalSets = cell(maxCardinality, 1);
for i = 1 : maxCardinality
    LethalSets{i} = model.rxns(LethalSetsIdx{i});
end
formatSpec = 'Analysis is completed. Total LPs solved: %.0f! \n';
fprintf(formatSpec, sum(LPSolved{2}) + LPSolved{1})
end