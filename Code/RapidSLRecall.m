function [LethalSets, LPSolved] = RapidSLRecall...
    (model, maxCardinality, cutOff, eliList, Mode, grRateMS, cplxModel, cplxModelL1Norm)
% RapidSLRecall finds all SLs with a predefined maximum cardinality. To
% perform parallel computaion, the first level of the DFS algorithm is
% implemented separately in "RapidSL" function.
%
% USAGE:
%
%    [LethalSets, LPSolved] = RapidSLRecall...
%       (model, maxCardinality, cutOff, eliList, grRateMS, cplxModel, cplxModelL1Norm)
%
% INPUT:
%    model:           COBRA model structure.
%    maxCardinality:  Maximum desired cardinality of a synthetic lethal set
%    cutOff:          Threshold of lethality.
%    eliList:         The list of reactions or genes that should be
%                     excluded from the analysis. ([] shows that no
%                     reaction or gene is excluded).
%    Mode:            The mode of the lethality analysis: 'Rxn' for
%                     reaction-based analysis and 'Gene' for gene-based
%                     analysis. (Default = 'Rxn').
%    grRateMS:        Maximum growth rate of the mutant strain
%                     (Default: [~, grRateMS] = buildCplexModel(model)).
%    cplxModel:       The CPLEX model of the parent mutant strain.
%    cplxModelL1Norm: The CPLEX model corresponding to the Taxicap Norm
%                     problem. 
%
% OUTPUTS:
%    LethalSets:      List of synthetic lethal targets - cell.
%    LPSolved:        Number of linear programming problems solved.
%
% .. Author:
%       - Mehdi Dehghan Manshadi 06/2018

%% Preparing the Cplex models
cplxModel.Model.lb = model.lb;
cplxModel.Model.ub = model.ub;
cplxModelL1Norm.Model.lhs(end) = grRateMS;
cplxModelL1Norm.Model.lb(1:length(model.rxns)) = model.lb;
cplxModelL1Norm.Model.ub(1:length(model.rxns)) = model.ub;
cplxModelL1Norm.DisplayFunc = [];

%% First Step: Finding the Seed Space
Sol = cplxModelL1Norm.solve();
Flux = Sol.x(1 : length(model.rxns));
Jnz = find(~eq(Flux, 0))';
if strcmp(Mode, 'Rxn')
    Jnz = setdiff(Jnz, eliList);
    targetList = Jnz;
elseif strcmp(Mode, 'Gene')
    [~, Gnz] = find(model.rxnGeneMat(Jnz, :));
    Gnz = setdiff(Gnz, eliList);
    Gnz = unique(Gnz');
    targetList = Gnz;
end

%% Second Step: Searching Within the Seed Space
if ~isempty(targetList)
    [LethalSets, NonLethalSets, GrowthRates, LPSolved1] = SearchWithinSeedSpaceRecall(model, maxCardinality, targetList, cutOff, Mode, cplxModel, grRateMS);
    % Branching of the DFS
    NoCases = zeros(1, maxCardinality - 1);
    if maxCardinality > 1
        for i = 1 : maxCardinality - 1
            if ~isempty(NonLethalSets{i})
                NoCases(i) = length(NonLethalSets{i}(:, 1));
            else
                NoCases(i) = 0;
            end
        end
    end
    if maxCardinality > 1 && sum(NoCases) > 0
        SNoCases = cumsum(NoCases);
        LPSolved2 = zeros(sum(NoCases), 1);
        for p = 1 : sum(NoCases)
            s = find(p <= SNoCases, 1);
            if s == 1
                i = p;
            else
                i = p - SNoCases(s - 1);
            end
            
            if strcmp(Mode, 'Rxn')
                constrainedRxns = NonLethalSets{s}(i, :);
                Set = constrainedRxns;
                eliListRecall = unique([eliList, Jnz]);
            elseif strcmp(Mode, 'Gene')
                constrainedGenes = NonLethalSets{s}(i, :);
                Set = constrainedGenes;
                constrainedRxns = evaluateRules(model, constrainedGenes);
                eliListRecall = unique([eliList, Gnz]);
            end
            
            gr = GrowthRates{s}(i);
            strain = model;
            strain.lb(constrainedRxns) = 0;
            strain.ub(constrainedRxns) = 0;
            [NewLethalSets, LPs] = RapidSLRecall(strain, maxCardinality - s, cutOff, eliListRecall, Mode, gr, cplxModel, cplxModelL1Norm);
            LPSolved2(p) = LPs;
            for q = s + 1 : maxCardinality
                Lethal_p = LethalSets{q};
                Nrows = size(NewLethalSets{q - s});
                Lethal_q = [repmat(Set, [Nrows(1), 1]), NewLethalSets{q - s}];
                LethalSets(q) = {[Lethal_p; Lethal_q]};
            end
        end
    else
        LPSolved2 = 0;
    end
else
    for i = 1 : maxCardinality
        LethalSets(i, 1) = {[]};
    end
    LPSolved2 = 0;
    LPSolved1 = 0;
end
LPSolved = sum(LPSolved2) + LPSolved1 + 1;
end
