function [lethalSets, nonLethalSets, growthRates, LPSolved] = SearchWithinSeedSpaceRecall...
    (model, maxCardinality, targetList, cutOff, Mode, cplxModel, grRateMS)
% SearchWithinSeedSpaceRecall finds all synthetic lethal sets with a
% maximum cardinality within the "Seed Space" using FBA. 
%
% USAGE:
%
%    [lethalSets, nonLethalSets, growthRates, LPSolved] = SearchWithinSeedSpaceRecall...
%    (model, maxCardinality, targetList, cutOff, Mode, cplxModel, grRateMS)
%
% INPUT:
%    model:           COBRA model structure.
%    maxCardinality:  Maximum desired cardinality of a synthetic lethal
%                     set.
%    targetList:      The list of targets related to the "Seed Space".
%
% OPTIONAL INPUTS:
%    cutOff:          Threshold of lethality (Default = 0.01 * maximum
%                     growth rate of the wild-type strain).
%    Mode:            The mode of the lethality analysis: 'Rxn' for
%                     reaction-based analysis and 'Gene' for gene-based
%                     analysis. (Default = 'Rxn').
%    cplxModel:       CPLEX model obtained from the COBRA "model"
%                     (Default = buildCplexModel(model)).
%    grRateMS:        Maximum growth rate of the mutant strain
%                     (Default: [~, grRateMS] = buildCplexModel(model)).
%
% OUTPUTS:
%    LethalSets:      List of synthetic lethal sets - cell.
%    NonLethalSets:   List of non-lethal sets - cell.
%    GrowthRates:     Maximum growth rates of the corresponding mutants of
%                     the non-lethal sets.
%    LPSolved:        Number of linear programming problems solved.
%
% .. Author:
%       - Mehdi Dehghan Manshadi 07/2021

if (nargin < 5) || isempty(Mode)
    Mode = 'Rxn';
end

if nargin < 7
    [cplxModel, grRateMS] = buildCplexModel(model);
end

if (nargin < 4) || isempty(cutOff)
    cutOff = 0.01*grRateMS;
end

if strcmp(Mode, 'Rxn')
    nTargets = length(model.rxns);
elseif strcmp(Mode, 'Gene')
    nTargets = length(model.genes);
end

if nTargets < 256
    Combination = nchoosek(uint8(targetList), 1);
elseif nTargets < 65536
    Combination = nchoosek(uint16(targetList), 1);
else
    Combination = nchoosek((targetList), 1);
end
% pre-allocation of the variables:
% lethalSets: set of lethal targets.
% nonLethalSets: set of non-lethal reactions or genes.
% growthRates: maximum growth rates of the mutants.
% LPSolved: number of linear programming problems solved.
lethalSets = cell(maxCardinality, 1);
nonLethalSets = cell(maxCardinality - 1, 1);
growthRates = cell(maxCardinality - 1, 1);
LPSolved = zeros(maxCardinality, 1);
StuStart = 1;
Start = 1;
for Stu = StuStart : min(maxCardinality, length(targetList))
    nTotalPairs = length(Combination(:, 1));
    LPSolved(Stu) = nTotalPairs;
    
    
    % deletion analysis:
    lethalSetIndx = zeros(nTotalPairs, 1);
    isNonLethalSet = zeros(nTotalPairs, 1);
    grRateMS = zeros(nTotalPairs, 1);
    for targetSet = Start : nTotalPairs
        if strcmp(Mode, 'Rxn')
            constrainedRxns = Combination(targetSet, :);
        elseif strcmp(Mode, 'Gene')
            geneSet = Combination(targetSet, :);
            constrainedRxns = evaluateRules(model, geneSet);
        end
        [grRate, ~] = optMod(cplxModel, constrainedRxns, model);
        if grRate < cutOff      % checking if the set is lethal
            lethalSetIndx(targetSet) = 1;
        else                    % checking if the set is not lethal
            isNonLethalSet(targetSet) = 1;
            grRateMS(targetSet) = grRate;
        end
    end
    
    lethal =  Combination(lethalSetIndx == 1, :);    % separating the lethal sets.
    lethalSets(Stu) = {lethal};
    nonLethalCombinations = Combination(isNonLethalSet == 1, :); % separating the non-lethal sets.
    grRateMS = grRateMS(isNonLethalSet == 1);           % storing the maximum growth rates of the mutant strains.
    
    % omiting the list of single lethal sets from the list of candidates:
    if Stu == 1
        targetList = Combination(lethalSetIndx == 0);
        if isempty(targetList)
            break
        end
    end
    
    % storing the non-lethal sets and the corresponding maximum growth
    % rates, excluding the last iteration:
    if Stu ~= maxCardinality
        nonLethalSets(Stu) = {nonLethalCombinations};
        growthRates(Stu) = {grRateMS};
        if length(targetList) < Stu + 1
            break
        end
    end
    
    % clearing these variables from the workspace:
    clear lethal Combination nonLethalCombinations grRateMS
    
    % creating the combinations of the next iteration
    if Stu ~=(maxCardinality)
        if length(targetList) == 1
            if nTargets < 256
                Combination = uint8(targetList);
            elseif nTargets < 65536
                Combination = uint16(targetList);
            else
                Combination = targetList;
            end
        else
            if nTargets < 256
                Combination = nchoosek(uint8(targetList), (Stu + 1));
            elseif nTargets < 65536
                Combination = nchoosek(uint16(targetList), (Stu + 1));
            else
                Combination = nchoosek((targetList), (Stu + 1));
            end
        end
        
        % excluding the combinations which are supersets of the lethal
        % sets:
        for i = 2 : Stu
            lethal = lethalSets{i};
            for j = 1 : length(lethal(:, 1))
                exclude = sum(ismember(Combination, lethal(j,:)), 2) >= i;
                Combination(exclude, :) = [];
            end
            warning off
        end
        if isempty(Combination)
            break
        end
    end
end
LPSolved = sum(LPSolved);
end
