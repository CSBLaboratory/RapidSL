function constrRxnIDs = evaluateRules(model,geneIDs)
% evaluateRules obtains the affected reactions from GPR rules of COBRA model
%
% USAGE:
%
%    constrRxnIDs = evaluateRules(model,geneIDs)
%
% INPUT:
%    model:           COBRA model structure contaning "rules" field.
%
% OUTPUTS:
%    constrRxnIDs:    The list of constrained reactions.
%
% .. Author:
%       - Mehdi Dehghan Manshadi 07/2021 based on github.com/RamanLab/FastSL
rxnInd = find(any(model.rxnGeneMat(:, geneIDs), 2));
if (~isempty(rxnInd))
    x = true(length(model.genes), 1);
    x(geneIDs) = false;
    constrainRxn = false(length(rxnInd), 1);
    for rxnNo = 1:length(rxnInd)
        if (~eval(model.rules{rxnInd(rxnNo)}))
            constrainRxn(rxnNo) = true;
        end
    end
    if any(constrainRxn)
        constrRxnIDs = rxnInd(constrainRxn);
    else
        constrRxnIDs = [];
    end
else
    constrRxnIDs = [];
end
end
