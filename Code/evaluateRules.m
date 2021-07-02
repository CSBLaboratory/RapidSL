function constrRxnIDs = evaluateRules(model,geneIDs)
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
