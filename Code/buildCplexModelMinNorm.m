function [cplxModel, Flux] = buildCplexModelMinNorm(model, grRate)
% buildCplexModelMinNorm constructs the CPLEX model structure from the
% COBRA model to calculate the minimum of the Taxicab Norm of fluxes
% corresponding to the maximum growth rate of the strain.
%
% USAGE:
%
%    [cplxModel, Flux] = buildCplexModelMinNorm(model, grRate)
%
% INPUT:
%    model:           COBRA model structure
%
% OUTPUTS:
%    cplxModel:       CPLEX model structure.
%    Flux:            Obtained fluxes through reactions corresponding to
%                     the optimal solution.
%
% .. Author:
%       - Mehdi Dehghan Manshadi 07/2021 based on github.com/RamanLab/FastSL

[nMets, nRxns] = size(model.S);
LPproblem2.A = [model.S sparse(nMets, 2*nRxns);
    speye(nRxns, nRxns) speye(nRxns, nRxns) sparse(nRxns, nRxns);
    -speye(nRxns, nRxns) sparse(nRxns, nRxns) speye(nRxns, nRxns);
    model.c' sparse(1, 2*nRxns)];
LPproblem2.c = [zeros(nRxns, 1); ones(2*nRxns, 1)];
LPproblem2.lb = [model.lb; zeros(2*nRxns, 1)];
LPproblem2.ub = [model.ub; 10000*ones(2*nRxns, 1)];
LPproblem2.b = [model.b; zeros(2*nRxns, 1); grRate];

cplxModel = Cplex();
cplxModel.Model.A = sparse(LPproblem2.A);
cplxModel.Model.obj = LPproblem2.c;

cplxModel.Model.rhs(1:nMets) = LPproblem2.b(1:nMets);
cplxModel.Model.rhs((nMets+1):(nMets+2*nRxns+1)) = Inf*((nMets+1):(nMets+2*nRxns+1));

cplxModel.Model.lhs(1:nMets) = LPproblem2.b(1:nMets);
cplxModel.Model.lhs((nMets+1):(nMets+2*nRxns+1)) = LPproblem2.b((nMets+1):(nMets+2*nRxns+1));

cplxModel.Model.lb = LPproblem2.lb;
cplxModel.Model.ub = LPproblem2.ub;
cplxModel.Model.sense = 'minimize';
cplxModel.DisplayFunc = [];
cplxModel.Param.simplex.display.Cur = 0;
solution = cplxModel.solve();
Flux = solution.x(1:nRxns);
end
