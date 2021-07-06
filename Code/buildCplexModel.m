function [cplxModel, grRate] = buildCplexModel(model)
% buildCplexModel constructs the CPLEX model structure from the COBRA model
% structure. 
%
% USAGE:
%
%    [cplxModel, grRate] = buildCplexModel(model)
%
% INPUT:
%    model:           COBRA model structure
%
% OUTPUTS:
%    cplxModel:       CPLEX model structure.
%    grRate:          Maximum growth rate of the strain
%
% .. Author:
%       - Mehdi Dehghan Manshadi 07/2021 based on github.com/RamanLab/FastSL

cplxModel = Cplex();
cplxModel.DisplayFunc = [];
cplxModel.Model.A = sparse(model.S);
cplxModel.Model.obj = model.c;
cplxModel.Model.rhs = model.b;
cplxModel.Model.lhs = model.b;
cplxModel.Model.sense = 'maximize';
cplxModel.Model.lb = model.lb;
cplxModel.Model.ub = model.ub;
cplxModel.Param.barrier.display.Cur = 0;
cplxModel.Param.simplex.display.Cur = 0;
solution = cplxModel.solve();
grRate = solution.objval;
end
