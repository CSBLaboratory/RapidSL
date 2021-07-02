function [sol, Flux] = optMod(cplxModel, rxnIDs, model, rhs)
% optMod changes the bounds of a reaction list, rxnIDs, to zero and
% finds the optimal solution for the input model. Also it calculates the
% minimum of the Taxicap Norm of the model by specifying the maximum growth
% rate. 
%
% USAGE:
%
%    [sol, solx] = optMod(cplxModel, rxnIDs, model, rhs)
%
% INPUT:
%    cplxModel:       CPLEX model structure obtained by
%                     "buildCplexModel" function.
%    rxnIDs:          Reaction IDs that should be removed from the
%                     cplexModel.
%
% OPTIONAL INPUTS:
%    sol:             Optimal solution for the CPLEX model.
%    Flux:            Obtained fluxes through reactions corresponding to
%                     the optimal solution.
%
% OUTPUTS:
%    model:           COBRA model structure including reaction names
%    rhs:             rhs parameter of the CPLEX model (i.e. the maximum
%                     growth rate of the strain.
%
% .. Author:
%       - Mehdi Dehghan Manshadi 06/2018 based on github.com/RamanLab/FastSL


if exist('rhs', 'var')
    if ~isempty(rhs)
        orig = cplxModel.Model.lhs(end);
        cplxModel.Model.lhs(end) = rhs;
    end
end
cplxModel.Model.lb(rxnIDs) = 0;
cplxModel.Model.ub(rxnIDs) = 0;
cplxModel.Param.simplex.display.Cur = 0;
cplxModel.DisplayFunc = [];

soln  =  cplxModel.solve();
cplxModel.Model.lb(rxnIDs) = model.lb(rxnIDs);
cplxModel.Model.ub(rxnIDs) = model.ub(rxnIDs);

if exist('rhs', 'var')
    if ~isempty(rhs)
        cplxModel.Model.lhs(end) = orig;
    end
end

if soln.status == 0
    sol = 0;
    Flux = [];
    if ~exist('rhs', 'var')
        buildCplexModel(model);
    else
        [~, grRate] = buildCplexModel(model);
        buildCplexModelMinNorm(model, grRate);
    end
else
    sol = model.c'*soln.x(1:length(model.rxns));
    Flux = soln.x(1:length(model.rxns))';
end
end