% Copyright (C) 2008-today The SG++ project
% This file is part of the SG++ project. For conditions of distribution and
% use, please see the copyright notice provided with SG++ or at
% sgpp.sparsegrids.org
sgpp.LoadJSGPPLib.loadJSGPPLib();

% import all packages
%import sgpp.*;

% or, better, include only the ones needed
import sgpp.DataVector;
import sgpp.GridGenerator;
import sgpp.GridStorage;
import sgpp.Grid;
import sgpp.GridIndex;
import sgpp.jsgpp;
import sgpp.OperationEval;

% define the function f
f = @(x0,x1) (16.0*(x0-1.0)*x0 * (x1-1.0)*x1);

% create a two-dimensional piecewise bilinear grid
dim = 2;
fprintf('dimensionality:         %u\n', dim);
grid = Grid.createLinearGrid(dim);
gridStorage = grid.getStorage();

% create regular grid, level 3
level = 3;
gridGen = grid.getGenerator();
gridGen.regular(level);
fprintf('number of grid points:  %u\n', gridStorage.size());

% create coefficient vector
alpha = DataVector(gridStorage.size());
alpha.setAll(0);
fprintf('length of alpha vector: %u\n', alpha.getSize());

% set function values in alpha
for i = 0:gridStorage.size()-1
    gp = gridStorage.get(i);
    alpha.set(i,f(gp.getCoord(0), gp.getCoord(1)));
end

fprintf('alpha before hierarchization: %s\n', char(alpha.toString()));

% hierarchize
operationHierarchisation = jsgpp.createOperationHierarchisation(grid); 
operationHierarchisation.doHierarchisation(alpha);
fprintf('alpha after hierarchization:  %s\n', char(alpha.toString()));

% evaluate
p = DataVector(dim);
p.set(0, 0.52);
p.set(1, 0.73);
opEval = jsgpp.createOperationEval(grid);

fprintf('u(0.52, 0.73) = %.4f\n', opEval.eval(alpha,p));
