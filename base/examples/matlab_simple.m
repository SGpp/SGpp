%INTERFACEJAVA Summary of this function goes here
%   Detailed explanation goes here

%import all packages
%import sgpp.*;

%or, better, include only the ones needed
import sgpp.DataVector;
import sgpp.GridGenerator;
import sgpp.GridStorage;
import sgpp.Grid;
import sgpp.GridIndex;
import sgpp.jsgpp;
import sgpp.OperationEval;

% define the function f
f = @(x0,x1) (16.0*(x0-1)*x0 * (x1-1)*x1);

sgpp.LoadJSGPPLib.loadJSGPPLib();

% create a two-dimensional piecewise bi- linear grid
dim = 2
grid = Grid.createLinearGrid(dim);
gridStorage = grid.getStorage();
display(gridStorage.dim());

% create regular grid, level 3
level = 3
gridGen = grid.createGridGenerator();
gridGen.regular(level);
Gridlevel = level
print = ['number of grid points = ',num2str(gridStorage.size())];
display(print);

% create coefficient vector
alpha = DataVector(gridStorage.size());
alpha.setAll(0);
print = ['length of alpha vector = ', num2str(alpha.getSize())];
display(print);

% set function values in alpha
for i = 0 : gridStorage.size()-1
    gp = gridStorage.get(i);
    alpha.set(i,f(gp.abs(0),gp.abs(1)));
end

display(alpha.toString());

% hierarchize
operationHierarchisation = jsgpp.createOperationHierarchisation(grid); 
operationHierarchisation.doHierarchisation(alpha);
display(alpha.toString());

% evaluate
p = DataVector(dim);
p.set(0,0.52);
p.set(1,0.73);
opEval = jsgpp.createOperationEval(grid);

print = ['u(0.52,0.73) = ', num2str(opEval.eval(alpha,p))];
display(print);