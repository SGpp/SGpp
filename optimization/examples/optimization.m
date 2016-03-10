% load jsgpp
sgpp.LoadJSGPPLib.loadJSGPPLib();
% disable OpenMP multi-threading within jsgpp
% (interferes with SWIG's director feature)
sgpp.jsgpp.omp_set_num_threads(1);
% increase output verbosity
printer = sgpp.OptPrinter.getInstance();
printer.setVerbosity(2);
printLine = @() fprintf(['----------------------------------------' ...
                         '----------------------------------------\n']);

fprintf('sgpp::optimization example program started.\n\n');

% objective function
f = ExampleFunction();
d = f.getNumberOfParameters();
p = 3;
N = 30;
gamma = 0.95;

grid = sgpp.Grid.createModBsplineGrid(d, p);
gridGen = sgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, gamma);

%% GRID GENERATION
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printLine();
fprintf('Generating grid...\n\n');

if ~gridGen.generate()
    error('Grid generation failed, exiting.');
end

%% HIERARCHIZATION
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printLine();
fprintf('Hierarchizing...\n\n');
functionValues = gridGen.getFunctionValues();
coeffs = sgpp.DataVector(functionValues.getSize());
hierSLE = sgpp.OptHierarchisationSLE(grid);
sleSolver = sgpp.OptAutoSLESolver();

% solve linear system
if ~sleSolver.solve(hierSLE, functionValues, coeffs)
    error('Solving failed, exiting.');
end

%% OPTIMIZATION OF THE SMOOTH INTERPOLANT
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printLine();
fprintf('Optimizing smooth interpolant...\n\n');
ft = sgpp.OptInterpolantScalarFunction(grid, coeffs);
ftGradient = sgpp.OptInterpolantScalarFunctionGradient(grid, coeffs);
gradientMethod = sgpp.OptGradientDescent(ft, ftGradient);
x0 = sgpp.DataVector(d);

% determine best grid point as starting point
gridStorage = grid.getStorage();

% index of grid point with minimal function value
x0Index = 0;
fX0 = functionValues.get(0);

for i = 1:functionValues.getSize()-1
    if functionValues.get(i) < fX0
        fX0 = functionValues.get(i);
        x0Index = i;
    end
end

for t = 0:d-1
    x0.set(t, gridStorage.get(x0Index).getCoord(t));
end

ftX0 = ft.eval(x0);

fprintf(['x0 = ' char(x0.toString()) '\n']);
fprintf(['f(x0) = ' num2str(fX0, 6) ', ft(x0) = ' num2str(ftX0, 6) '\n\n']);

gradientMethod.setStartingPoint(x0);
gradientMethod.optimize();
xOpt = gradientMethod.getOptimalPoint();
ftXOpt = gradientMethod.getOptimalValue();
fXOpt = f.eval(xOpt);

fprintf(['\nxOpt = ' char(xOpt.toString()) '\n']);
fprintf(['f(xOpt) = ' num2str(fXOpt, 6) ...
         ', ft(xOpt) = ' num2str(ftXOpt, 6) '\n\n']);

%% NELDER-MEAD OPTIMIZATION OF OBJECTIVE FUNCTION
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printLine();
fprintf('Optimizing objective function (for comparison)...\n\n');

nelderMead = sgpp.OptNelderMead(f, 1000);
nelderMead.optimize();
xOptNM = nelderMead.getOptimalPoint();
fXOptNM = nelderMead.getOptimalValue();
ftXOptNM = ft.eval(xOptNM);

fprintf(['\nxOptNM = ' char(xOptNM.toString()) '\n']);
fprintf(['f(xOptNM) = ' num2str(fXOptNM, 6) ...
         ', ft(xOptNM) = ' num2str(ftXOptNM, 6) '\n\n']);

printLine();

fprintf('\nsgpp::optimization example program terminated.\n');

