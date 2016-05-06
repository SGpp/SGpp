%% \page example_optimization_m optimization.m
%%
%% On this page, we look at an example application of the sgpp::optimization module.
%% Identical versions of the example are given in all languages
%% currently supported by SG++: C++, Python, Java, and MATLAB.
%%
%% The example interpolates a bivariate test function like the \ref example_tutorial_cpp example.
%% However, we use B-splines here instead to obtain a smoother interpolant.
%% The resulting sparse grid function is then minimized with the method of steepest descent.
%% For comparison, we also minimize the objective function with Nelder-Mead's method.
%%
%% For instructions on how to use SG++ within MATLAB, please see \ref installation.
%% However, for this example to work, some extra steps are necessary.
%% In the following, we assume that we want to run the example on Linux.
%%
%% Please note that in order to get sgpp::optimization to work with MATLAB,
%% you have to disable support for Armadillo and UMFPACK when compiling SG++,
%% i.e. set USE_ARMADILLO and USE_UMFPACK to "no".
%% This is due to incompatible BLAS and LAPACK libraries
%% of Armadillo/UMFPACK and MATLAB
%% (MATLAB uses instead MKL versions of LAPACK and BLAS
%% with different pointer sizes of 64 bits).
%% You can somehow override MATLAB's choice of libraries with
%% the environmental variables BLAS_VERSION and LAPACK_VERSION,
%% but this is strongly discouraged as MATLAB itself may produce
%% unexpected wrong results (e.g., <tt>det [1 2; 3 4] = 2</tt>).
%% Static linking to Armadillo and UMFPACK would be
%% a possible solution to circumvent this problem.
%%
%% To make SG++ usable within MATLAB, please follow
%% \ref linux_using_matlab_jsgpp "the instructions provided here (for Linux)".
%% However, to run the MATLAB example, you additionally need to compile
%% the class \c ExampleFunction into a \c .jar file
%% (before starting MATLAB):
%% \verbatim
%% javac -cp .:/PATH_TO_SGPP/lib/jsgpp/jsgpp.jar ExampleFunction.java
%% jar -cf ExampleFunction.jar ExampleFunction.class
%% \endverbatim
%% Under Windows, note that you have to replace the <tt>:</tt> delimiter in
%% the first line by <tt>;</tt>
%% (<tt>:</tt> is reserved for the use after the drive letter).
%% After setting the environment variables and the
%% \c librarypath.txt correctly
%% (\ref linux_using_matlab_jsgpp "as described here (for Linux)"),
%% you can start MATLAB.
%% In order to run the example, you need to add the example objective
%% function to MATLAB's Java search path via
%% \verbatim
%% javaaddpath('/PATH_TO_SGPP/optimization/examples/ExampleFunction.jar');
%% \endverbatim
%% (i.e., right after adding \c jsgpp.jar to MATLAB's Java path).
%%
%% Now, you should be able to run the MATLAB example,
%% which we describe in the rest of this page.
%% \dontinclude optimization.m
%%
%% At the beginning of the program, we have to load the shared library object file.
%% We can do so by using <tt>sgpp.LoadJSGPPLib.loadJSGPPLib</tt>.
%% Also, we disable OpenMP within jsgpp since it interferes with SWIG's director feature.
sgpp.LoadJSGPPLib.loadJSGPPLib();
sgpp.jsgpp.omp_set_num_threads(1);

fprintf('sgpp::optimization example program started.\n\n');
% increase output verbosity
printer = sgpp.OptPrinter.getInstance();
printer.setVerbosity(2);
printLine = @() fprintf([repmat('-', 1, 80) '\n']);

%% Here, we set define some parameters: objective function, dimensionality,
%% B-spline degree, maximal number of grid points, and adaptivity.
% objective function
f = ExampleFunction();
% dimension of domain
d = f.getNumberOfParameters();
% B-spline degree
p = 3;
% maximal number of grid points
N = 30;
% adaptivity of grid generation
gamma = 0.95;

%% First, we define a grid with modified B-spline basis functions and
%% an iterative grid generator, which can generate the grid adaptively.
grid = sgpp.Grid.createModBsplineGrid(d, p);
gridGen = sgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, gamma);

%% With the iterative grid generator, we generate adaptively a sparse grid.
printLine();
fprintf('Generating grid...\n\n');

if ~gridGen.generate()
    error('Grid generation failed, exiting.');
end

%% Then, we hierarchize the function values to get hierarchical B-spline
%% coefficients of the B-spline sparse grid interpolant
%% \f$\tilde{f}\colon [0, 1]^d \to \mathbb{R}\f$.
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

%% We define the interpolant \f$\tilde{f}\f$ and its gradient
%% \f$\nabla\tilde{f}\f$ for use with the gradient method (steepest descent).
%% Of course, one can also use other optimization algorithms from
%% sgpp::optimization::optimizer.
printLine();
fprintf('Optimizing smooth interpolant...\n\n');
ft = sgpp.OptInterpolantScalarFunction(grid, coeffs);
ftGradient = sgpp.OptInterpolantScalarFunctionGradient(grid, coeffs);
gradientDescent = sgpp.OptGradientDescent(ft, ftGradient);
x0 = sgpp.DataVector(d);

%% The gradient method needs a starting point.
%% We use a point of our adaptively generated sparse grid as starting point.
%% More specifically, we use the point with the smallest
%% (most promising) function value and save it in x0.
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

%% We apply the gradient method and print the results.
gradientDescent.setStartingPoint(x0);
gradientDescent.optimize();
xOpt = gradientDescent.getOptimalPoint();
ftXOpt = gradientDescent.getOptimalValue();
fXOpt = f.eval(xOpt);

fprintf(['\nxOpt = ' char(xOpt.toString()) '\n']);
fprintf(['f(xOpt) = ' num2str(fXOpt, 6) ...
         ', ft(xOpt) = ' num2str(ftXOpt, 6) '\n\n']);

%% For comparison, we apply the classical gradient-free Nelder-Mead method
%% directly to the objective function \f$f\f$.
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

%% The example program outputs the following results:
%% \verbinclude optimization.output.txt
%%
%% We see that both the gradient-based optimization of the smooth sparse grid
%% interpolant and the gradient-free optimization of the objective function
%% find reasonable approximations of the minimum, which lies at
%% \f$(3\pi/16, 3\pi/14) \approx (0.58904862, 0.67319843)\f$.
