/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "SG++-Doxygen-Documentation", "index.html", [
    [ "SG++: General Sparse Grid Toolbox", "index.html", "index" ],
    [ "Copyright", "copyright.html", null ],
    [ "Developer Manual", "development.html", [
      [ "Overview", "development.html#development_overview", null ],
      [ "Eclipse", "development.html#development_eclipse", [
        [ "Setup", "development.html#development_eclipse_setup", null ],
        [ "Fix C++11 header indexing", "development.html#development_eclipse_header", null ],
        [ "Configure SCons build", "development.html#development_eclipse_scons", [
          [ "Overwrite Eclipse Build", "development.html#development_eclipse_scons_overwrite", null ],
          [ "Install SConsolidator", "development.html#development_eclipse_scons_sconsolidator", null ]
        ] ]
      ] ],
      [ "Testing", "development.html#development_testing", null ],
      [ "Documentations, Styleguide, and Doxygen", "development.html#development_doxygen", [
        [ "Usage", "development.html#development_doxygen_usage", null ]
      ] ],
      [ "Coding", "development.html#development_coding", [
        [ "Comments", "development.html#development_coding_comments", null ],
        [ "Style", "development.html#development_coding_style", null ],
        [ "Naming", "development.html#development_coding_naming", null ],
        [ "Don'ts", "development.html#development_coding_donts", null ]
      ] ],
      [ "Creating new modules", "development.html#development_newmodules", null ]
    ] ],
    [ "Usage Examples", "examples.html", [
      [ "C++ Examples", "examples_cpp.html", [
        [ "Module sgpp::base", "examples_cpp.html#examples_cpp_module_base", null ],
        [ "Module sgpp::combigrid", "examples_cpp.html#examples_cpp_module_combigrid", null ],
        [ "Module sgpp::datadriven", "examples_cpp.html#examples_cpp_module_datadriven", null ],
        [ "Module sgpp::optimization", "examples_cpp.html#examples_cpp_module_optimization", null ],
        [ "Module sgpp::solver", "examples_cpp.html#examples_cpp_module_solver", null ],
        [ "benchmark_gridInteraction.cpp", "example_benchmark_gridInteraction_cpp.html", null ],
        [ "Using the DataMatrix object", "example_dataMatrixSerializeDemo_cpp.html", null ],
        [ "Using the DataVector object", "example_dataVectorSerializeDemo_cpp.html", null ],
        [ "Detect the configuration of OpenCL platforms", "example_detectPlatformConfiguration_cpp.html", null ],
        [ "Interaction-Term aware sparse grids.", "example_gridInteractionExample_cpp.html", null ],
        [ "Generalised Sparse Grids", "example_gridTExample_cpp.html", null ],
        [ "Using JSON", "example_json_cpp.html", null ],
        [ "Spatially-Dimension-Adaptive Refinement in C++", "example_predictiveRefinement_cpp.html", null ],
        [ "Quadrature in C++", "example_quadrature_cpp.html", null ],
        [ "Refinement Example", "example_refinement_cpp.html", null ],
        [ "tutorial.cpp (Start Here)", "example_tutorial_cpp.html", null ],
        [ "Grid unserialization", "example_unserializeGrid_cpp.html", null ],
        [ "List of different Grid Types", "GridTypes.html", null ],
        [ "Combigrid Example Dimensional Adaptivity (C++)", "example_combigrid_adaptive_cpp.html", null ],
        [ "Combigrid Example (C++)", "example_combigrid_cpp.html", null ],
        [ "examplePCE.cpp", "example_examplePCE_cpp.html", null ],
        [ "Learner Classification Test", "example_learnerClassificationTest_cpp.html", null ],
        [ "Regression Learner", "example_learnerRegressionTest_cpp.html", null ],
        [ "Learner SGDE OnOff", "example_learnerSGDEOnOffTest_cpp.html", null ],
        [ "learner SGDE", "example_learnerSGDETest_cpp.html", null ],
        [ "Learner SGD", "example_learnerSGDTest_cpp.html", null ],
        [ "new_sgde.cpp", "example_new_sgde_cpp.html", null ],
        [ "optimize_kde_bandwidth.cpp", "example_optimize_kde_bandwidth_cpp.html", null ],
        [ "Constrained Optimization", "example_constrainedOptimization_cpp.html", null ],
        [ "Fuzzy Extension Principle (C++)", "example_fuzzy_cpp.html", null ],
        [ "Optimization Example (C++)", "example_optimization_cpp.html", null ],
        [ "FISTA Solver", "example_fistaExample_cpp.html", null ]
      ] ],
      [ "Python Examples", "examples_py.html", [
        [ "Module sgpp::base", "examples_py.html#examples_py_module_base", null ],
        [ "Module sgpp::combigrid", "examples_py.html#examples_py_module_combigrid", null ],
        [ "Module sgpp::datadriven", "examples_py.html#examples_py_module_datadriven", null ],
        [ "Module sgpp::optimization", "examples_py.html#examples_py_module_optimization", null ],
        [ "Module sgpp::pde", "examples_py.html#examples_py_module_pde", null ],
        [ "Using the DataMatrix object", "example_dataMatrixSerializeDemo_py.html", null ],
        [ "Using the DataVector object", "example_dataVectorSerializeDemo_py.html", null ],
        [ "Generalised Sparse Grids", "example_gridTExample_py.html", null ],
        [ "Spatially-Dimension-Adaptive Refinement of ANOVA Components in Python", "example_predictiveANOVARefinement_py.html", null ],
        [ "Spatially-Dimension-Adaptive Refinement in Python", "example_predictiveRefinement_py.html", null ],
        [ "Quadrature in Python", "example_quadrature_py.html", null ],
        [ "refinement.py", "example_refinement_py.html", null ],
        [ "Dimension-Adaptive Refinement in Python", "example_subspaceRefinement_py.html", null ],
        [ "tutorial.py (Start Here)", "example_tutorial_py.html", null ],
        [ "Combigrid Example Dimensional Adaptivity (Python)", "example_combigrid_adaptive_py.html", null ],
        [ "Combigrid Example (Python)", "example_combigrid_py.html", null ],
        [ "Generalised Sparse Grids", "example_generalisedGridsest_py.html", null ],
        [ "learnerExample.py", "example_learnerExample_py.html", null ],
        [ "learnerSGDETest.py", "example_learnerSGDETest_py.html", null ],
        [ "positive_density.py", "example_positive_density_py.html", null ],
        [ "test_Rosenblatt.py", "example_test_Rosenblatt_py.html", null ],
        [ "Fuzzy Extension Principle (Python)", "example_fuzzy_py.html", null ],
        [ "Optimization Example (Python)", "example_optimization_py.html", null ],
        [ "splineResponseSurface_example.py", "example_splineResponseSurface_example_py.html", null ],
        [ "LTwoDotTest.py", "example_LTwoDotTest_py.html", null ]
      ] ],
      [ "Java Examples", "examples_java.html", [
        [ "Module sgpp::base", "examples_java.html#examples_java_module_base", null ],
        [ "Module sgpp::datadriven", "examples_java.html#examples_java_module_datadriven", null ],
        [ "Module sgpp::optimization", "examples_java.html#examples_java_module_optimization", null ],
        [ "Refinement Example", "example_refinement_java.html", null ],
        [ "tutorial.java (Start Here)", "example_tutorial_java.html", null ],
        [ "Learner SGDE", "example_example_learnerSGDE_java.html", null ],
        [ "Optimization Example (Java)", "example_optimization_java.html", null ]
      ] ],
      [ "MATLAB Examples", "examples_m.html", [
        [ "Module sgpp::base", "examples_m.html#examples_m_module_base", null ],
        [ "Module sgpp::optimization", "examples_m.html#examples_m_module_optimization", null ],
        [ "tutorial.m (Start Here)", "example_tutorial_m.html", null ],
        [ "Optimization Example (MATLAB)", "example_optimization_m.html", null ]
      ] ]
    ] ],
    [ "Integrate Dakota", "install_dakota.html", null ],
    [ "SGDE Miner", "example_SGDEMinerFromConfigFile_py.html", null ],
    [ "Deprecated List", "deprecated.html", null ],
    [ "Todo List", "todo.html", null ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", "namespacemembers_dup" ],
        [ "Functions", "namespacemembers_func.html", "namespacemembers_func" ],
        [ "Variables", "namespacemembers_vars.html", "namespacemembers_vars" ],
        [ "Typedefs", "namespacemembers_type.html", null ],
        [ "Enumerations", "namespacemembers_enum.html", null ],
        [ "Enumerator", "namespacemembers_eval.html", null ]
      ] ]
    ] ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", "functions_vars" ],
        [ "Typedefs", "functions_type.html", null ],
        [ "Enumerations", "functions_enum.html", null ],
        [ "Enumerator", "functions_eval.html", null ],
        [ "Related Symbols", "functions_rela.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Functions", "globals_func.html", null ],
        [ "Variables", "globals_vars.html", null ],
        [ "Typedefs", "globals_type.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"ANOVAHashRefinement_8cpp.html",
"DBMatDMSOrthoAdapt_8cpp.html",
"DensityEstimationConfiguration_8hpp.html",
"GridFactory_8cpp.html",
"LearnedKnowledgeFormatter_8py.html",
"NaturalBsplineBasis_8hpp.html#a7387c9b142ebc897c68f5835d17af91f",
"OperationEvalModPolyClenshawCurtisNaive_8hpp.html",
"OperationLaplacePolyClenshawCurtisBoundary_8hpp.html",
"OperationMultipleEvalSubspaceSimple__multImpl_8cpp.html",
"OperationUPCombinationGrid_8cpp.html",
"ScalarFunctionHessian_8hpp.html",
"WaveletBoundaryGrid_8hpp.html",
"classjson_1_1IDNode.html#af3978bcf51ac3d7f1f0ec7f875b9d214",
"classpython_1_1controller_1_1CheckpointController_1_1CheckpointController.html#a127c368c750f80a2035cc35f26753d99",
"classpython_1_1learner_1_1LearnerBuilder_1_1LearnerBuilder.html#a932bad4e4782258d3654b5e29a8d54ce",
"classpython_1_1learner_1_1Types_1_1SolverTypes.html",
"classpython_1_1uq_1_1analysis_1_1asgc_1_1ASGCAnalysisBuilder_1_1ASGCAnalysisBuilder.html#a7bc927811550c1d10a076a5e9d5cc362",
"classpython_1_1uq_1_1dists_1_1Beta_1_1Beta.html#a7ee541d8f7ace0b78789b92ad0309f06",
"classpython_1_1uq_1_1dists_1_1LibAGFDist_1_1LibAGFDist.html#a6302f0850bc63e15ae4a4510727765bf",
"classpython_1_1uq_1_1dists_1_1TLognormal_1_1TLognormal.html#a79ac9052804c019aee63caf6623aad2b",
"classpython_1_1uq_1_1learner_1_1Learner_1_1Learner.html#a7761430b699d4d07e2728fbc96645f33",
"classpython_1_1uq_1_1learner_1_1builder_1_1RegressorSpecificationDescriptor_1_1FoldingDescriptor.html#a9ff9529d2a6a5e79a6ea11a74e91df02",
"classpython_1_1uq_1_1models_1_1Model_1_1Model.html#a8487f7ac4280cdd263ac7f6b57b99ac2",
"classpython_1_1uq_1_1operations_1_1forcePositivity_1_1localFullGridSearch_1_1LocalFullGridCandidates.html#af4478f793cef2e545e3fa90eebf902e1",
"classpython_1_1uq_1_1parameters_1_1ParameterBuilder_1_1ParameterBuilder.html#ad9ab35ce196929c9f5b19c0916e09d22",
"classpython_1_1uq_1_1quadrature_1_1HashQuadrature_1_1HashQuadratureMap.html#a87f51965f6b806e5a757bde31aef08c8",
"classpython_1_1uq_1_1refinement_1_1RefinementManagerDescriptor_1_1RefineCurrentNodesDescriptor.html#aef00d3970617b08be0216cde341c9c4d",
"classpython_1_1uq_1_1refinement_1_1RefinementStrategy_1_1VarianceOptRanking.html#a57bb5d37fca845a5b5b7e0021d2842ca",
"classpython_1_1uq_1_1sampler_1_1asgc_1_1ASGCSamplerSpecification_1_1ASGCSamplerSpecification.html#a99857fa7617ae1301dce6c7788a7b6ca",
"classpython_1_1uq_1_1uq__setting_1_1UQSettingManager_1_1UQSettingManager.html#a44097df8976f47f7ef3f48d7ba24a12d",
"classpython_1_1utils_1_1GzipSerializer_1_1GzipSerializer.html",
"classsgpp_1_1base_1_1BoundingBox.html#a7bdddf4d393e2dd3dc4e5eebfa772a07",
"classsgpp_1_1base_1_1BsplineModifiedClenshawCurtisBasis.html#a656e9bbabd420169d49d18e26602b7d2",
"classsgpp_1_1base_1_1DataMatrix.html#a086d1900db9169266fc337a1d2e059e4",
"classsgpp_1_1base_1_1DataMatrixSP.html#af0116d1f8c411f92582ff940b9b9f8ef",
"classsgpp_1_1base_1_1DehierarchisationFundamentalNakSplineBoundary.html#a07a31a823ec693f452ab2e8e64009767",
"classsgpp_1_1base_1_1DehierarchisationPolyBoundary.html#a97378713bde24873da6d6324e899520a",
"classsgpp_1_1base_1_1EvalCuboidGenerator.html",
"classsgpp_1_1base_1_1FundamentalSplineModifiedBasis.html#a40c29c01a6a37edc54a20ea76d584dd0",
"classsgpp_1_1base_1_1Grid.html#afc2b2af07299be3aec944711f2decc45",
"classsgpp_1_1base_1_1HashGridIterator.html#acc5fe8ca780c0bdef4599c85a1849347",
"classsgpp_1_1base_1_1HashRefinement.html#aafeea5daa612c1722d9a5ba7dc06f5ec",
"classsgpp_1_1base_1_1HierarchisationModLinearClenshawCurtis.html#af35815f38fbf9b935ea633257aa8e4e0",
"classsgpp_1_1base_1_1HierarchisationSLE.html#ae9e76db3afe531012f14d63a4321e202",
"classsgpp_1_1base_1_1L0BoundaryGridGenerator.html#a5af933dbf84b080376ded6169f2555e6",
"classsgpp_1_1base_1_1LinearPeriodicBasis.html",
"classsgpp_1_1base_1_1ModPolyClenshawCurtisGrid.html#a593a0ec2dd0e1ee3e8b51a76b6d8e889",
"classsgpp_1_1base_1_1NakBsplineExtendedBasis.html#abe5257fa78947645ef68e27d060a2e19",
"classsgpp_1_1base_1_1OCLClonedBuffer.html#af3f2f1634b4a76cda01bf2e19c19ee79",
"classsgpp_1_1base_1_1OperationConvertPrewavelet.html#a94b03ca5e6ccde2f5709ca8c16595931",
"classsgpp_1_1base_1_1OperationEvalGradientModFundamentalSplineNaive.html#a49e9d2c663073f29c5575b2c265b301c",
"classsgpp_1_1base_1_1OperationEvalHessianBsplineBoundaryNaive.html#a0fb94435901d1cb7f49299145efd20f7",
"classsgpp_1_1base_1_1OperationEvalHessianWaveletNaive.html#a787162e0342edb1f7f1e532a5446b35e",
"classsgpp_1_1base_1_1OperationEvalModNakBsplineNaive.html",
"classsgpp_1_1base_1_1OperationEvalPartialDerivativeFundamentalNakSplineNaive.html#a0a9f0d5750b2167c5b20c41bcecb9122",
"classsgpp_1_1base_1_1OperationEvalPolyBoundaryNaive.html",
"classsgpp_1_1base_1_1OperationFirstMomentModPolyClenshawCurtis.html#a40f16c4952b2901600cd0658105e0dd8",
"classsgpp_1_1base_1_1OperationHierarchisationModPolyClenshawCurtis.html#a27d831ca43e0b5618cee4503439ec67b",
"classsgpp_1_1base_1_1OperationMultipleEvalLinearBoundaryNaive.html#a3b5ab31abe3033229360f25f57cab2ec",
"classsgpp_1_1base_1_1OperationMultipleEvalNakBsplineModifiedNaive.html#ad312088c3ea845cfd828d69394772c2c",
"classsgpp_1_1base_1_1OperationQuadratureMC.html#a9147ab68a9d010bbf9d91e238ebd7a77",
"classsgpp_1_1base_1_1OperationSecondMomentModBspline.html#ad3b22881ee89972fea80fa905d2c8cc0",
"classsgpp_1_1base_1_1OperationWeightedQuadratureNakPBspline.html#a6dab8854635d98dad457c678630f03f8",
"classsgpp_1_1base_1_1PolyGrid.html#a42af2e3a73a4011932f25e4edff1abc1",
"classsgpp_1_1base_1_1Printer.html#ac5749a4b8896907c615c0c170522c34b",
"classsgpp_1_1base_1_1ScaledScalarFunctionGradient.html#a9a3103cb9e17c065622d49b33d595e62",
"classsgpp_1_1base_1_1StencilHierarchisationModLinear.html#accbad2bbb4c73f3ddd5e9662f9f3ceba",
"classsgpp_1_1base_1_1VectorFunctionGradient.html#ad6c382ed5a9c2b75864666e432128089",
"classsgpp_1_1base_1_1WeaklyFundamentalNakSplineModifiedBasisDeriv2.html#abfe9ac3c9758749e63fbd277edce329c",
"classsgpp_1_1base_1_1operation__exception.html",
"classsgpp_1_1combigrid_1_1FullGrid.html#a48d6ddb7cf9f444df565f4d0126e7337",
"classsgpp_1_1combigrid_1_1OperationPoleHierarchisationGeneral_1_1HierarchisationGeneralSLE.html",
"classsgpp_1_1datadriven_1_1AlgorithmAdaBoostBase.html#af1f54362d8b5876ac6d931c3a2726945",
"classsgpp_1_1datadriven_1_1ClassificationRefinementFunctor.html#a89de666a38ced5191d0076f08ca476f7",
"classsgpp_1_1datadriven_1_1DBMatOffline.html#a9afa11ec1544b947dc2de033468bd25f",
"classsgpp_1_1datadriven_1_1DBMatOnlineDE.html#aa071b7b9a2b9437ed66d98c5ab5d5103",
"classsgpp_1_1datadriven_1_1DataBasedRefinementFunctor.html#a3a7f376999d8812eacd89459bce56647",
"classsgpp_1_1datadriven_1_1DataShufflingFunctorSequential.html",
"classsgpp_1_1datadriven_1_1DensityEstimationMinerFactory.html#a1cb56aa8f157b11fa8cc2fe56a4261bb",
"classsgpp_1_1datadriven_1_1EpanechnikovKernel.html#af0da6218a2e57820fa6a2502cd34c708",
"classsgpp_1_1datadriven_1_1GridPointBasedRefinementFunctor.html#a87f2ebc16d85fca3d910d0c2936c9a6b",
"classsgpp_1_1datadriven_1_1LearnerBase.html#a244802c78e9ff5ad3a221f6b3d8187b6",
"classsgpp_1_1datadriven_1_1LearnerSGDE.html#a485e7aadf4f92e2ea707f0444b33e98c",
"classsgpp_1_1datadriven_1_1LearnerSVM.html#afd342bed9fae503aaef01baa6644424e",
"classsgpp_1_1datadriven_1_1ModelFittingBase.html#ad3ca8a1b41bc2dbb01c173ee7fedebde",
"classsgpp_1_1datadriven_1_1ModelFittingDensityEstimationCombi.html#a50f68e3b938af9d429a5e2c260a306de",
"classsgpp_1_1datadriven_1_1NearestNeighbors.html#ac148ec80d11919488b46415b169709ae",
"classsgpp_1_1datadriven_1_1OperationInverseRosenblattTransformation1DModBsplineClenshawCurtis.html#aea621aea593bc085fd00a254788b56ee",
"classsgpp_1_1datadriven_1_1OperationInverseRosenblattTransformationPolyBoundary.html#adc08fc93d8eda0bc352081777347fe81",
"classsgpp_1_1datadriven_1_1OperationMultiEvalCuda.html#a7ddafc2663c68ab5fd49fd8f2c8580e4",
"classsgpp_1_1datadriven_1_1OperationMultiEvalStreamingModOCLFastMultiPlatform.html#ab299305c00cb7316d32e07bfa851139c",
"classsgpp_1_1datadriven_1_1OperationMultipleEvalMatrix.html",
"classsgpp_1_1datadriven_1_1OperationRosenblattTransformationBspline.html#af41b9f08f1bacbc450088b58d53b34ed",
"classsgpp_1_1datadriven_1_1OperationTestLinearBoundary.html#a6dd8237dae84b5f2f33a5ca7b066af43",
"classsgpp_1_1datadriven_1_1PrimalDualSVM.html#a5a870ec5f803c18b22ffa018fc27f2f3",
"classsgpp_1_1datadriven_1_1SortedDataset.html#a363a046ae7c5ca5cef8c0d2ace982d53",
"classsgpp_1_1datadriven_1_1StreamingModOCLFastMultiPlatform_1_1KernelMultTranspose.html#ab840b52b197ff3fa815a49d7cb76edb9",
"classsgpp_1_1datadriven_1_1SubspaceNodeCombined.html#a1fad7f511e02aaf95b7a9d3b642af42cacb4fb1757fb37c43cded35d3eb857c43",
"classsgpp_1_1datadriven_1_1VisualizerDummy.html#a45a68a9053cb09aef2cf22e2a11ba6b8",
"classsgpp_1_1datadriven_1_1clusteringmpi_1_1OperationDummy.html#a790fa0f255f3869e94ad13a5535b4b1b",
"classsgpp_1_1optimization_1_1FuzzyInterval.html#a24f731fd570b51007efed161136bfd76",
"classsgpp_1_1optimization_1_1OperationMultipleHierarchisation.html#a6afa33944368945c5d54cc5893583747",
"classsgpp_1_1optimization_1_1OperationMultipleHierarchisationModNakBspline.html#ae947535ab9cf184cc8b918ec19635161",
"classsgpp_1_1optimization_1_1ResponseSurfaceVector.html#af5d77bbe236ae1fea1bf34db87f78587",
"classsgpp_1_1optimization_1_1optimizer_1_1AugmentedLagrangian.html#aac77b2f28f60410af78e50871b3c9417",
"classsgpp_1_1optimization_1_1optimizer_1_1LevenbergMarquardt.html#a2de51102b3d2b19865c7e913b62f365c",
"classsgpp_1_1optimization_1_1optimizer_1_1Rprop.html#a12334548098540574c5a0ba0caeb92eb",
"classsgpp_1_1optimization_1_1test__problems_1_1Branin02.html",
"classsgpp_1_1optimization_1_1test__problems_1_1G04EqualityConstraint.html#ac57bd7145626bb7f7df7b56bf422fd31",
"classsgpp_1_1optimization_1_1test__problems_1_1G09EqualityConstraint.html#a8f2c4e4d3150011256400df053666d61",
"classsgpp_1_1optimization_1_1test__problems_1_1Griewank.html#a354337235045fb7e66222a1a468f2b6e",
"classsgpp_1_1optimization_1_1test__problems_1_1SHCB.html#a7e5580adeb572bde1b9ab9c71e1f5fda",
"classsgpp_1_1optimization_1_1test__problems_1_1TremblingParabolaObjective.html",
"classsgpp_1_1pde_1_1LaplaceEnhancedDownBBLinearBoundary.html",
"classsgpp_1_1pde_1_1OperationLaplaceExplicitLinear.html#a3e0c747bedf7695a51dd94654b92f4d1",
"classsgpp_1_1pde_1_1OperationMatrixLTwoDotExplicitLinear.html#a75d2c7f6b1f16590a010d736559c9c25",
"classsgpp_1_1pde_1_1OperationParabolicPDESolverSystemDirichlet.html#a341c2e66018a6b3b6421dabed311d655",
"classsgpp_1_1pde_1_1PoissonEquationSolver.html#a1546ee7b0507268a297acf29ba61433c",
"classsgpp_1_1pde_1_1UpDownOneOpDimWithShadow.html#a78fcc3547397ed3a9cc6995012fb3e69",
"classsgpp_1_1solver_1_1ConjugateGradientsSP.html#ae69e87b4b2005fd29aadc3ad7f7de669",
"classsgpp_1_1solver_1_1StepsizeControl.html#ac8c613e25e1a88cdf66e266797ac5302",
"dir_05bb2443d48d50300800ebf900242aa1.html",
"example__dataMatrixSerializeDemo__cpp_8doxy.html",
"functions_vars_h.html",
"namespacegenerate__datasets.html",
"namespacepython_1_1plotDeltas3d.html#af552b24d750d2a6063e318673f34af6a",
"namespacepython_1_1uq_1_1operations_1_1sparse__grid.html#a0f62fdfc8463e1f7936edaaf075e3d5b",
"namespacepython_1_1utils_1_1data__projections.html#adcdaef4f46a30a07c7e83826ec25e941",
"namespacesgpp_1_1datadriven.html#a265b15991588657514902559d1cdc555",
"namespacesgpp_1_1optimization_1_1file__io.html#a77813e2ef46cb457bcb309c5483ffcb0",
"sg__projections_8py.html#a1435dce1b5ac53cacee5dc44409586db",
"structsgpp_1_1datadriven_1_1CrossvalidationConfiguration.html#ac6f48ebd334333dfe79b99fbb30734da",
"structsgpp_1_1datadriven_1_1RefinementResultNetworkMessage.html#a896422e8da0fffcf371c36bf654c27dd"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';