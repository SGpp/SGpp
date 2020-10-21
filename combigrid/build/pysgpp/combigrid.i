// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%ignore sgpp::combigrid::IndexVectorRange::begin;
%ignore sgpp::combigrid::IndexVectorRange::end;
%ignore sgpp::combigrid::OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE;
%shared_ptr(sgpp::combigrid::OperationEvalFullGrid);

%include "std_unique_ptr.i"
wrap_unique_ptr(RelevanceCalculatorUniquePtr, sgpp::combigrid::RelevanceCalculator);
wrap_unique_ptr(PriorityEstimatorUniquePtr, sgpp::combigrid::PriorityEstimator);

%include "combigrid/src/sgpp/combigrid/LevelIndexTypes.hpp"

%include "combigrid/src/sgpp/combigrid/basis/HeterogeneousBasis.hpp"

%include "combigrid/src/sgpp/combigrid/grid/FullGrid.hpp"
%include "combigrid/src/sgpp/combigrid/grid/CombinationGrid.hpp"

%include "combigrid/src/sgpp/combigrid/operation/OperationEvalCombinationGrid.hpp"
%include "combigrid/src/sgpp/combigrid/operation/OperationEvalFullGrid.hpp"
%include "combigrid/src/sgpp/combigrid/operation/OperationPole.hpp"
%include "combigrid/src/sgpp/combigrid/operation/OperationPoleHierarchisationGeneral.hpp"
%include "combigrid/src/sgpp/combigrid/operation/OperationPoleHierarchisationLinear.hpp"
%include "combigrid/src/sgpp/combigrid/operation/OperationPoleNodalisationBspline.hpp"
%include "combigrid/src/sgpp/combigrid/operation/OperationPoleNodalisationLinear.hpp"
%include "combigrid/src/sgpp/combigrid/operation/OperationUPFullGrid.hpp"
%include "combigrid/src/sgpp/combigrid/operation/OperationUPCombinationGrid.hpp"

%include "combigrid/src/sgpp/combigrid/tools/IndexVectorRange.hpp"
%include "combigrid/src/sgpp/combigrid/tools/LevelVectorTools.hpp"

%include "combigrid/src/sgpp/combigrid/adaptive/PriorityEstimator.hpp"
%include "combigrid/src/sgpp/combigrid/adaptive/AveragingPriorityEstimator.hpp"
// an alias to use maps used for calculating priorities
%template(map_levelvector_real) std::map<sgpp::combigrid::LevelVector, double>;
%include "combigrid/src/sgpp/combigrid/adaptive/RelevanceCalculator.hpp"
%include "combigrid/src/sgpp/combigrid/adaptive/WeightedRelevanceCalculator.hpp"
// this doesn't work due to insufficient unique_ptr support in SWIG
// %include "combigrid/src/sgpp/combigrid/adaptive/AdaptiveCombinationGridGenerator.hpp"

namespace std {
  %template(BasisVector) vector<sgpp::base::Basis<sgpp::combigrid::level_t, sgpp::combigrid::index_t>*>;
  %template(FullGridVector) vector<sgpp::combigrid::FullGrid>;
  %template(LevelVector) vector<sgpp::combigrid::level_t>;
  %template(LevelVectorVector) vector<sgpp::combigrid::LevelVector>;
  %template(OperationPoleVector) vector<sgpp::combigrid::OperationPole*>;
}
