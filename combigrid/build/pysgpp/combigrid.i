// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%ignore sgpp::combigrid::IndexVectorRange::begin;
%ignore sgpp::combigrid::IndexVectorRange::end;
%shared_ptr(sgpp::combigrid::OperationEvalFullGrid);

%include "combigrid/src/sgpp/combigrid/LevelIndexTypes.hpp"
%include "combigrid/src/sgpp/combigrid/HeterogeneousBasis.hpp"
%include "combigrid/src/sgpp/combigrid/FullGrid.hpp"
%include "combigrid/src/sgpp/combigrid/IndexVectorRange.hpp"
%include "combigrid/src/sgpp/combigrid/CombinationGrid.hpp"

%include "combigrid/src/sgpp/combigrid/OperationEvalCombinationGrid.hpp"
%include "combigrid/src/sgpp/combigrid/OperationEvalFullGrid.hpp"
%include "combigrid/src/sgpp/combigrid/OperationPole.hpp"
%include "combigrid/src/sgpp/combigrid/OperationPoleHierarchisationGeneral.hpp"
%include "combigrid/src/sgpp/combigrid/OperationPoleHierarchisationLinear.hpp"
%include "combigrid/src/sgpp/combigrid/OperationPoleNodalisationBspline.hpp"
%include "combigrid/src/sgpp/combigrid/OperationPoleNodalisationLinear.hpp"
%include "combigrid/src/sgpp/combigrid/OperationUPFullGrid.hpp"
%include "combigrid/src/sgpp/combigrid/OperationUPCombinationGrid.hpp"

namespace std {
  %template(BasisVector) vector<sgpp::base::Basis<sgpp::combigrid::level_t, sgpp::combigrid::index_t>*>;
  %template(FullGridVector) vector<sgpp::combigrid::FullGrid>;
  %template(LevelVector) vector<sgpp::combigrid::level_t>;
  %template(LevelVectorVector) vector<sgpp::combigrid::LevelVector>;
  %template(OperationPoleVector) vector<sgpp::combigrid::OperationPole*>;
}
