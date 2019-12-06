// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%rename(isEqual) sgpp::combigrid::HeterogeneousBasis::operator==;
%ignore sgpp::combigrid::HeterogeneousBasis::operator!=;
%rename(isEqual) sgpp::combigrid::FullGrid::operator==;
%ignore sgpp::combigrid::FullGrid::operator!=;
%ignore sgpp::combigrid::FullGrid::getLevel;
%ignore sgpp::combigrid::IndexVectorRange::begin;
%ignore sgpp::combigrid::IndexVectorRange::end;

%include "combigrid/src/sgpp/combigrid/LevelIndexTypes.hpp"

%include "combigrid/src/sgpp/combigrid/basis/HeterogeneousBasis.hpp"

%include "combigrid/src/sgpp/combigrid/grid/FullGrid.hpp"
%include "combigrid/src/sgpp/combigrid/grid/IndexVectorRange.hpp"
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

namespace std {
  %template(BasisVector) vector<sgpp::base::Basis<sgpp::combigrid::level_t, sgpp::combigrid::index_t>*>;
  %template(FullGridVector) vector<sgpp::combigrid::FullGrid>;
  %template(LevelVector) vector<sgpp::combigrid::level_t>;
  %template(LevelVectorVector) vector<sgpp::combigrid::LevelVector>;
  %template(OperationPoleVector) vector<sgpp::combigrid::OperationPole*>;
}
