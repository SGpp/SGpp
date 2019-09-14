// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%ignore sgpp::combigrid::HeterogeneousBasis::HeterogeneousBasis(size_t, base::Basis<level_t, index_t>&);
%ignore sgpp::combigrid::IndexVectorRange::begin;
%ignore sgpp::combigrid::IndexVectorRange::end;
%ignore sgpp::combigrid::OperationUPCombinationGrid::OperationUPCombinationGrid(const CombinationGrid&, OperationPole&);
%shared_ptr(sgpp::combigrid::OperationEvalFullGrid);

%include "combigrid/src/sgpp/combigrid/LevelIndexTypes.hpp"
%include "combigrid/src/sgpp/combigrid/HeterogeneousBasis.hpp"
%include "combigrid/src/sgpp/combigrid/FullGrid.hpp"
%include "combigrid/src/sgpp/combigrid/IndexVectorRange.hpp"
%include "combigrid/src/sgpp/combigrid/CombinationGrid.hpp"

%include "combigrid/src/sgpp/combigrid/OperationEvalCombinationGrid.hpp"
%include "combigrid/src/sgpp/combigrid/OperationEvalFullGrid.hpp"
%include "combigrid/src/sgpp/combigrid/GeneralOperation.hpp"
%include "combigrid/src/sgpp/combigrid/OperationPole.hpp"
%include "combigrid/src/sgpp/combigrid/OperationPoleNodalisationBspline.hpp"
%include "combigrid/src/sgpp/combigrid/OperationPoleNodalisationLinear.hpp"
%include "combigrid/src/sgpp/combigrid/OperationUPFullGrid.hpp"
%include "combigrid/src/sgpp/combigrid/OperationUPCombinationGrid.hpp"

namespace std {
  %template(FullGridVector) vector<sgpp::combigrid::FullGrid>;
  %template(LevelVectorVector) vector<sgpp::combigrid::LevelVector>;
  %template(BasisVector) vector<sgpp::base::Basis<sgpp::combigrid::level_t, sgpp::combigrid::index_t>*>;
  %template(OperationPoleVector) vector<sgpp::combigrid::OperationPole*>;
}
