// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractBasisCoefficientsStorage.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/SLECoefficientsStorage.hpp>

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {} /* namespace combigrid */
} /* namespace sgpp */
