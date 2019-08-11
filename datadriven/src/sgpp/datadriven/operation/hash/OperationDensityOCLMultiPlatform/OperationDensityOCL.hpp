// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/KernelMult.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/KernelB.hpp>

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// Base class for density multiplication operation
class OperationDensity: public base::OperationMatrix {
 public:
  OperationDensity()  {
  }
  /// Execute one matrix-vector multiplication with the density matrix
  virtual void mult(base::DataVector& alpha, base::DataVector& result) = 0;
  /// Use before calling partial_mult directly
  virtual void initialize_alpha(double *alpha) = 0;
  /// Execute a partial (startindex to startindex+chunksize) multiplication with the density matrix
  virtual void start_partial_mult(int start_id, int chunksize) = 0;
  virtual void finish_partial_mult(double *result, int start_id, int chunksize) = 0;
  /// Generates the right hand side vector for the density equation
  virtual void generateb(base::DataMatrix &dataset, sgpp::base::DataVector &b,
                         size_t start_id = 0,  size_t chunksize = 0) = 0;
  virtual void initialize_dataset(base::DataMatrix &dataset) = 0;
  virtual void start_rhs_generation(size_t start_id,  size_t chunksize) = 0;
  virtual void finalize_rhs_generation(sgpp::base::DataVector &b,
                         size_t start_id,  size_t chunksize) = 0;
  /// Generate the default parameters in die json configuration
  static void load_default_parameters(base::OCLOperationConfiguration *parameters) {
  if (parameters->contains("INTERNAL_PRECISION") == false) {
    std::cout << "Warning! No internal precision setting detected."
              << " Using double precision from now on!" << std::endl;
    parameters->addIDAttr("INTERNAL_PRECISION", "double");
  }

  if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    DensityOCLMultiPlatform::KernelDensityMult<float>::augmentDefaultParameters(*parameters);
    DensityOCLMultiPlatform::KernelDensityB<float>::augmentDefaultParameters(*parameters);
  } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    DensityOCLMultiPlatform::KernelDensityMult<double>::augmentDefaultParameters(*parameters);
    DensityOCLMultiPlatform::KernelDensityB<double>::augmentDefaultParameters(*parameters);
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"OperationDensityOCLMultiPlatform\": "
                << " invalid value for parameter \"INTERNAL_PRECISION\"";
    throw base::operation_exception(errorString.str().c_str());
  }
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
