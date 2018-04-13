// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONCONVERTPREWAVELET_HPP
#define OPERATIONCONVERTPREWAVELET_HPP

#include <sgpp/base/operation/hash/OperationConvert.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 *
 */
class OperationConvertPrewavelet : public OperationConvert {
 public:
  /**
   * Constructor of OperationHierarchisationPrewavelet
   *
     * An adaptive grid with prewavelet ansatz functions requires for operations
   * using the up-down algorithm shadow points. These shadow points a needed just
   * for data transport, thus they do not have an influence on the final function.
   * Please refer to sgpp::pde::UpDownOneOpDimWithShadow for more information.
   *
     * @param storage Pointer to the grid's gridstorage obejct
   * @param shadowstorage shadow points (see detailed description)
   */
  OperationConvertPrewavelet(GridStorage& storage, GridStorage& shadowstorage) :
    storage(storage), shadowstorage(shadowstorage) {
  }

  /**
   * Destructor
   */
  ~OperationConvertPrewavelet() override {
  }

  void doConvertToLinear(DataVector& alpha) override;
  void doConvertFromLinear(DataVector& alpha) override;

 protected:
  /// reference to the grid's GridStorage object
  GridStorage& storage;
  GridStorage& shadowstorage;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONCONVERTPREWAVELET_HPP */
