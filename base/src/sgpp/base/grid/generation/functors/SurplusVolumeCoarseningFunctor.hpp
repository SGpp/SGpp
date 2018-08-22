// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SURPLUSVOLUMECOARSENINGFUNCTOR_HPP
#define SURPLUSVOLUMECOARSENINGFUNCTOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * A coarsening functor, removing points according to the minimal absolute values in a DataVector provided,
 * weighted with the corresponding basis function's surplus, i.e., with @f$2^{-|\vec{l}|_1} = 2^{\sum_{k=1}^d l_d}@f$.
 */
class SurplusVolumeCoarseningFunctor : public CoarseningFunctor {
 public:
  /**
   * Constructor.
   *
   * @param alpha DataVector that is basis for coarsening decisions. The i-th entry corresponds to the i-th grid point.
   * @param removements_num Number of grid points which should be removed (if possible - there could be less removable grid points)
   * @param threshold The absolute value of the entries have to be less or equal than the threshold to be considered for coarsening
   */
  SurplusVolumeCoarseningFunctor(DataVector& alpha, size_t removements_num = 1,
                                 double threshold = 0.0);

  /**
   * Destructor
   */
  ~SurplusVolumeCoarseningFunctor() override;

  double operator()(GridStorage& storage, size_t seq) override;

  double start() const override;

  size_t getRemovementsNum() const override;

  double getCoarseningThreshold() const override;

 protected:
  /// pointer to the vector that stores the alpha values
  DataVector& alpha;

  /// number of grid points to remove
  size_t removements_num;

  /**
   * threshold, only the points with greater to equal absolute values of the
   * refinement criterion (e.g. alpha or error) will be refined
   */
  double threshold;
};

}  // namespace base
}  // namespace sgpp

#endif /* SURPLUSVOLUMECOARSENINGFUNCTOR_HPP */
