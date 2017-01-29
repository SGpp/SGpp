// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeLinear.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

void OperationDensityMarginalizeLinear::doMarginalize(base::DataVector& alpha, base::Grid*& mg,
                                                      base::DataVector& malpha, unsigned int mdim) {
  /**
   * Note: Because of adaptively refined sparse grids, we cannot simply
   * generate a regular grid. Thus, we need to add point after point
   * to the new grid mg
   */
  // create grid of dimensions d - 1 of the same type
  base::GridStorage* gs = &this->grid->getStorage();

  if (gs->getDimension() < 2)
    throw sgpp::base::operation_exception(
        "OperationDensityMarginalize is not possible for less than 2 dimensions");

  mg = base::Grid::createLinearGrid(gs->getDimension() - 1);
  base::GridStorage* mgs = &mg->getStorage();

  // run through grid g and add points to mg
  sgpp::base::GridPoint mgp(mgs->getDimension());

  for (unsigned int i = 0; i < gs->getSize(); i++) {
    sgpp::base::GridPoint& gp = gs->getPoint(i);

    for (unsigned int d = 0; d < gs->getDimension(); d++) {
      // skip direction in which we marginalize
      if (d == mdim) {
        continue;
      } else {
        if (d < mdim) {
          mgp.set(d, gp.getLevel(d), gp.getIndex(d));
        } else {
          mgp.set(d - 1, gp.getLevel(d), gp.getIndex(d));
        }
      }
    }

    if (!mgs->isContaining(mgp)) mgs->insert(mgp);
  }

  mgs->recalcLeafProperty();

  /**
   * Compute coefficients for marginalized density
   * Each coefficient has to be weighted with the integral of
   * the basis functions in direction mdim
   */
  malpha.resize(mgs->getSize());
  malpha.setAll(0.0);
  unsigned int mdimLevel = 0;
  size_t mseqNr;

  for (size_t seqNr = 0; seqNr < gs->getSize(); seqNr++) {
    sgpp::base::GridPoint& gp = gs->getPoint(seqNr);

    for (unsigned int d = 0; d < gs->getDimension(); d++) {
      if (d == mdim)
        mdimLevel = gp.getLevel(d);
      else if (d < mdim)
        mgp.set(d, gp.getLevel(d), gp.getIndex(d));
      else
        mgp.set(d - 1, gp.getLevel(d), gp.getIndex(d));
    }

    if (!mgs->isContaining(mgp))
      throw sgpp::base::operation_exception(
          "Key not found! This should not happen! There is something seriously wrong!");

    // get index in alpha vector for current basis function
    mseqNr = mgs->getSequenceNumber(mgp);
    /**
     * Attention:
     * The integral of one basis functions changes for if another
     * type of basis is used!
     */
    // update corresponding coefficient
    malpha[mseqNr] += alpha[seqNr] * pow(2.0, -static_cast<double>(mdimLevel));
  }
}
}  // namespace datadriven
}  // namespace sgpp
