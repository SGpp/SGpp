/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Benjamin Peherstorfer (pehersto@in.tum.de)

#include <sgpp/datadriven/basis/linear/noboundary/operation/OperationDensityMarginalizeLinear.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

namespace sg {
  namespace datadriven {

    void OperationDensityMarginalizeLinear::doMarginalize(base::DataVector& alpha, base::Grid*& mg, base::DataVector& malpha, unsigned int mdim) {
      /**
       * Note: Because of adaptively refined sparse grids, we cannot simply
       * generate a regular grid. Thus, we need to add point after point
       * to the new grid mg
       */
      //create grid of dimensions d - 1 of the same type
      base::GridStorage* gs = this->grid->getStorage();

      if (gs->dim() < 2)
        throw sg::base::operation_exception("OperationDensityMarginalize is not possible for less than 2 dimensions");

      mg = base::Grid::createLinearGrid(gs->dim() - 1);
      base::GridStorage* mgs = mg->getStorage();

      //run through grid g and add points to mg
      sg::base::GridIndex* gp;
      sg::base::GridIndex mgp(mgs->dim());

      for (unsigned int i = 0; i < gs->size(); i++) {
        gp = gs->get(i);

        for (unsigned int d = 0; d < gs->dim(); d++) {
          //skip direction in which we marginalize
          if (d == mdim)
            continue;
          else if (d < mdim) {
            mgp.set(d, gp->getLevel(d), gp->getIndex(d));
          } else {
            mgp.set(d - 1, gp->getLevel(d), gp->getIndex(d));
          }
        }

        if (!mgs->has_key(&mgp))
          mgs->insert(mgp);
      }

      mgs->recalcLeafProperty();

      /**
       * Compute coefficients for marginalized density
       * Each coefficient has to be weighted with the integral of
       * the basis functions in direction mdim
       */
      malpha.resize(mgs->size());
      malpha.setAll(0.0);
      unsigned int mdimLevel = 0;
      size_t mseqNr;

      for (size_t seqNr = 0; seqNr < gs->size(); seqNr++) {
        gp = gs->get(seqNr);

        for (unsigned int d = 0; d < gs->dim(); d++) {
          if (d == mdim)
            mdimLevel = gp->getLevel(d);
          else if (d < mdim)
            mgp.set(d, gp->getLevel(d), gp->getIndex(d));
          else
            mgp.set(d - 1, gp->getLevel(d), gp->getIndex(d));
        }

        if (!mgs->has_key(&mgp))
          throw sg::base::operation_exception("Key not found! This should not happen! There is something seriously wrong!");

        //get index in alpha vector for current basis function
        mseqNr = mgs->seq(&mgp);
        /**
         * Attention:
         * The integral of one basis functions changes for if another
         * type of basis is used!
         */
        //update corresponding coefficient
        malpha[mseqNr] += alpha[seqNr] * pow(2.0, -static_cast<double>(mdimLevel));
      }
    }
  }
}

