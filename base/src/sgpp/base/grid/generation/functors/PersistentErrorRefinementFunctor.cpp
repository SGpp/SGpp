// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <limits>

#include <sgpp/base/grid/generation/functors/PersistentErrorRefinementFunctor.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    PersistentErrorRefinementFunctor::PersistentErrorRefinementFunctor(
      DataVector* alpha, Grid* grid, size_t refinements_num, float_t threshold) :
      alpha(alpha), refinements_num(refinements_num), threshold(threshold), grid(
        grid), trainDataset(
          NULL), classes(NULL), errors(NULL), accum(NULL) {

    }
    PersistentErrorRefinementFunctor::~PersistentErrorRefinementFunctor() {
    }

    /*
     * current[j] = \sum_{i=0}^{N} (r_i + y_i) * \phi_j(x_i)
     * accum[j] = BETA * accum[j] + current[j]
     * functor value = -alpha_j * accum[j]
     */
    float_t PersistentErrorRefinementFunctor::operator()(GridStorage* storage,
        size_t seq) {

      if (trainDataset == NULL || classes == NULL) {
        throw base::application_exception(
          "Training dataset or classes not set");
      }

      const float_t BETA = 0.1;
      const size_t MIN_SUPPORT = 5;

      size_t numData = trainDataset->getNrows();
      //size_t dim = trainDataset->getNcols();

      DataVector phi_x(numData);

      SGPP::base::DataVector singleAlpha(alpha->getSize());
      singleAlpha.setAll(0.0);
      singleAlpha.set(seq, 1.0);
      SGPP::op_factory::createOperationMultipleEval(*grid, *trainDataset)->mult(
        singleAlpha, phi_x);

      if (phi_x.getNumberNonZero() < MIN_SUPPORT) {
        return start(); // threshold is 0.0
      }

      // Make sure that the error vector is as large as
      // the coefficient vector
      size_t numCoeff = alpha->getSize();

      if (accum == NULL) {
        accum = new SGPP::base::DataVector(numCoeff);
        accum->setAll(start());
      } else if (accum->getSize() != numCoeff) {
        accum->resizeZero(numCoeff);
      }

      // Calculate current

      // Calculate
      // current[j] = \sum_{i=0}^{N} (r_i + y_i) * \phi_j(x_i)

      float_t current_j = 0;
      float_t tmp = 0;

      for (size_t i = 0; i < numData; i++) {
        /* (r_i + y_i) * \phi_j(x_i) */
        tmp = (errors->get(i)) * phi_x[i];
        current_j += tmp * tmp;

      }

      // Accumulation
      /*for (size_t i = 0; i < numCoeff; i++) {
       accum->set(i, accum->get(i) * BETA + current->get(i));
       }*/
      accum->set(seq, accum->get(seq) * (1 - BETA) + BETA * current_j * fabs(alpha->get(seq)));

      float_t func_val = accum->get(seq);

      //std::cout << "Functor value (of " << seq << "): " << func_val << std::endl;
      return func_val;
    }

    float_t PersistentErrorRefinementFunctor::start() {
      return this->threshold;
    }

    size_t PersistentErrorRefinementFunctor::getRefinementsNum() {
      return this->refinements_num;
    }

    float_t PersistentErrorRefinementFunctor::getRefinementThreshold() {
      return this->threshold;
    }

    void PersistentErrorRefinementFunctor::setTrainDataset(
      DataMatrix* trainDataset_) {
      trainDataset = trainDataset_;
    }

    void PersistentErrorRefinementFunctor::setClasses(DataVector* classes_) {
      classes = classes_;
    }
    void PersistentErrorRefinementFunctor::setErrors(DataVector* errors_) {
      errors = errors_;
    }

  }
}
