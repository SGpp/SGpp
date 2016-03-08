// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "BasisEval.hpp"

double basisEvalDx(sgpp::base::SBasis& basis, sgpp::base::level_t l,
                          sgpp::base::index_t i, double x) {
  sgpp::base::SBsplineBase* bsplineBasis =
      dynamic_cast<sgpp::base::SBsplineBase*>(&basis);
  sgpp::base::SBsplineBoundaryBase* bsplineBoundaryBasis =
      dynamic_cast<sgpp::base::SBsplineBoundaryBase*>(&basis);
  sgpp::base::SBsplineClenshawCurtisBase* bsplineClenshawCurtisBasis =
      dynamic_cast<sgpp::base::SBsplineClenshawCurtisBase*>(&basis);
  sgpp::base::SBsplineModifiedBase* bsplineModifiedBasis =
      dynamic_cast<sgpp::base::SBsplineModifiedBase*>(&basis);
  sgpp::base::SBsplineModifiedClenshawCurtisBase* bsplineModifiedClenshawCurtisBasis =
      dynamic_cast<sgpp::base::SBsplineModifiedClenshawCurtisBase*>(&basis);
  sgpp::base::SFundamentalSplineBase* fundamentalSplineBasis =
      dynamic_cast<sgpp::base::SFundamentalSplineBase*>(&basis);
  sgpp::base::SFundamentalSplineModifiedBase* fundamentalSplineModifiedBasis =
      dynamic_cast<sgpp::base::SFundamentalSplineModifiedBase*>(&basis);
  sgpp::base::SWaveletBase* waveletBasis =
      dynamic_cast<sgpp::base::SWaveletBase*>(&basis);
  sgpp::base::SWaveletBoundaryBase* waveletBoundaryBasis =
      dynamic_cast<sgpp::base::SWaveletBoundaryBase*>(&basis);
  sgpp::base::SWaveletModifiedBase* waveletModifiedBasis =
      dynamic_cast<sgpp::base::SWaveletModifiedBase*>(&basis);

  if (bsplineBasis != nullptr) {
    return bsplineBasis->evalDx(l, i, x);
  } else if (bsplineBoundaryBasis != nullptr) {
    return bsplineBoundaryBasis->evalDx(l, i, x);
  } else if (bsplineClenshawCurtisBasis != nullptr) {
    return bsplineClenshawCurtisBasis->evalDx(l, i, x);
  } else if (bsplineModifiedBasis != nullptr) {
    return bsplineModifiedBasis->evalDx(l, i, x);
  } else if (bsplineModifiedClenshawCurtisBasis != nullptr) {
    return bsplineModifiedClenshawCurtisBasis->evalDx(l, i, x);
  } else if (fundamentalSplineBasis != nullptr) {
    return fundamentalSplineBasis->evalDx(l, i, x);
  } else if (fundamentalSplineModifiedBasis != nullptr) {
    return fundamentalSplineModifiedBasis->evalDx(l, i, x);
  } else if (waveletBasis != nullptr) {
    return waveletBasis->evalDx(l, i, x);
  } else if (waveletBoundaryBasis != nullptr) {
    return waveletBoundaryBasis->evalDx(l, i, x);
  } else if (waveletModifiedBasis != nullptr) {
    return waveletModifiedBasis->evalDx(l, i, x);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Invalid basis."));
    return NAN;
  }
}

double basisEvalDxDx(sgpp::base::SBasis& basis, sgpp::base::level_t l,
                            sgpp::base::index_t i, double x) {
  sgpp::base::SBsplineBase* bsplineBasis =
      dynamic_cast<sgpp::base::SBsplineBase*>(&basis);
  sgpp::base::SBsplineBoundaryBase* bsplineBoundaryBasis =
      dynamic_cast<sgpp::base::SBsplineBoundaryBase*>(&basis);
  sgpp::base::SBsplineClenshawCurtisBase* bsplineClenshawCurtisBasis =
      dynamic_cast<sgpp::base::SBsplineClenshawCurtisBase*>(&basis);
  sgpp::base::SBsplineModifiedBase* bsplineModifiedBasis =
      dynamic_cast<sgpp::base::SBsplineModifiedBase*>(&basis);
  sgpp::base::SBsplineModifiedClenshawCurtisBase* bsplineModifiedClenshawCurtisBasis =
      dynamic_cast<sgpp::base::SBsplineModifiedClenshawCurtisBase*>(&basis);
  sgpp::base::SFundamentalSplineBase* fundamentalSplineBasis =
      dynamic_cast<sgpp::base::SFundamentalSplineBase*>(&basis);
  sgpp::base::SFundamentalSplineModifiedBase* fundamentalSplineModifiedBasis =
      dynamic_cast<sgpp::base::SFundamentalSplineModifiedBase*>(&basis);
  sgpp::base::SWaveletBase* waveletBasis =
      dynamic_cast<sgpp::base::SWaveletBase*>(&basis);
  sgpp::base::SWaveletBoundaryBase* waveletBoundaryBasis =
      dynamic_cast<sgpp::base::SWaveletBoundaryBase*>(&basis);
  sgpp::base::SWaveletModifiedBase* waveletModifiedBasis =
      dynamic_cast<sgpp::base::SWaveletModifiedBase*>(&basis);

  if (bsplineBasis != nullptr) {
    return bsplineBasis->evalDxDx(l, i, x);
  } else if (bsplineBoundaryBasis != nullptr) {
    return bsplineBoundaryBasis->evalDxDx(l, i, x);
  } else if (bsplineClenshawCurtisBasis != nullptr) {
    return bsplineClenshawCurtisBasis->evalDxDx(l, i, x);
  } else if (bsplineModifiedBasis != nullptr) {
    return bsplineModifiedBasis->evalDxDx(l, i, x);
  } else if (bsplineModifiedClenshawCurtisBasis != nullptr) {
    return bsplineModifiedClenshawCurtisBasis->evalDxDx(l, i, x);
  } else if (fundamentalSplineBasis != nullptr) {
    return fundamentalSplineBasis->evalDxDx(l, i, x);
  } else if (fundamentalSplineModifiedBasis != nullptr) {
    return fundamentalSplineModifiedBasis->evalDxDx(l, i, x);
  } else if (waveletBasis != nullptr) {
    return waveletBasis->evalDxDx(l, i, x);
  } else if (waveletBoundaryBasis != nullptr) {
    return waveletBoundaryBasis->evalDxDx(l, i, x);
  } else if (waveletModifiedBasis != nullptr) {
    return waveletModifiedBasis->evalDxDx(l, i, x);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Invalid basis."));
    return NAN;
  }
}
