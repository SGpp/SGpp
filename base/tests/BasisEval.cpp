// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include "BasisEval.hpp"

#include <boost/test/unit_test.hpp>


double basisEvalDx(sgpp::base::SBasis& basis, sgpp::base::level_t l, sgpp::base::index_t i,
                   double x) {
  sgpp::base::SBsplineBase* bsplineBasis = dynamic_cast<sgpp::base::SBsplineBase*>(&basis);
  sgpp::base::SBsplineBoundaryBase* bsplineBoundaryBasis =
      dynamic_cast<sgpp::base::SBsplineBoundaryBase*>(&basis);
  sgpp::base::SBsplineClenshawCurtisBase* bsplineClenshawCurtisBasis =
      dynamic_cast<sgpp::base::SBsplineClenshawCurtisBase*>(&basis);
  sgpp::base::SBsplineModifiedBase* bsplineModifiedBasis =
      dynamic_cast<sgpp::base::SBsplineModifiedBase*>(&basis);
  sgpp::base::SBsplineModifiedClenshawCurtisBase* bsplineModifiedClenshawCurtisBasis =
      dynamic_cast<sgpp::base::SBsplineModifiedClenshawCurtisBase*>(&basis);
  sgpp::base::SFundamentalNakSplineBase* fundamentalNakSplineBasis =
      dynamic_cast<sgpp::base::SFundamentalNakSplineBase*>(&basis);
  sgpp::base::SFundamentalSplineBase* fundamentalSplineBasis =
      dynamic_cast<sgpp::base::SFundamentalSplineBase*>(&basis);
  sgpp::base::SFundamentalSplineModifiedBase* fundamentalSplineModifiedBasis =
      dynamic_cast<sgpp::base::SFundamentalSplineModifiedBase*>(&basis);
  sgpp::base::SWeaklyFundamentalNakSplineBase* weaklyFundamentalNakSplineBasis =
      dynamic_cast<sgpp::base::SWeaklyFundamentalNakSplineBase*>(&basis);
  sgpp::base::SWeaklyFundamentalNakSplineModifiedBase* weaklyFundamentalNakSplineModifiedBasis =
      dynamic_cast<sgpp::base::SWeaklyFundamentalNakSplineModifiedBase*>(&basis);
  sgpp::base::SWeaklyFundamentalSplineBase* weaklyFundamentalSplineBasis =
      dynamic_cast<sgpp::base::SWeaklyFundamentalSplineBase*>(&basis);
  sgpp::base::SNakBsplineBase* nakBsplineBasis = dynamic_cast<sgpp::base::SNakBsplineBase*>(&basis);
  sgpp::base::SNakBsplineModifiedBase* nakBsplineModifiedBasis =
      dynamic_cast<sgpp::base::SNakBsplineModifiedBase*>(&basis);
  sgpp::base::SWaveletBase* waveletBasis = dynamic_cast<sgpp::base::SWaveletBase*>(&basis);
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
  } else if (fundamentalNakSplineBasis != nullptr) {
    return fundamentalNakSplineBasis->evalDx(l, i, x);
  } else if (fundamentalSplineBasis != nullptr) {
    return fundamentalSplineBasis->evalDx(l, i, x);
  } else if (fundamentalSplineModifiedBasis != nullptr) {
    return fundamentalSplineModifiedBasis->evalDx(l, i, x);
  } else if (weaklyFundamentalNakSplineBasis != nullptr) {
    sgpp::base::SWeaklyFundamentalNakSplineBaseDeriv1 basisDeriv1(
        weaklyFundamentalNakSplineBasis->getDegree());
    return basisDeriv1.eval(l, i, x);
  } else if (weaklyFundamentalNakSplineModifiedBasis != nullptr) {
    sgpp::base::SWeaklyFundamentalNakSplineModifiedBaseDeriv1 basisDeriv1(
        weaklyFundamentalNakSplineModifiedBasis->getDegree());
    return basisDeriv1.eval(l, i, x);
  } else if (weaklyFundamentalSplineBasis != nullptr) {
    sgpp::base::SWeaklyFundamentalSplineBaseDeriv1 basisDeriv1(
        weaklyFundamentalSplineBasis->getDegree());
    return basisDeriv1.eval(l, i, x);
  } else if (nakBsplineBasis != nullptr) {
    sgpp::base::SNakBsplineBaseDeriv1 basisDeriv1(nakBsplineBasis->getDegree());
    return basisDeriv1.eval(l, i, x);
  } else if (nakBsplineModifiedBasis != nullptr) {
    sgpp::base::SNakBsplineModifiedBaseDeriv1 basisDeriv1(nakBsplineModifiedBasis->getDegree());
    return basisDeriv1.eval(l, i, x);
    // return nakBsplineModifiedBasis->evalDx(l, i, x);
  } else if (waveletBasis != nullptr) {
    return waveletBasis->evalDx(l, i, x);
  } else if (waveletBoundaryBasis != nullptr) {
    return waveletBoundaryBasis->evalDx(l, i, x);
  } else if (waveletModifiedBasis != nullptr) {
    return waveletModifiedBasis->evalDx(l, i, x);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Invalid basis."));
  }
}

double basisEvalDxDx(sgpp::base::SBasis& basis, sgpp::base::level_t l, sgpp::base::index_t i,
                     double x) {
  sgpp::base::SBsplineBase* bsplineBasis = dynamic_cast<sgpp::base::SBsplineBase*>(&basis);
  sgpp::base::SBsplineBoundaryBase* bsplineBoundaryBasis =
      dynamic_cast<sgpp::base::SBsplineBoundaryBase*>(&basis);
  sgpp::base::SBsplineClenshawCurtisBase* bsplineClenshawCurtisBasis =
      dynamic_cast<sgpp::base::SBsplineClenshawCurtisBase*>(&basis);
  sgpp::base::SBsplineModifiedBase* bsplineModifiedBasis =
      dynamic_cast<sgpp::base::SBsplineModifiedBase*>(&basis);
  sgpp::base::SBsplineModifiedClenshawCurtisBase* bsplineModifiedClenshawCurtisBasis =
      dynamic_cast<sgpp::base::SBsplineModifiedClenshawCurtisBase*>(&basis);
  sgpp::base::SFundamentalNakSplineBase* fundamentalNakSplineBasis =
      dynamic_cast<sgpp::base::SFundamentalNakSplineBase*>(&basis);
  sgpp::base::SFundamentalSplineBase* fundamentalSplineBasis =
      dynamic_cast<sgpp::base::SFundamentalSplineBase*>(&basis);
  sgpp::base::SFundamentalSplineModifiedBase* fundamentalSplineModifiedBasis =
      dynamic_cast<sgpp::base::SFundamentalSplineModifiedBase*>(&basis);
  sgpp::base::SWeaklyFundamentalNakSplineBase* weaklyFundamentalNakSplineBasis =
      dynamic_cast<sgpp::base::SWeaklyFundamentalNakSplineBase*>(&basis);
  sgpp::base::SWeaklyFundamentalNakSplineModifiedBase* weaklyFundamentalNakSplineModifiedBasis =
      dynamic_cast<sgpp::base::SWeaklyFundamentalNakSplineModifiedBase*>(&basis);
  sgpp::base::SWeaklyFundamentalSplineBase* weaklyFundamentalSplineBasis =
      dynamic_cast<sgpp::base::SWeaklyFundamentalSplineBase*>(&basis);
  sgpp::base::SNakBsplineBase* nakBsplineBasis = dynamic_cast<sgpp::base::SNakBsplineBase*>(&basis);
  sgpp::base::SNakBsplineModifiedBase* nakBsplineModifiedBasis =
      dynamic_cast<sgpp::base::SNakBsplineModifiedBase*>(&basis);
  sgpp::base::SWaveletBase* waveletBasis = dynamic_cast<sgpp::base::SWaveletBase*>(&basis);
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
  } else if (fundamentalNakSplineBasis != nullptr) {
    return fundamentalNakSplineBasis->evalDxDx(l, i, x);
  } else if (fundamentalSplineBasis != nullptr) {
    return fundamentalSplineBasis->evalDxDx(l, i, x);
  } else if (fundamentalSplineModifiedBasis != nullptr) {
    return fundamentalSplineModifiedBasis->evalDxDx(l, i, x);
  } else if (weaklyFundamentalNakSplineBasis != nullptr) {
    sgpp::base::SWeaklyFundamentalNakSplineBaseDeriv2 basisDeriv2(
        weaklyFundamentalNakSplineBasis->getDegree());
    return basisDeriv2.eval(l, i, x);
  } else if (weaklyFundamentalNakSplineModifiedBasis != nullptr) {
    sgpp::base::SWeaklyFundamentalNakSplineModifiedBaseDeriv2 basisDeriv2(
        weaklyFundamentalNakSplineModifiedBasis->getDegree());
    return basisDeriv2.eval(l, i, x);
  } else if (weaklyFundamentalSplineBasis != nullptr) {
    sgpp::base::SWeaklyFundamentalSplineBaseDeriv2 basisDeriv2(
        weaklyFundamentalSplineBasis->getDegree());
    return basisDeriv2.eval(l, i, x);
  } else if (nakBsplineBasis != nullptr) {
    sgpp::base::SNakBsplineBaseDeriv2 basisDeriv2(nakBsplineBasis->getDegree());
    return basisDeriv2.eval(l, i, x);
  } else if (nakBsplineModifiedBasis != nullptr) {
    sgpp::base::SNakBsplineModifiedBaseDeriv2 basisDeriv2(nakBsplineModifiedBasis->getDegree());
    return basisDeriv2.eval(l, i, x);
  } else if (waveletBasis != nullptr) {
    return waveletBasis->evalDxDx(l, i, x);
  } else if (waveletBoundaryBasis != nullptr) {
    return waveletBoundaryBasis->evalDxDx(l, i, x);
  } else if (waveletModifiedBasis != nullptr) {
    return waveletModifiedBasis->evalDxDx(l, i, x);
  } else {
    BOOST_THROW_EXCEPTION(std::runtime_error("Invalid basis."));
  }
}
