#include "BasisEval.hpp"

using namespace SGPP;
using namespace SGPP::base;

SGPP::float_t basisEvalDx(
  SBasis& basis, level_t l, index_t i, SGPP::float_t x) {
  SBsplineBase* bsplineBasis =
    dynamic_cast<SBsplineBase*>(&basis);
  SBsplineBoundaryBase* bsplineBoundaryBasis =
    dynamic_cast<SBsplineBoundaryBase*>(&basis);
  SBsplineClenshawCurtisBase* bsplineClenshawCurtisBasis =
    dynamic_cast<SBsplineClenshawCurtisBase*>(&basis);
  SBsplineModifiedBase* bsplineModifiedBasis =
    dynamic_cast<SBsplineModifiedBase*>(&basis);
  SBsplineModifiedClenshawCurtisBase* bsplineModifiedClenshawCurtisBasis =
    dynamic_cast<SBsplineModifiedClenshawCurtisBase*>(&basis);
  SFundamentalSplineBase* fundamentalSplineBasis =
    dynamic_cast<SFundamentalSplineBase*>(&basis);
  SFundamentalSplineModifiedBase* fundamentalSplineModifiedBasis =
    dynamic_cast<SFundamentalSplineModifiedBase*>(&basis);
  SWaveletBase* waveletBasis =
    dynamic_cast<SWaveletBase*>(&basis);
  SWaveletBoundaryBase* waveletBoundaryBasis =
    dynamic_cast<SWaveletBoundaryBase*>(&basis);
  SWaveletModifiedBase* waveletModifiedBasis =
    dynamic_cast<SWaveletModifiedBase*>(&basis);

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

SGPP::float_t basisEvalDxDx(
  SBasis& basis, level_t l, index_t i, SGPP::float_t x) {
  SBsplineBase* bsplineBasis =
    dynamic_cast<SBsplineBase*>(&basis);
  SBsplineBoundaryBase* bsplineBoundaryBasis =
    dynamic_cast<SBsplineBoundaryBase*>(&basis);
  SBsplineClenshawCurtisBase* bsplineClenshawCurtisBasis =
    dynamic_cast<SBsplineClenshawCurtisBase*>(&basis);
  SBsplineModifiedBase* bsplineModifiedBasis =
    dynamic_cast<SBsplineModifiedBase*>(&basis);
  SBsplineModifiedClenshawCurtisBase* bsplineModifiedClenshawCurtisBasis =
    dynamic_cast<SBsplineModifiedClenshawCurtisBase*>(&basis);
  SFundamentalSplineBase* fundamentalSplineBasis =
    dynamic_cast<SFundamentalSplineBase*>(&basis);
  SFundamentalSplineModifiedBase* fundamentalSplineModifiedBasis =
    dynamic_cast<SFundamentalSplineModifiedBase*>(&basis);
  SWaveletBase* waveletBasis =
    dynamic_cast<SWaveletBase*>(&basis);
  SWaveletBoundaryBase* waveletBoundaryBasis =
    dynamic_cast<SWaveletBoundaryBase*>(&basis);
  SWaveletModifiedBase* waveletModifiedBasis =
    dynamic_cast<SWaveletModifiedBase*>(&basis);

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
