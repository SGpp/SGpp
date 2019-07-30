// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BASIS_EVAL_HPP
#define BASIS_EVAL_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalNakSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineModifiedBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineModifiedBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/NaturalBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include "../src/sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBoundaryBasis.hpp"

double basisEvalDx(sgpp::base::SBasis& basis, sgpp::base::level_t l, sgpp::base::index_t i,
                   double x);
double basisEvalDxDx(sgpp::base::SBasis& basis, sgpp::base::level_t l, sgpp::base::index_t i,
                     double x);

#endif /* BASIS_EVAL_HPP */
