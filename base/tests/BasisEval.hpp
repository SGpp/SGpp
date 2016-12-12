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
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeNotAKnotSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeNotAKnotSplineBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeNotAKnotSplineBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeNotAKnotSplineBasisDeriv3.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeNotAKnotSplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeNotAKnotSplineModifiedBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeNotAKnotSplineModifiedBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeNotAKnotSplineModifiedBasisDeriv3.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeSplineBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeSplineBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/LagrangeSplineBasisDeriv3.hpp>
#include <sgpp/base/operation/hash/common/basis/NaturalBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NotAKnotBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NotAKnotBsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

double basisEvalDx(sgpp::base::SBasis& basis, sgpp::base::level_t l,
                          sgpp::base::index_t i, double x);
double basisEvalDxDx(sgpp::base::SBasis& basis, sgpp::base::level_t l,
                            sgpp::base::index_t i, double x);

#endif /* BASIS_EVAL_HPP */
