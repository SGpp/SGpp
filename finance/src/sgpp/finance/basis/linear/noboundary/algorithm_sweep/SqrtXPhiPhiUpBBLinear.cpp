// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqrtXPhiPhiUpBBLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

SqrtXPhiPhiUpBBLinear::SqrtXPhiPhiUpBBLinear(sgpp::base::GridStorage* storage)
    : storage(storage), boundingBox(storage->getBoundingBox()) {}

SqrtXPhiPhiUpBBLinear::~SqrtXPhiPhiUpBBLinear() {}

void SqrtXPhiPhiUpBBLinear::operator()(sgpp::base::DataVector& source,
                                       sgpp::base::DataVector& result, grid_iterator& index,
                                       size_t dim) {
  double q = boundingBox->getIntervalWidth(dim);
  double t = boundingBox->getIntervalOffset(dim);

  bool useBB = false;

  if (q != 1.0 || t != 0.0) {
    useBB = true;
  }

  // get boundary values
  double fl = 0.0;
  double fr = 0.0;

  if (useBB) {
    recBB(source, result, index, dim, fl, fr, q, t);
  } else {
    rec(source, result, index, dim, fl, fr);
  }
}

void SqrtXPhiPhiUpBBLinear::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                grid_iterator& index, size_t dim, double& fl, double& fr) {
  size_t seq = index.seq();

  fl = fr = 0.0;
  double fml = 0.0;
  double fmr = 0.0;

  sgpp::base::GridStorage::index_type::level_type current_level;
  sgpp::base::GridStorage::index_type::index_type current_index;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fl, fml);
    }

    index.stepRight(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fmr, fr);
    }

    index.up(dim);
  }

  index.get(dim, current_level, current_index);

  double fm = fml + fmr;

  double alpha_value = source[seq];

  double i = static_cast<double>(current_index);
  double h = (1.0 / static_cast<double>(1 << current_level));

  double c = sqrt(i * h + h);
  double d = sqrt(i * h);
  double e = sqrt(i * h - h);

  double hSq = pow(h, 2);
  double hCu = pow(h, 3);
  double iSq = pow(i, 2);
  double iCu = pow(i, 3);

  double flTemp = 2 * c * hCu - 4 * d * iCu * hCu - 8 * e * hCu * i + 2 * e * iCu * hCu +
                   e * iSq * hCu + 5 * e * hCu + 6 * c * hCu * i + 6 * c * iSq * hCu +
                   2 * c * iCu * hCu - 7 * d * iSq * hCu;

  flTemp = (4.0 / 105.0) * (1 / (hSq)) * flTemp;

  double frTemp = 5 * c * hCu + 4 * d * iCu * hCu - 7 * d * iSq * hCu - 6 * e * hCu * i -
                   2 * e * iCu * hCu + 6 * e * iSq * hCu + 2 * e * hCu - 2 * c * iCu * hCu +
                   c * iSq * hCu + 8 * c * i * hCu;

  frTemp = (4.0 / 105.0) * (1 / (hSq)) * frTemp;

  fl = ((fm / 2.0) + (alpha_value * (flTemp))) + fl;
  fr = ((fm / 2.0) + (alpha_value * (frTemp))) + fr;

  // transposed operations:
  result[seq] = fm;
}

void SqrtXPhiPhiUpBBLinear::recBB(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                  grid_iterator& index, size_t dim, double& fl, double& fr,
                                  double q, double t) {
  size_t seq = index.seq();

  fl = fr = 0.0;
  double fml = 0.0;
  double fmr = 0.0;

  sgpp::base::GridStorage::index_type::level_type current_level;
  sgpp::base::GridStorage::index_type::index_type current_index;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      recBB(source, result, index, dim, fl, fml, q, t);
    }

    index.stepRight(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      recBB(source, result, index, dim, fmr, fr, q, t);
    }

    index.up(dim);
  }

  index.get(dim, current_level, current_index);

  double fm = fml + fmr;

  double alpha_value = source[seq];
  double h = (1.0 / static_cast<double>(1 << current_level));

  // transposed operations:
  result[seq] = fm;

  double i = static_cast<double>(current_index);

  double c = sqrt(i * h * q + t + q * h);
  double d = sqrt(i * h * q + t);
  double e = sqrt(i * h * q + t - q * h);

  double qSq = pow(q, 2);
  double qCu = pow(q, 3);
  double hSq = pow(h, 2);
  double hCu = pow(h, 3);
  double tSq = pow(t, 2);
  double tCu = pow(t, 3);
  double iSq = pow(i, 2);
  double iCu = pow(i, 3);

  double flTemp =
      2 * c * qCu * hCu - 14 * d * i * hSq * qSq * t + 6 * c * tSq * i * h * q +
      6 * c * iSq * hSq * qSq * t + 12 * c * i * hSq * qSq * t + 2 * e * tCu - 4 * d * tCu -
      4 * d * iCu * hCu * qCu - 8 * e * qCu * hCu * i - 8 * e * qSq * hSq * t + e * tSq * q * h +
      2 * e * iCu * hCu * qCu + e * iSq * hCu * qCu + 5 * e * qCu * hCu +
      6 * e * iSq * hSq * qSq * t + 6 * e * i * h * q * tSq + 2 * e * i * hSq * qSq * t -
      12 * d * iSq * hSq * qSq * t - 12 * d * i * h * q * tSq + 2 * c * tCu + 6 * c * tSq * q * h +
      6 * c * qCu * hCu * i + 6 * c * qSq * hSq * t + 6 * c * iSq * hCu * qCu +
      2 * c * iCu * hCu * qCu - 7 * d * iSq * hCu * qCu - 7 * d * tSq * q * h;

  flTemp = (4.0 / 105.0) * (1 / (qSq * hSq)) * flTemp;

  double frTemp =
      5 * c * qCu * hCu - 6 * c * tSq * i * h * q - 6 * c * iSq * hSq * qSq * t +
      2 * c * i * hSq * qSq * t - 2 * e * tCu + 4 * d * tCu + 4 * d * iCu * hCu * qCu -
      7 * d * iSq * hCu * qCu - 7 * d * tSq * q * h - 6 * e * qCu * hCu * i -
      6 * e * qSq * hSq * t + 6 * e * tSq * q * h - 2 * e * iCu * hCu * qCu +
      6 * e * iSq * hCu * qCu + 2 * e * qCu * hCu - 14 * d * i * hSq * qSq * t -
      6 * e * iSq * hSq * qSq * t - 6 * e * i * h * q * tSq + 12 * e * i * hSq * qSq * t +
      12 * d * iSq * hSq * qSq * t + 12 * d * i * h * q * tSq - 2 * c * tCu + c * tSq * q * h -
      2 * c * iCu * hCu * qCu + c * iSq * hCu * qCu + 8 * c * t * qSq * hSq + 8 * c * i * hCu * qCu;

  frTemp = (4.0 / 105.0) * (1 / (qSq * hSq)) * frTemp;

  fl = ((fm / 2.0) + (alpha_value * (flTemp))) + fl;
  fr = ((fm / 2.0) + (alpha_value * (frTemp))) + fr;
}

}  // namespace finance
}  // namespace sgpp
