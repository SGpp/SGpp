// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqrtXPhiPhiDownBBLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

SqrtXPhiPhiDownBBLinear::SqrtXPhiPhiDownBBLinear(sgpp::base::GridStorage* storage)
    : storage(storage), boundingBox(storage->getBoundingBox()) {}

SqrtXPhiPhiDownBBLinear::~SqrtXPhiPhiDownBBLinear() {}

void SqrtXPhiPhiDownBBLinear::operator()(sgpp::base::DataVector& source,
                                         sgpp::base::DataVector& result, grid_iterator& index,
                                         size_t dim) {
  double q = this->boundingBox->getIntervalWidth(dim);
  double t = this->boundingBox->getIntervalOffset(dim);

  bool useBB = false;

  if (q != 1.0 || t != 0.0) {
    useBB = true;
  }

  if (useBB) {
    recBB(source, result, index, dim, 0.0, 0.0, q, t);
  } else {
    rec(source, result, index, dim, 0.0, 0.0);
  }
}

void SqrtXPhiPhiDownBBLinear::rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                  grid_iterator& index, size_t dim, double fl, double fr) {
  size_t seq = index.seq();

  double alpha_value = source[seq];

  sgpp::base::level_t l;
  sgpp::base::index_t i_idx;

  index.get(dim, l, i_idx);

  double i = static_cast<double>(i_idx);
  int l_int = static_cast<int>(l);

  double h = (1.0 / (static_cast<double>(1 << (l_int))));

  double a = fl;
  double b = fr;

  double c = sqrt(i * h + h);
  double d = sqrt(i * h);
  double e = sqrt(i * h - h);

  double hSq = pow(h, 2);
  double hCu = pow(h, 3);
  double iSq = pow(i, 2);
  double iCu = pow(i, 3);

  double diagResultTemp = e * iCu * hCu - 3 * e * iSq * hCu + 3 * e * i * hCu + 7 * d * iSq * hCu -
                           c * iCu * hCu - 3 * c * iSq * hCu - 3 * c * i * hCu - e * hCu - c * hCu;
  diagResultTemp = (-16.0 / 105.0) * (1.0 / (hSq)) * diagResultTemp;

  double downResultTemp =
      2 * c * hCu * a + 5 * c * hCu * b + 4 * d * iCu * hCu * b - 7 * d * iSq * hCu * b -
      4 * d * iCu * hCu * a - 8 * e * hCu * a * i - 6 * e * hCu * b * i - 2 * e * iCu * hCu * b +
      6 * e * iSq * hCu * b + 2 * e * iCu * hCu * a + e * iSq * hCu * a + 5 * e * hCu * a +
      2 * e * hCu * b + 6 * c * hCu * a * i - 2 * c * iCu * hCu * b + c * iSq * hCu * b +
      6 * c * iSq * hCu * a + 2 * c * iCu * hCu * a + 8 * c * i * hCu * b - 7 * d * iSq * hCu * a;

  downResultTemp = (4.0 / 105.0) * (1.0 / (hSq)) * downResultTemp;

  // integration
  result[seq] = downResultTemp + ((diagResultTemp)*alpha_value);

  // dehierarchisation
  double fm = ((fl + fr) / 2.0) + alpha_value;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fl, fm);
    }

    index.stepRight(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fm, fr);
    }

    index.up(dim);
  }
}

void SqrtXPhiPhiDownBBLinear::recBB(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                    grid_iterator& index, size_t dim, double fl, double fr,
                                    double q, double t) {
  size_t seq = index.seq();

  double alpha_value = source[seq];

  sgpp::base::level_t l;
  sgpp::base::index_t i_idx;

  index.get(dim, l, i_idx);

  double i = static_cast<double>(i_idx);
  int l_int = static_cast<int>(l);

  double h = (1.0 / (static_cast<double>(1 << (l_int))));

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

  double a = fl;
  double b = fr;

  double diagResultTemp =
      e * iCu * hCu * qCu - 3 * e * iSq * hCu * qCu + 3 * e * i * hCu * qCu - 3 * e * tSq * q * h +
      3 * e * t * qSq * hSq + 7 * d * tSq * q * h + 7 * d * iSq * hCu * qCu - c * iCu * hCu * qCu -
      3 * c * iSq * hCu * qCu - 3 * c * i * hCu * qCu - 3 * c * tSq * q * h -
      3 * c * t * qSq * hSq + 3 * e * iSq * hSq * qSq * t + 3 * e * i * h * q * tSq -
      6 * e * i * hSq * qSq * t + 14 * d * t * qSq * hSq * i - 3 * c * iSq * hSq * qSq * t -
      3 * c * i * h * q * tSq - 6 * c * i * hSq * qSq * t + e * tCu - c * tCu - e * qCu * hCu -
      c * qCu * hCu;

  diagResultTemp = (-16.0 / 105.0) * (1.0 / (qSq * hSq)) * diagResultTemp;

  double downResultTemp =
      2 * c * qCu * hCu * a + 5 * c * qCu * hCu * b - 14 * d * i * hSq * qSq * a * t +
      6 * c * a * tSq * i * h * q - 6 * c * b * tSq * i * h * q + 6 * c * iSq * hSq * qSq * a * t +
      12 * c * i * hSq * qSq * a * t - 6 * c * iSq * hSq * qSq * b * t +
      2 * c * i * hSq * qSq * b * t + 2 * e * a * tCu - 2 * e * b * tCu + 4 * d * b * tCu -
      4 * d * a * tCu + 4 * d * iCu * hCu * qCu * b - 7 * d * iSq * hCu * qCu * b -
      4 * d * iCu * hCu * qCu * a - 7 * d * tSq * b * q * h - 8 * e * qCu * hCu * a * i -
      8 * e * qSq * hSq * a * t - 6 * e * qCu * hCu * b * i - 6 * e * qSq * hSq * b * t +
      6 * e * b * tSq * q * h + e * a * tSq * q * h - 2 * e * iCu * hCu * qCu * b +
      6 * e * iSq * hCu * qCu * b + 2 * e * iCu * hCu * qCu * a + e * iSq * hCu * qCu * a +
      5 * e * qCu * hCu * a + 2 * e * qCu * hCu * b - 14 * d * i * hSq * qSq * b * t +
      6 * e * iSq * hSq * qSq * a * t + 6 * e * i * h * q * a * tSq +
      2 * e * i * hSq * qSq * a * t - 6 * e * iSq * hSq * qSq * b * t -
      6 * e * i * h * q * b * tSq + 12 * e * i * hSq * qSq * b * t +
      12 * d * iSq * hSq * qSq * b * t + 12 * d * i * h * q * b * tSq -
      12 * d * iSq * hSq * qSq * a * t - 12 * d * i * h * q * a * tSq + 2 * c * a * tCu -
      2 * c * b * tCu + 6 * c * a * tSq * q * h + c * b * tSq * q * h + 6 * c * qCu * hCu * a * i +
      6 * c * qSq * hSq * a * t - 2 * c * iCu * hCu * qCu * b + c * iSq * hCu * qCu * b +
      6 * c * iSq * hCu * qCu * a + 2 * c * iCu * hCu * qCu * a + 8 * c * t * b * qSq * hSq +
      8 * c * i * hCu * qCu * b - 7 * d * iSq * hCu * qCu * a - 7 * d * tSq * a * q * h;

  downResultTemp = (4.0 / 105.0) * (1.0 / (qSq * hSq)) * downResultTemp;

  // integration
  result[seq] = downResultTemp + ((diagResultTemp)*alpha_value);

  // dehierarchisation
  double fm = ((fl + fr) / 2.0) + alpha_value;

  if (!index.hint()) {
    index.leftChild(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      recBB(source, result, index, dim, fl, fm, q, t);
    }

    index.stepRight(dim);

    if (!storage->isValidSequenceNumber(index.seq())) {
      recBB(source, result, index, dim, fm, fr, q, t);
    }

    index.up(dim);
  }
}

}  // namespace finance
}  // namespace sgpp
