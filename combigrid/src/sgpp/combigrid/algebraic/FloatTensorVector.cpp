// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>

#include <cmath>

namespace sgpp {
namespace combigrid {

void FloatTensorVector::ensureDim(size_t dim) {
  if (d < dim) {
    d = dim;
    if (d == 0 && dim == 1) {
      return;
    }
    auto newStorage = std::make_shared<TreeStorage<FloatScalarVector>>(dim);
    auto it = values->getStoredDataIterator();
    for (; it->isValid(); it->moveToNext()) {
      MultiIndex index = it->getMultiIndex();
      for (size_t i = index.size(); i < dim; ++i) {
        index.push_back(0);
      }
      newStorage->set(index, it->value());
    }
    values = newStorage;
  }
}

FloatTensorVector::FloatTensorVector(const std::shared_ptr<TreeStorage<FloatScalarVector>>& values)
    : d(values->getNumDimensions()), values(values) {}

FloatTensorVector::FloatTensorVector(size_t d) : d(d), values(nullptr) {
  if (d > 0) {
    values = std::make_shared<TreeStorage<FloatScalarVector>>(d);
  } else {
    values = std::make_shared<TreeStorage<FloatScalarVector>>(1);
  }
}

FloatTensorVector::FloatTensorVector(FloatScalarVector scalar)
    : d(0), values(std::make_shared<TreeStorage<FloatScalarVector>>(1)) {
  values->set(MultiIndex{0}, scalar);
}

FloatTensorVector::FloatTensorVector(const FloatTensorVector& other) : d(other.d) {
  values = std::make_shared<TreeStorage<FloatScalarVector>>(other.values->getNumDimensions());
  auto it = other.values->getStoredDataIterator();
  for (; it->isValid(); it->moveToNext()) {
    values->set(it->getMultiIndex(), it->value());
  }
}

FloatTensorVector& FloatTensorVector::operator=(const FloatTensorVector& other) {
  if (&other != this) {
    d = other.d;
    values = std::make_shared<TreeStorage<FloatScalarVector>>(other.values->getNumDimensions());
    auto it = other.values->getStoredDataIterator();
    for (; it->isValid(); it->moveToNext()) {
      values->set(it->getMultiIndex(), it->value());
    }
  }

  return *this;
}

FloatScalarVector FloatTensorVector::get(MultiIndex i) {
  if (d == 0) {
    return values->get(MultiIndex{0});
  }
  return values->get(i);
}

FloatScalarVector& FloatTensorVector::at(MultiIndex i) {
  if (d == 0) {
    return values->get(MultiIndex{0});
  }
  return values->get(i);
}

FloatScalarVector& FloatTensorVector::operator[](MultiIndex i) {
  if (d == 0) {
    return values->get(MultiIndex{0});
  }
  return values->get(i);
}

FloatScalarVector FloatTensorVector::operator[](size_t i) {
  std::cerr << "FloatTensorVector: operator[] cannot be used with a size_t parameter. Use a "
               "MultiIndex instead. Returning zero."
            << std::endl;
  return FloatScalarVector::zero();
}

void FloatTensorVector::add(const FloatTensorVector& other) {
  ensureDim(other.d);
  auto it = other.values->getStoredDataIterator();
  for (; it->isValid(); it->moveToNext()) {
    MultiIndex index = it->getMultiIndex();
    for (size_t i = index.size(); i < values->getNumDimensions(); ++i) {
      index.push_back(0);
    }
    values->get(index).add(it->value());
  }
}

void FloatTensorVector::sub(const FloatTensorVector& other) {
  ensureDim(other.d);
  auto it = other.values->getStoredDataIterator();
  for (; it->isValid(); it->moveToNext()) {
    MultiIndex index = it->getMultiIndex();
    for (size_t i = index.size(); i < values->getNumDimensions(); ++i) {
      index.push_back(0);
    }
    values->get(index).sub(it->value());
  }
}

void FloatTensorVector::componentwiseMult(const FloatTensorVector& other) {
  // compute tensor product
  if (d == 0 && other.d == 0) {
    values->get(MultiIndex{0}).componentwiseMult(other.values->get(MultiIndex{0}));
  } else if (other.d == 0) {
    scalarMult(other.values->get(MultiIndex{0}).value());
  } else if (d == 0) {
    double scalar = values->get(MultiIndex{0}).value();
    *this = other;
    scalarMult(scalar);
  } else {
    size_t dnew = d + other.d;
    auto newValues = std::make_shared<TreeStorage<FloatScalarVector>>(dnew);
    auto it1 = values->getStoredDataIterator();
    for (; it1->isValid(); it1->moveToNext()) {
      auto it2 = other.values->getStoredDataIterator();
      for (; it2->isValid(); it2->moveToNext()) {
        MultiIndex index = it1->getMultiIndex();
        MultiIndex index2 = it2->getMultiIndex();
        index.insert(index.end(), index2.begin(), index2.end());
        FloatScalarVector scalar = it1->value();
        scalar.componentwiseMult(it2->value());
        newValues->set(index, scalar);
      }
    }
    d = dnew;
    values = newValues;
  }
}

void FloatTensorVector::scalarMult(const double& factor) {
  auto it = values->getStoredDataIterator();
  for (; it->isValid(); it->moveToNext()) {
    it->value().scalarMult(factor);
  }
}

double FloatTensorVector::norm() const {
  double sum = 0.0;
  for (auto it = values->getStoredDataIterator(); it->isValid(); it->moveToNext()) {
    double coeff = it->value().value();
    sum += coeff * coeff;
  }
  return std::sqrt(sum);
}

} /* namespace combigrid */
} /* namespace sgpp */
