// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <list>
#include <type_traits>
#include <vector>

// // cf.
// //
// https://stackoverflow.com/questions/9407367/determine-if-a-type-is-an-stl-container-at-compile-time/31105859#31105859
// namespace is_one_d_stl_container_impl {
// template <typename T>
// struct is_one_d_stl_container : std::false_type {};
// template <typename T, std::size_t N>
// struct is_one_d_stl_container<std::array<T, N>> : std::true_type {};
// template <typename... Args>
// struct is_one_d_stl_container<std::vector<Args...>> : std::true_type {};
// template <typename... Args>
// struct is_one_d_stl_container<std::list<Args...>> : std::true_type {};
// }  // namespace is_one_d_stl_container_impl

namespace std {

// plus operator for std::arrays -- used by std::plus
template <class T, std::size_t N>
std::array<T, N> operator+(const std::array<T, N>& a, const std::array<T, N>& b) {
  std::array<T, N> result;
  for (size_t i = 0; i < N; ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}

// scalar multiplication operator for std::arrays
template <class T, std::size_t N>
std::array<T, N> operator*(const double& a, const std::array<T, N>& b) {
  std::array<T, N> result;
  for (size_t i = 0; i < N; ++i) {
    result[i] = a * b[i];
  }
  return result;
}
// minus operator for std::arrays
template <class T, std::size_t N>
std::array<T, N> operator-(const std::array<T, N>& a, const std::array<T, N>& b) {
  return a + (-1.) * b;
}

// plus operator for std::vectors -- used by std::plus
template <class T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b) {
  std::vector<T> result;
  assert(a.size() == b.size());
  for (size_t i = 0; i < a.size(); ++i) {
    result.emplace_back(a[i] + b[i]);
  }
  return result;
}

// scalar multiplication operator for std::vectors
template <class T>
std::vector<T> operator*(const double& a, const std::vector<T>& b) {
  std::vector<T> result;
  for (size_t i = 0; i < b.size(); ++i) {
    result.emplace_back(a * b[i]);
  }
  return result;
}

// minus operator for std::vectors
template <class T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b) {
  return a + (-1.) * b;
}

}  // namespace std

namespace sgpp {
namespace combigrid {

template <typename T>
void setQuietNan(T& t) {
  std::fill(t.begin(), t.end(), std::numeric_limits<typename T::value_type>::quiet_NaN());
}

template <>
void setQuietNan<double>(double& t) {
  t = std::numeric_limits<double>::quiet_NaN();
}

template <typename T>
T getQuietNan() {
  T t;
  setQuietNan(t);
  return t;
}

template <typename T>
bool isNan(const T& t) {
  return std::any_of(t.begin(), t.end(), [](double d) { return std::isnan(d); });
}

template <>
bool isNan<double>(const double& t) {
  return std::isnan(t);
}

template <typename T>
void setZero(T& t) {
  std::fill(t.begin(), t.end(), 0.);
}

template <>
void setZero<double>(double& t) {
  t = 0.;
}

template <typename T>
T getZero() {
  T t;
  setZero(t);
  return t;
}

}  // namespace combigrid
}  // namespace sgpp