// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATAVECTORDEFINITION_HPP
#define DATAVECTORDEFINITION_HPP

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {
/**
 * This struct is needed for exporting a DataVector
 * to another address space, so it contains all
 * information that is needed to reconstruct a
 * DataVector object
 *
 * The space required by a DataVector object is:
 * (size+unused)*sizeof(double)
 */
struct DataVectorDefinition {
  /// Array to store the data
  double* data;
  /// Number of elements of the data vector
  size_t size;
  /// Number of additional rows for which memory has already been reserved
  size_t unused;
  /**
   * Number of elements by which the reserved memory is increased,
   * if adding an element would exceed the storage reserved so far.
   */
  size_t inc_elems;
};

}  // namespace base
}  // namespace sgpp

#endif /* DATAVECTORDEFINITION_HPP */
