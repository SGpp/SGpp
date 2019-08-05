// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/sle/system/SLE.hpp>
#include <sgpp/globaldef.hpp>

#include <memory>

namespace sgpp {
namespace base {

/**
 * Abstract class for "cloneable" linear systems.
 * This class is needed in the case that matrix entry lookups are not
 * possible concurrently
 * (e.g. for hierarchization systems with Clenshaw-Curtis grids).
 */
class CloneableSLE : public SLE {
 public:
  /**
   * Constructor.
   */
  CloneableSLE() : SLE() {}

  /**
   * Destructor.
   */
  ~CloneableSLE() override {}

  /**
   * Pure virtual method for cloning the linear system.
   * It should generate a pointer to the cloned object and it's used for
   * parallel computations
   * (e.g. the getMatrixEntry() method might not be thread-safe).
   *
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<CloneableSLE>& clone) const = 0;

  /**
   * @return whether this system derives from CloneableSLE or not (true)
   */
  bool isCloneable() const override { return true; }
};
}  // namespace base
}  // namespace sgpp
