// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POLYGRID_HPP
#define POLYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * grid with polynomial base functions
 */
class PolyGrid : public Grid {
 protected:
  explicit PolyGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with polynomial base functions
   *
   * @param dim the dimension of the grid
   * @param degree the max. polynom's degree
   */
  PolyGrid(size_t dim, size_t degree);

  /**
   * Destructor
   */
  ~PolyGrid() override;

  sgpp::base::GridType getType() override;
  const SBasis& getBasis() override;
  void serialize(std::ostream& ostr, int version = SERIALIZATION_VERSION) override;

  GridGenerator& getGenerator() override;

  static Grid* unserialize(std::istream& istr);
  size_t getDegree() const;

 protected:
  /// grid generator
  StandardGridGenerator generator;
  /// max. polynom's degree
  size_t degree;
  std::unique_ptr<SPolyBase> basis_;
};

}  // namespace base
}  // namespace sgpp

#endif /* POLYGRID_HPP */
