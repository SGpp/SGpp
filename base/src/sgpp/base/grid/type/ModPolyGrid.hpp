// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODPOLYGRID_HPP
#define MODPOLYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * grid with modified polynomial base functions
 */
class ModPolyGrid : public Grid {
 protected:
  explicit ModPolyGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with modified polynomial base functions
   *
   * @param dim the dimension of the grid
   * @param degree the max. polynom's degree
   */
  ModPolyGrid(size_t dim, size_t degree);

  /**
   * Destructor
   */
  ~ModPolyGrid() override;

  sgpp::base::GridType getType() override;
  void serialize(std::ostream& ostr, int version = SERIALIZATION_VERSION) override;

  SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);
  virtual size_t getDegree() const;

 protected:
  /// grid generator
  StandardGridGenerator generator;
  /// max. polynom's degree
  size_t degree;
  std::unique_ptr<SPolyModifiedBase> basis_;
};

}  // namespace base
}  // namespace sgpp

#endif /* MODPOLYGRID_HPP */
