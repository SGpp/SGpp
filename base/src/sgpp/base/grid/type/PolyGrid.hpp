// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POLYGRID_HPP
#define POLYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/globaldef.hpp>


namespace SGPP {
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

  SGPP::base::GridType getType() override;
  const SBasis& getBasis() override;
  void serialize(std::ostream& ostr) override;

  GridGenerator* createGridGenerator() override;

  static Grid* unserialize(std::istream& istr);
  size_t getDegree() const;

 protected:
  /// max. polynom's degree
  size_t degree;
  const SPolyBase* basis_;
};

}  // namespace base
}  // namespace SGPP

#endif /* POLYGRID_HPP */
