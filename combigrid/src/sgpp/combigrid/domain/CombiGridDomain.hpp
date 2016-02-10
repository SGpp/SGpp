/* ****************************************************************************
 * Copyright (C) 2011 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Janos Benk (benk@in.tum.de)
#ifndef COMBIGRIDDOMAIN_HPP_
#define COMBIGRIDDOMAIN_HPP_

#include <sgpp/combigrid/domain/AbstractStretchingMaker.hpp>
#include <sgpp/combigrid/domain/CombiDomain1D.hpp>
#include "../utils/combigrid_utils.hpp"

namespace combigrid {

/** grid domain*/
class GridDomain {
 public:
  /** Constructor for homogeneously stretched griddomain*/
  GridDomain(int dim, const std::vector<int>& levels,
             const std::vector<double>& min, const std::vector<double>& max,
             AbstractStretchingMaker& stretchingMaker);
  /** constructor for heterogeneously stretched griddomain*/
  GridDomain(int dim, const std::vector<int>& levels,
             const std::vector<double>& min, const std::vector<double>& max,
             std::vector<AbstractStretchingMaker*> stretchingMakers);

  /* copy constructor*/
  GridDomain(const GridDomain& domain);

  /**
   * Destructor!
   */
  virtual ~GridDomain() { ; }

  int getDim() const { return (int)_axisDomains.size(); }
  /** transform from real coordinate into unit coordinates
   * @param coords [IN/OUT]
   * @param levels_in [IN] the required levels
   * @param boundaryFlag [IN] for each dimensions if there are boundary points*/
  void transformRealToUnit(std::vector<double>& coords,
                           const std::vector<int>& levels_in,
                           const std::vector<bool>& boundaryFlag) const;

  /** return 1D axis, can be used for back transformation for each dimension
   * @param d [IN] the dimension */
  const Domain1D& get1DDomain(int d) const { return _axisDomains[d]; }
  /*Tell, to whoever is asking, which stretching type does this domain define...
   *
   **/

  Stretching getStretchingType() const { return _stretching_type; }

  inline std::vector<double> getMin() const { return _min; }

  inline std::vector<double> getMax() const { return _max; }

  void printDomain();

 private:
  /** dimension of the domain */
  int dim_;

  /** Domain needs to know about the underlying stretching of its grid points */
  Stretching _stretching_type;

  /** array to store the maping for each axis */
  std::vector<Domain1D> _axisDomains;

  std::vector<double> _min;
  std::vector<double> _max;
};
}

#endif /* COMBIGRIDDOMAIN_HPP_ */
