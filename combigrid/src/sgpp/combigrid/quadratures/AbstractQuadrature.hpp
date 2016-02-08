/*
 * AbstractQuadrature.hpp
 *
 *  Created on: 21 Jul 2014
 *      Author: kenny
 */

#ifndef ABSTRACTQUADRATURE_HPP_
#define ABSTRACTQUADRATURE_HPP_

#include <sgpp/combigrid/domain/CombiGridDomain.hpp>
#include <sgpp/combigrid/combigrid/CombiGrid.hpp>
#include <sgpp/combigrid/utils/combigrid_ultils.hpp>

namespace combigrid {

template<typename _Tp>

class AbstractQuadratureRule {
 public:

  /**
   * Pure virtual, interface function which is inherited from all quadrature rules
   * residing in the "quadrature" folder. Multi-dimensional quadratures on combi-grid
   * are implemented employing the "smolyaks'" construction, where the numerical integral
   * on the combigrid is simply a linear combination of the numerical integral values of the
   * underlying fullgrids, where the coefficients are namely the stadnard combination technique
   * coefficients.
   *
   * @param grids - an implementation of the tamplate Combigrid containing the data we wish to intergrate on
   * @param f - an optional function pointer can be handed over if the user wishes to apply the numerical quadrature
   * on a specific function "f" . if no such pointer is supplied, or equivalently a NULL value is given as a 3rd argument
   * of the integrate funciton, the quadrature rule is applied to whatever information is stored on the points of the
   * underlying full grids...
   */
  virtual _Tp integrate(CombiGrid<_Tp>* grids,
                        _Tp (*f)(std::vector<double>) = NULL) = 0;

  virtual ~AbstractQuadratureRule() {
    ;
  }

 protected:
  void getGridValues(FullGrid<_Tp>* grid, bool badstretching,
                     AbstractStretchingMaker* stretching, std::vector<_Tp>& f_values,
                     _Tp (*f)(std::vector<double>)) {

    int dim = grid->getDimension();
    std::vector<Domain1D> domains;
    GridDomain* domain = grid->getDomain();

    if (domain != NULL) {
      for (int d = 0; d < dim; d++) {
        domains.push_back(domain->get1DDomain(d));
      }

    }

    unsigned int num_elem =
      grid->getNrElements(); // get the total number of grid points

    // precomupte the function values depending on whether or not the user wants to integrate an external function
    // by the function pointer or he/she wants to sample the values stored on the grid...
    f_values.clear();
    f_values.resize(num_elem, 0.0);
    std::vector<int> indices(dim, 0);
    std::vector<double> coords(dim, 0);

    if (f == NULL) {
      if (badstretching) {
        std::vector<std::vector<double> > stretchings(dim,
            std::vector<double>());
        std::vector<std::vector<double> > jacobs(dim,
            std::vector<double>());

        for (int d = 0; d < dim; d++)
          stretching->get1DStretching(grid->getLevels()[d],
                                      domain->getMin()[d], domain->getMax()[d],
                                      stretchings[d], jacobs[d]);

        for (unsigned int j = 0; j < num_elem; j++) {

          grid->getVectorIndex(j, indices);
          _Tp jacobian_j = 1.0;

          /***
           * Here we add the jacobian in the whole computation... as for infinite and semi-infinite intervals we are
           * not dealing with linear transformations, the jacobian HAS to be computed at each grid point.... Luckily,
           * it is already available as data stored in the Domai1D objects.
           */
          for (unsigned int d = 0; d < indices.size(); d++) {
            jacobian_j *= (_Tp) jacobs[d][indices[d]];
            coords[d] = stretchings[d][indices[d]];
          }

          f_values[j] = grid->eval(coords) * jacobian_j;
        }

      } else { // not bad stretching

        for (unsigned int j = 0; j < num_elem; j++) {
          grid->getVectorIndex(j, indices);
          _Tp jacobian_j = 1.0;

          /***
           * Here we add the jacobian in the whole computation... as for infinite and semi-infinite intervals we are
           * not dealing with linear transformations, the jacobian HAS to be computed at each grid point.... Luckily,
           * it is already available as data stored in the Domai1D objects.
           */
          for (unsigned int d = 0; d < indices.size(); d++) {
            jacobian_j *=
              (_Tp) domains[d].axisJacobian()[indices[d]];
          }

          f_values[j] = grid->getElementVector()[j] *
                        jacobian_j; // simply pick out the value on the grid...
        }

      }

    } else { // f !=NULL

      if (badstretching) {

        std::vector<std::vector<double> > stretchings(dim,
            std::vector<double>());
        std::vector<std::vector<double> > jacobs(dim,
            std::vector<double>());

        for (int d = 0; d < dim; d++)
          stretching->get1DStretching(grid->getLevels()[d],
                                      domain->getMin()[d], domain->getMax()[d],
                                      stretchings[d], jacobs[d]);

        for (unsigned int j = 0; j < num_elem; j++) {
          grid->getVectorIndex(j, indices);
          _Tp jacobian_j = 1.0;

          for (unsigned int d = 0; d < indices.size(); d++) {
            jacobian_j *= (_Tp) jacobs[d][indices[d]];
            coords[d] = stretchings[d][indices[d]];
          }

          f_values[j] = f(coords) * jacobian_j;
        }

      } else { // not bad stretching
        for (unsigned int j = 0; j < num_elem; j++) {
          grid->getVectorIndex(j, indices);
          grid->getCoords(j, coords);
          _Tp jacobian_j = 1.0;

          for (unsigned int d = 0; d < indices.size(); d++) {
            jacobian_j *=
              (_Tp) domains[d].axisJacobian()[indices[d]];
          }

          f_values[j] = f(coords) *
                        jacobian_j; // simply pick out the value on the grid...
        }
      }

    }

  }

};
}

#endif /* ABSTRACTQUADRATURE_HPP_ */
