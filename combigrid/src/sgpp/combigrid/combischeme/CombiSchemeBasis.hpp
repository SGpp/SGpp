/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef COMBICOMBISCHEMEBASIS_HPP_
#define COMBICOMBISCHEMEBASIS_HPP_

#include <sgpp/combigrid/utils/combigrid_ultils.hpp>
//#include <sgpp/combigrid/combigridkernel/CombiGridKernel.hpp>

namespace combigrid {

  /** base class for any combi scheme. From this class should all the scheme classes derived e.g. S-CT, TS-CT, ... <br>
   * */
  class CombiSchemeBasis {

    public:

      /** Empty constructor */
      CombiSchemeBasis() {
        dim_ = 0;
        levels_vector_.resize(0);
        levels_.resize(dim_, 0);
        cofficients_.resize(0);
      }

      /** Ctor */
      CombiSchemeBasis(int dim, int level) {
        dim_ = dim;
        levels_vector_.resize(0);
        levels_.resize(dim_, level);
        cofficients_.resize(0);
      }

      /** Ctor */
      CombiSchemeBasis(int dim, const std::vector<int>& levels) {
        dim_ = dim;
        levels_vector_.resize(0);
        levels_ = levels;
        cofficients_.resize(0);
      }

      /** return the dimension of the combi scheme */
      inline int getDim() const {
        return dim_;
      }

      /** number of subsapces */
      inline int getNrSapces() const {
        return static_cast<int>(levels_vector_.size());
      }

      /** returns the level vector for one subspace */
      inline const std::vector<int>& getLevel(int i) const {
        return levels_vector_[i];
      }

      inline const std::vector<std::vector<int> >& getLevels() const {
        return levels_vector_;
      }

      /** the levels of the full grid which we extrapolate with the combi scheme*/
      inline const std::vector<int>& getMaxLevel() const {
        return levels_;
      }

      /** returns the coefficient for one subspace */
      inline double getCoef(int i) const {
        return cofficients_[i];
      }

      /** returns the coefficient for one subspace */
      inline std::vector<double> getCoef() const {
        return cofficients_;
      }

      /** Method adding a new full grid.
       * Coefficients and levels are updated and the indices of the changed
       * levels are returned.
       */
      std::vector<int> updateScheme(std::vector<std::vector<int> > levelsNew,
                                    std::vector<double> coef);

      void setCoef(std::vector<double> newCoef);

      void setCoef(int i, double newCoef);

    protected:

      /** function which removes the duplicate spaces */
      void removeDuplicates();

      /** the dimension of the scheme */
      int dim_;

      /** the level vector for each space */
      std::vector<std::vector<int> > levels_vector_;

      /** the coefficients for the spaces */
      std::vector<double> cofficients_;
      /** the levels of the full grid which we*/
      std::vector<int> levels_;

  };

}

#endif /* COMBICOMBISCHEMEBASIS_HPP_ */
