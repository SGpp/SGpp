// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRIDLEVELVECTOR_HPP_
#define COMBIGRIDLEVELVECTOR_HPP_

#include <sgpp/combigrid/utils/combigrid_ultils.hpp>

namespace combigrid {

  /**Levelvector allowing computations done in hegland.2003.adaptive equation (16)
   * for the computation of required grid to create a valid combination solution.
   * Represents the P_alpha in that paper. It is thus also including the respective
   * coefficient of the full grid level vector. */
  class CombigridLevelVector {
    public:
      /**Creates CombigridLevelVector containing the levels of one grid. The
       * dimension of the grid is determined by the size of the input vector */
      CombigridLevelVector(std::vector<int> level);

      /**Creates an CombiGridLevelVector of a certain dimension, with all levels
       *  having the maximum size. This can be used as a unity operator
       *  for the multiplication.*/
      CombigridLevelVector(int dim);
      /** Creates an empty CombigridLevelVector*/
      CombigridLevelVector() {
        levelVec_.resize(0);
        coef_.resize(0);
      }
      ;
      /** Creates an CombigridLevelVector with a certain coefficient. Might only
       * be used in rare cases.*/
      CombigridLevelVector(std::vector<int> level, double coef);
      /** Creates a complete set of fullgrids including their coefficients.*/
      CombigridLevelVector(std::vector<std::vector<int> > input,
                           std::vector<double> coef);

      /** Returns all full Grid level vectors*/
      std::vector<std::vector<int> > getLevelVec() const {
        return levelVec_;
      }

      /** Returns a single Full Grid Level vector containing the level of the
       * full grid in each dimension.
       */
      std::vector<int> getLevelVecSingle(int i) const {
        return levelVec_[i];
      }

      /** Returns a vector of coefficients in the respective ordering*/
      std::vector<double> getCoef() const {
        return coef_;
      }

      /** Returns the dimensionality of the grids*/
      int getDim() const {
        return int( levelVec_[0].size() );
      }

      /** Returns the number of full grids */
      int getN() const {
        return int( levelVec_.size() );
      }

      CombigridLevelVector& operator=(const CombigridLevelVector& rhs);

      /** Multiplication operator, which returns the minimum in each dimension.*/
      const CombigridLevelVector operator*(const CombigridLevelVector& b) const;

      /** Operator unifying two CombigridLevelVectors. Full Grids of the some size
       * get their coefficients added.
       */
      const CombigridLevelVector operator+(const CombigridLevelVector& b) const;
      /** Operator unifying two CombigridLevelVectors. Full Grids of the some size
       * get their coefficients subtracted.
       */
      const CombigridLevelVector operator-(const CombigridLevelVector& b) const;

      /** Double level vectors get unified.
       * */
      void doAddition();
      /** Print the CombigridLevelVector */
      void printLevelVec();

      /** Splits a CombigridLevelVector in separate level vectors*/
      std::vector<CombigridLevelVector> split();

      /**
       * Function creating the levels of a combigrid containing all subgrids
       * provided in the Combigridlevelvectors
       */
      static CombigridLevelVector getCombiLevels(
        std::vector<CombigridLevelVector> input);
      /**
       * Function creating the levels of a combigrid containing all subgrids
       * provided in the Combigridlevelvectors
       */
      static CombigridLevelVector getCombiLevels(
        std::vector<std::vector<int> > input);
      /**
       * Function creating the levels of a combigrid containing all subgrids
       * provided in the Combigridlevelvectors
       */
      static CombigridLevelVector getCombiLevels(CombigridLevelVector input);

      /** Add new active full grid. The changes of the combischeme are returned */
      CombigridLevelVector getChanges(std::vector<int> newFullGridLevels);

      /** Add a new FullGrid to the current CombigridLevelVec*/
      void update(std::vector<int> newFullGridLevels);

    protected:
      std::vector<std::vector<int> > levelVec_;
      std::vector<double> coef_;
      static const int LEVELMAX;
  };

}

#endif /* COMBIGRIDLEVELVECTOR_HPP_ */