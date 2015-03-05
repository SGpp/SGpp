// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATORFG_HPP_
#define OPERATORFG_HPP_

#include "combigrid.hpp"
#include <sgpp/combigrid/utils/combigrid_ultils.hpp>
#include <sgpp/combigrid/combigrid/AbstractCombiGrid.hpp>

namespace combigrid {

  /** abstract class which is the interface between a problem specific and the multigird (or other)
   * solving algorithm.<br>
   * It specifies an abstract operator acting on a full grid. This operator is problem specific
   * and here only the interface is defined, through which a solver (e.g. multigrid) can use the operator.
   * The main interacting component of the interface are the vectors (unknown vector and right hand side vector).
   * It is IMPORTANT to note that there can be more than one unknowns per node. <br>
   * In this way the vectors have [nrNodes*nrSpace] elements. <br>
   * The index CONVENTION is : [nodeIndex*nrSpace + spaceIndex] <br>
   * IT IS IMPORTANT THAT ALL THE OPERATORS USE THIS CONVENTION !!! <br>
   * IT IS IMPORTANT WE NUMBER THE UNKOWNS ACCORDING TO THE CONVENTION IN THE "combigrid::FullGridD" CLASS !!! <br>*/
  class OperatorFG {

    public:

      /** Ctor
       * @param fg [IN] the full grid on which the operator operates
       * @param nrSpace [IN] number of unknowns per node */
      OperatorFG(const FullGridD* fg , int nrSpace ): fg_(fg), nrSpace_(nrSpace) {}

      /** method to create a new operator which acts on a given full grid
       * @param fg [IN] */
      virtual OperatorFG* factory(const FullGridD* fg) const = 0;

      virtual ~OperatorFG() {
        delete fg_;
      }

      /** method to get the right hand side vector for one level
       * @param rhs [OUT] vector with the right hand side values
       * @param nrSpace [OUT] number of unknowns per node */
      virtual void getRHS(std::vector<double>& rhs , int& nrSpace) const = 0;

      /** method to multiply one vector with the operator matrix.
       * @param inVect [IN] input vector
       * @param outVect [OUT] output vector (here will be the result A*inVect = outVect)*/
      virtual void multiplyVector(std::vector<double>& inVect , std::vector<double>& outVect) const = 0;

      /** method to make Smoothing iterations (the number of iterations is specified by the user).
       * @param nrIt [IN] number of smoothing iterations to execute
       * @param u [IN/OUT] the unknown vector
       * @param rhs [IN] the right hand side for the iterations */
      virtual void doSmoothing(int nrIt ,
                               std::vector<double>& u, std::vector<double>& rhs) const = 0;

      // todo: introduce function for initial guess for solutions

      /** returns the pointer to the full grid */
      inline const FullGridD* getFG() const {
        return fg_;
      }

      /** returns the number of spaces (number of unknowns per node)*/
      inline int getNrSpace() const {
        return nrSpace_;
      }

    private:

      /** full grid on which the operator operates */
      const FullGridD* fg_;

      /** number of spaces; number of unknowns per node <br>.
       * This can vary from case to case, from problem to problem. */
      int nrSpace_;


  };
}

#endif /* OPERATORFG_HPP_ */