/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP
#define OPERATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP

#include "pde/operation/OperationEllipticPDESolverSystem.hpp"
#include "base/grid/common/DirichletUpdateVector.hpp"
#include "base/grid/common/DirichletGridConverter.hpp"

#include <string>

namespace sg {
  namespace pde {

    /**
     * Defines a System that is used to solve elliptic partial
     * differential equations. So an instance of this class has to pass to
     * any SLE Solver used in SGpp.
     *
     * \f$L \vec{u} = rhs\f$
     *
     * L: space discretization (L-Operator)
     * rhs: right hand side
     *
     * This class is a specialized version of OperationEllipticPDESolverSystem which
     * exploits Dirichlet boundary conditions. Since there are no degrees of freedom
     * on on the boundaries the iterative solver (CG or BiCGSTAB) has only to take
     * inner grid points into account.
     *
     * The inner grid is constructed during the constructor call!
     */
    class OperationEllipticPDESolverSystemDirichlet : public OperationEllipticPDESolverSystem {
      protected:
        /// Pointer to the alphas (ansatzfunctions' coefficients; inner points only)
        sg::base::DataVector* alpha_inner;
        /// Routine to modify the boundaries/inner points of the grid
        sg::base::DirichletUpdateVector* BoundaryUpdate;
        /// Class that allows a simple conversion between a grid with and a without boundary points
        sg::base::DirichletGridConverter* GridConverter;
        /// Pointer to the inner grid object
        sg::base::Grid* InnerGrid;
        /// rhs for the inner grid
        sg::base::DataVector* rhs_inner;

        /**
         * applies the PDE's system matrix, on complete grid - with boundaries
         *
         * @param alpha the coefficients of the sparse grid's ansatzfunctions
         * @param result reference to the sg::base::DataVector into which the result is written
         */
        virtual void applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;

        /**
         * applies the PDE's system matrix, on inner grid only
         *
         * @param alpha the coefficients of the sparse grid's ansatzfunctions
         * @param result reference to the sg::base::DataVector into which the result is written
         */
        virtual void applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;

        /**
         * Use this function in order to obtain the system for
         * solving an elliptical PDE on Sparse Grids with an extern
         * solver (e.g. Intel's MKL). The matrix is written into the
         * mtxString in Matrix Market format (http://math.nist.gov/MatrixMarket/formats.html).
         *
         * @param mtxString reference to string-object into which the serialized matrix is stored
         * @param complete indicates whether the matrix of the complete system (including boundaries) (=true) should be generated or not
         *
         * @return the number of non zeros in the system matrix
         */
        size_t getMatrix(std::string& mtxString, bool complete);

        /**
         * Use this function in order to obtain the system for
         * solving an elliptical PDE on Sparse Grids with an extern
         * solver (e.g. Intel's MKL). The matrix is written into the
         * mtxString in Matrix Market format (http://math.nist.gov/MatrixMarket/formats.html).
         *
         * Only the systemmatrix's diagonal is exported
         *
         * @param mtxString reference to string-object into which the serialized matrix is stored
         * @param complete indicates whether the matrix of the complete system (including boundaries) (=true) should be generated or not
         */
        void getMatrixDiagonal(std::string& mtxString, bool complete);

        /**
         * Use this function in order to obtain the system for
         * solving an elliptical PDE on Sparse Grids with an extern
         * solver (e.g. Intel's MKL). The matrix is written into the
         * mtxString in Matrix Market format (http://math.nist.gov/MatrixMarket/formats.html).
         *
         * The system's row sum is exported as a diagonal matrix
         *
         * @param mtxString reference to string-object into which the serialized matrix is stored
         * @param complete indicates whether the matrix of the complete system (including boundaries) (=true) should be generated or not
         */
        void getMatrixDiagonalRowSum(std::string& mtxString, bool complete);

      public:
        /**
         * Constructor
         *
         * @param SparseGrid the grid, for which the system should be solved
         * @param rhs the right hand side of the corresponding system
         */
        OperationEllipticPDESolverSystemDirichlet(sg::base::Grid& SparseGrid, sg::base::DataVector& rhs);

        /**
         * Destructor
         */
        virtual ~OperationEllipticPDESolverSystemDirichlet();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual sg::base::DataVector* generateRHS();

        /**
         * gets a pointer to the sparse grids coefficients used in the CG method to solve
         * one timestep. This is useful because (direchlet) boundaries can be skipped when
         * solving the system.
         *
         * @return alpha vector for CG method
         */
        virtual sg::base::DataVector* getGridCoefficientsForCG();

        /**
         * Gets the solution for the complete grid
         *
         * @param Solution sg::base::DataVector that must have a dimension equal to the bound's grid dimension, the result is written to Solution
         * @param SolutionInner Solution on the inner grid
         */
        virtual void getSolutionBoundGrid(sg::base::DataVector& Solution, sg::base::DataVector& SolutionInner);

        /**
         * Use this function in order to obtain the system for
         * solving an elliptical PDE on Sparse Grids with an extern
         * solver (e.g. Intel's MKL). The matrix is written into the
         * mtxString in Matrix Market format (http://math.nist.gov/MatrixMarket/formats.html).
         *
         * For this function the matrix including the boundary ansatzfunctions
         * is generated
         *
         * @param mtxString reference to string-object into which the serialized matrix is stored
         *
         * @return the number of non zeros in the system matrix
         */
        size_t getCompleteMatrix(std::string& mtxString);

        /**
         * Use this function in order to obtain the system for
         * solving an elliptical PDE on Sparse Grids with an extern
         * solver (e.g. Intel's MKL). The matrix is written into the
         * mtxString in Matrix Market format (http://math.nist.gov/MatrixMarket/formats.html).
         *
         * For this function the matrix excluding the boundary ansatzfunctions
         * is generated
         *
         * @param mtxString reference to string-object into which the serialized matrix is stored
         *
         * @return the number of non zeros in the system matrix
         */
        size_t getInnerMatrix(std::string& mtxString);

        /**
         * Use this function in order to obtain the system for
         * solving an elliptical PDE on Sparse Grids with an extern
         * solver (e.g. Intel's MKL). The matrix is written into the
         * mtxString in Matrix Market format (http://math.nist.gov/MatrixMarket/formats.html).
         *
         * For this function the matrix excluding the boundary ansatzfunctions
         * is generated
         *
         * Only the systemmatrix's diagonal is exported
         *
         * @param mtxString reference to string-object into which the serialized matrix is stored
         *
         * @return the number of non zeros in the system matrix
         */
        void getInnerMatrixDiagonal(std::string& mtxString);

        /**
         * Use this function in order to obtain the system for
         * solving an elliptical PDE on Sparse Grids with an extern
         * solver (e.g. Intel's MKL). The matrix is written into the
         * mtxString in Matrix Market format (http://math.nist.gov/MatrixMarket/formats.html).
         *
         * For this function the matrix excluding the boundary ansatzfunctions
         * is generated
         *
         * The systemmatrix's row sum is exported as a diagonal matrix
         *
         * @param mtxString reference to string-object into which the serialized matrix is stored
         *
         * @return the number of non zeros in the system matrix
         */
        void getInnerMatrixDiagonalRowSum(std::string& mtxString);
    };

  }
}

#endif /* OPERATIONELLIPTICPDESOLVERMATRIXDIRICHLET_HPP */
