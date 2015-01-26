/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Sam Maurus (MA thesis)

#ifndef HESTONPARABOLICPDESOLVERSYSTEMEUROAMER_HPP
#define HESTONPARABOLICPDESOLVERSYSTEMEUROAMER_HPP

#include "base/grid/Grid.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/GridStorage.hpp"
#include "pde/operation/OperationParabolicPDESolverSystemDirichlet.hpp"
#include "finance/tools/Hedging.hpp"

namespace sg {
  namespace finance {

    /**
     * This class implements the ParabolicPDESolverSystem for the BlackScholes
     * Equation.
     *
     * Here European or American Options with fix Dirichlet boundaries are solved.
     */
    class HestonParabolicPDESolverSystemEuroAmer : public sg::pde::OperationParabolicPDESolverSystemDirichlet {
      protected:

        /// the riskfree interest rate
        double r;

        /// Operator on the boundary grid corresponding to the A operator in the thesis
        sg::base::OperationMatrix* OpABound;

        /// Operator on the inner grid corresponding to the A operator in the thesis
        sg::base::OperationMatrix* OpAInner;

        /// Operator on the boundary grid corresponding to the B operator in the thesis
        sg::base::OperationMatrix* OpBBound;

        /// Operator on the inner grid corresponding to the B operator in the thesis
        sg::base::OperationMatrix* OpBInner;

        /// Operator on the boundary grid corresponding to the C operator in the thesis
        sg::base::OperationMatrix* OpCBound;

        /// Operator on the inner grid corresponding to the C operator in the thesis
        sg::base::OperationMatrix* OpCInner;

        /// Operator on the boundary grid corresponding to the D operator in the thesis
        sg::base::OperationMatrix* OpDBound;

        /// Operator on the inner grid corresponding to the D operator in the thesis
        sg::base::OperationMatrix* OpDInner;

        /// Operator on the boundary grid corresponding to the E operator in the thesis
        sg::base::OperationMatrix* OpEBound;

        /// Operator on the inner grid corresponding to the E operator in the thesis
        sg::base::OperationMatrix* OpEInner;

        /// Operator on the boundary grid corresponding to the F operator in the thesis
        sg::base::OperationMatrix* OpFBound;

        /// Operator on the inner grid corresponding to the F operator in the thesis
        sg::base::OperationMatrix* OpFInner;

        /// Operator on the boundary grid corresponding to the G operator in the thesis
        sg::base::OperationMatrix* OpGBound;

        /// Operator on the inner grid corresponding to the G operator in the thesis
        sg::base::OperationMatrix* OpGInner;

        /// Operator on the boundary grid corresponding to the H operator in the thesis
        sg::base::OperationMatrix* OpHBound;

        /// Operator on the inner grid corresponding to the H operator in the thesis
        sg::base::OperationMatrix* OpHInner;

        /// Operator on the boundary grid corresponding to the K operator in the thesis
        sg::base::OperationMatrix* OpKBound;

        /// Operator on the inner grid corresponding to the K operator in the thesis
        sg::base::OperationMatrix* OpKInner;

        /// Operator on the boundary grid corresponding to the X operator in the thesis
        sg::base::OperationMatrix* OpXBound;

        /// Operator on the inner grid corresponding to the X operator in the thesis
        sg::base::OperationMatrix* OpXInner;

        /// Operator on the boundary grid corresponding to the Y operator in the thesis
        sg::base::OperationMatrix* OpYBound;

        /// Operator on the inner grid corresponding to the Y operator in the thesis
        sg::base::OperationMatrix* OpYInner;

        /// Operator on the boundary grid corresponding to the W operator in the thesis
        sg::base::OperationMatrix* OpWBound;

        /// Operator on the inner grid corresponding to the W operator in the thesis
        sg::base::OperationMatrix* OpWInner;

        /// Operator on the boundary grid corresponding to the Z operator in the thesis
        sg::base::OperationMatrix* OpZBound;

        /// Operator on the inner grid corresponding to the Z operator in the thesis
        sg::base::OperationMatrix* OpZInner;

        /// Pointer to the vector containing the volatility of volatility values
        sg::base::DataVector* volvols;

        /// Pointer to the kappas
        sg::base::DataVector* kappas;

        /// Pointer to the thetas
        sg::base::DataVector* thetas;

        /// Pointer to the rhos
        sg::base::DataMatrix* hMatrix;

        /// Coefficient collection for the D operator. Only one custom up/down operator, so it's just a vector of coefficients.
        sg::base::DataVector* dCoeff;

        /// Coefficient collection for the E operator. Only one custom up/down operator, so it's just a vector of coefficients.
        sg::base::DataVector* eCoeff;

        /// Coefficient collection for the F operator. Only one custom up/down operator, so it's just a vector of coefficients.
        sg::base::DataVector* fCoeff;

        /// Coefficient collection for the G operator. Only one custom up/down operator, so it's just a vector of coefficients.
        sg::base::DataVector* gCoeff;

        /// Coefficient collection for the Z operator. Only one custom up/down operator, so it's just a vector of coefficients.
        sg::base::DataVector* zCoeff;

        /// Coefficient collection for the B operator. Two custom up/down operators, so it's a matrix of coefficients.
        sg::base::DataMatrix* bCoeff;

        /// Coefficient collection for the C operator. Two custom up/down operators, so it's a matrix of coefficients.
        sg::base::DataMatrix* cCoeff;

        /// Coefficient collection for the H operator. Two custom up/down operators, so it's a matrix of coefficients.
        sg::base::DataMatrix* hCoeff;

        /// Coefficient collection for the X operator. Two custom up/down operators, so it's a matrix of coefficients.
        sg::base::DataMatrix* xCoeff;

        /// Coefficient collection for the Y operator. Two custom up/down operators, so it's a matrix of coefficients.
        sg::base::DataMatrix* yCoeff;

        /// Coefficient collection for the W operator. Two custom up/down operators, so it's a matrix of coefficients.
        sg::base::DataMatrix* wCoeff;

        /// Up/down four op dims
        double**** kCoeff;

        /// Whether or not to use coarsening between timesteps in order to reduce gridsize
        bool useCoarsen;

        /// Adaptive mode during solving Heston Equation: coarsen, refine, coarsenNrefine
        std::string adaptSolveMode;

        /// Number of points the are coarsened in each coarsening-step
        int numCoarsenPoints;

        /// Threshold used to decide if a grid point should be deleted
        double coarsenThreshold;

        /// Threshold used to decide if a grid point should be refined
        double refineThreshold;

        /// Refine mode during solving the Heston Equation: classic or maxLevel
        std::string refineMode;

        /// MaxLevel max. Level of refinement
        sg::base::GridIndex::level_type refineMaxLevel;

        /// The algorithmic dimensions used in this system
        std::vector<size_t> HestonAlgoDims;

        /// The number of assets (half the number of problem dimensions)
        size_t nAssets;

        /// Stores the number of executed timesteps
        size_t nExecTimesteps;

        /// The strike of the current option
        double dStrike;

        /// The type of the current option
        std::string option_type;

        /// Store whether log coordinates are used
        bool b_log_transform;

        virtual void applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result);
        virtual void applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result);
        virtual void applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result);
        virtual void applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result);

        /**
         * Builds the coefficient object for the (non-log-transformed) D operator.
         * Note that only vanilla options are supported for this (non-transformed) case.
         */
        void buildDCoefficients();

        /**
         * Builds the coefficient object for the (non-log-transformed) F operator.
         * Note that only vanilla options are supported for this (non-transformed) case.
         */
        void buildFCoefficients();

        /**
         * Builds the coefficient object for the G operator.
         * Note that only vanilla options are supported for this (non-transformed) case.
         */
        void buildGCoefficients();

        /**
         * Builds the coefficient object for the X operator.
         * Note that only vanilla options are supported for this (non-transformed) case.
         */
        void buildXCoefficients();

        /**
         * Builds the coefficient object for the Y operator.
         * Note that only vanilla options are supported for this (non-transformed) case.
         */
        void buildYCoefficients();

        /**
         * Builds the coefficient object for the W operator.
         * Note that only vanilla options are supported for this (non-transformed) case.
         */
        void buildWCoefficients();

        /**
         * Builds the coefficient object for the Z operator.
         * Note that only vanilla options are supported for this (non-transformed) case.
         */
        void buildZCoefficients();

        /**
         * Builds the coefficients object for the B operator.
         * Operator B has two custom 1D operators, so the coefficients object is a matrix.
         * The operator has a non-zero coefficient for each pairing of a stock and its variance, e.g. elements (1,2), (3,4), (5,6) etc. (zero-based array elements (0,1), (2,3), (4,5) etc.)
         */
        void buildBCoefficientsLogTransform();


        /**
         * Builds the coefficients object for the C operator.
         * Operator C has two custom 1D operators, so the coefficients object is a matrix.
         * The operator has a non-zero coefficient for each pairing of a stock and its variance, e.g. elements (1,2), (3,4), (5,6) etc. (zero-based array elements (0,1), (2,3), (4,5) etc.)
         */
        void buildCCoefficientsLogTransform();

        /**
         * Builds the coefficients object for the D operator.
         * The D operator is oneOpDim, so there is only a vector of coefficients to set.
         * Only the variance dimensions are active in this operator, i.e. elements 2, 4, 6, 8 etc. (zero-based array elements 1, 3, 5, 7 etc.)
         */
        void buildDCoefficientsLogTransform();

        /**
         * Builds the coefficients object for the E operator.
         * The E operator is oneOpDim, so there is only a vector of coefficients to set.
         * Only the stock-price dimensions are active in this operator, i.e. elements 1, 3, 5, 7 etc. (zero-based array elements 0, 2, 4, 6 etc.)
         */
        void buildECoefficientsLogTransform();

        /**
         * Builds the coefficients object for the F operator.
         * The F operator is oneOpDim, so there is only a vector of coefficients to set.
         * Only the variance elements are active in this operator, i.e. elements 2, 4, 6, 8 etc. (zero-based array elements 1, 3, 5, 7 etc.)
         */
        void buildFCoefficientsLogTransform();

        /**
         * Builds the coefficients object for the G operator.
         * The G operator is oneOpDim, so there is only a vector of coefficients to set.
         * Only the variance elements are active in this operator, i.e. elements 2, 4, 6, 8 etc. (zero-based vector elements 1, 3, 5, 7 etc.)
         */
        void buildGCoefficientsLogTransform();

        /**
         * Builds the coefficients object for the H operator.
         * Operator H has two custom 1D operators, so the coefficients object is a matrix.
         * The operator has a non-zero coefficient for each pairing of a stock and its variance, e.g. elements (1,2), (3,4), (5,6) etc. (zero-based array elements (0,1), (2,3), (4,5) etc.)
         */
        void buildHCoefficientsLogTransform();

        /**
         * Builds the coefficients object for the K operator.
         * This operator is more exciting. It has four custom 1D operators, so the coefficients object is a 4d array.
         * The operator has a non-zero coefficient for each quad-pairing of a stock and its variance, and a DIFFERENT stock and its variance, e.g. for three assets: (1,2,3,4), (1,2,5,6), (3,4,5,6). (zero-based array elements (0,1,2,3), (0,1,4,5), (2,3,4,5))
         */
        void buildKCoefficientsLogTransform();

        /**
         * Utility method for creating a 4d array of equal size in each dimension. The array is allocated dynamically based on the provided size.
         * @param dimSize The number of elements in each dimension of the array.
         * @param array Pointer to the first element of the 4d array.
         */
        void create4dEqualDimSizeArray(size_t dimSize, double**** * array);


        /**
         * Utility method to free the memory allocated to a 4d array of equal size in each dimension.
         * @param dimSize The number of elements in each dimension of the array.
         * @param array Pointer to the first element of the 4d array.
         */
        void delete4dEqualDimSizeArray(size_t dimSize, double**** * array);

        /**
         * Utility method to set all values in a 4d array of equal size to one value (e.g. this method can be used to set all the values to zero).
         * @param dimSize The number of elements in each dimension of the array.
         * @param array Pointer to the first element of the 4d array.
         * @param value Value that each of the array elements will be set to.
         */
        void setAll4dEqualDimSizeArray(size_t dimSize, double**** * array, double value);

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid the sparse grid on which the solution is based
         * @param alpha basis function coefficients
         * @param thetas collection of theta values, one for each asset
         * @param volvols collection of the volatility of the volatility values, one for each asset
         * @param kappas collection of the kappa values, one for each asset
         * @param rho collection of the correlation values between the Wiener processes. The matrix size must be twice the number of assets to encompass the correlation between the stock price and variance processes.
         * @param r market risk-free interest rate
         * @param TimestepSize size of the timestep used in the ODE solver
         * @param OperationMode the solver used for solving the ODE, e.g. "CrNic"
         * @param dStrike strike price of the option
         * @param option_type option flavour specifier, e.g. "std_euro_call"
         * @param bLogTransform true if log-transformed stock-price coordinates are used, false otherwise. The variance coordinates are always linear.
         * @param useCoarsen true if coarsening should be used
         * @param coarsenThreshold threshold value from which coarsening is applied
         * @param adaptSolveMode specifies the type of adaptivity used
         * @param numCoarsenPoints specifies the number of points used in the coarsening algorithm
         * @param refineThreshold threshold value from which refinement is applied
         * @param refineMode mode of refinement
         * @param refineMaxLevel maximum sparse grid level above which no refinement is performed
         */
        HestonParabolicPDESolverSystemEuroAmer(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, sg::base::DataVector& thetas, sg::base::DataVector& volvols,
                                               sg::base::DataVector& kappas,
                                               sg::base::DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
                                               double dStrike, std::string option_type,
                                               bool bLogTransform = false, bool useCoarsen = false, double coarsenThreshold = 0.0, std::string adaptSolveMode = "none",
                                               int numCoarsenPoints = -1, double refineThreshold = 0.0, std::string refineMode = "classic", sg::base::GridIndex::level_type refineMaxLevel = 0);

        /**
         * Destructor. Just does the usual...releases all allocated memory.
         */
        virtual ~HestonParabolicPDESolverSystemEuroAmer();

        virtual void finishTimestep();

        /**
         * Coarsens and refines the sparse grid based on the relevant thresholds.
         * @param isLastTimestep true if the current timestep is the last. If this is true, integrated coarsening is not performed.
         */
        virtual void coarsenAndRefine(bool isLastTimestep = false);

        void startTimestep();
    };

  }

}

#endif /* HESTONPARABOLICPDESOLVERSYSTEMEUROAMER_HPP */
