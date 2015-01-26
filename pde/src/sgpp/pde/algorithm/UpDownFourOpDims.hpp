/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Sam Maurus (MA thesis)

#ifndef UPDOWNFOUROPDIMS_HPP
#define UPDOWNFOUROPDIMS_HPP

#include <vector>
#include <map>

#include "base/grid/GridStorage.hpp"
#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

#ifndef TASKS_PARALLEL_UPDOWN
#define TASKS_PARALLEL_UPDOWN 4
#endif

namespace sg {
  namespace pde {

    /**
     * Implements the Up/Down scheme with four dimensions with special operations: i,j,k,l
     *
     * @version $HEAD$
     */
    class UpDownFourOpDims: public sg::base::OperationMatrix {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         * @param coef 4d tensor that contains the constant coefficients of this operation
         */
        UpDownFourOpDims(sg::base::GridStorage* storage, double**** * coef);

        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        UpDownFourOpDims(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~UpDownFourOpDims();


        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

      protected:
        typedef sg::base::GridStorage::grid_iterator grid_iterator;

        /// Function pointer type. This is used in fnMap to map the particular dimension situation to the relevant method handler.
        typedef void (sg::pde::UpDownFourOpDims::*MFP)(sg::base::DataVector&, sg::base::DataVector&, size_t, size_t, size_t, size_t, size_t);

        /// Pointer to the grid's storage object
        sg::base::GridStorage* storage;
        /// Pointer to the coefficients of this bilinear form
        double**** coefs;
        /// algorithmic dimensions, operator is applied in this dimensions
        const std::vector<size_t> algoDims;
        /// number of algorithmic dimensions
        const size_t numAlgoDims_;
        /// max number of parallel stages (dimension recursive calls)
        static const size_t maxParallelDims_ = TASKS_PARALLEL_UPDOWN;

        /// Map of integer to function pointer. This is used to map the dimension situation to the relevant method handler.
        std::map<size_t, MFP> fnMap;

        /**
         * Utility method to generate the fnMap member for mappings.
         */
        void generateMap();

        /**
         * Recursive procedure for updown, parallel version using OpenMP 3
         *
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         * @param dim the current dimension
         * @param op_dim_one the dimension in which to use the first gradient
         * @param op_dim_two the dimension in which to use the second gradient
         * @param op_dim_three the dimension in which to use the third gradient
         * @param op_dim_four the dimension in which to use the fourth gradient
         */
        void updown(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        // Unidirectional
        void specialOpUnidirectional(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        // Singles
        void specialOpOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
        void specialOpTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
        void specialOpThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
        void specialOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
        void specialOpX(sg::base::DataVector& alpha, sg::base::DataVector& result, void (sg::pde::UpDownFourOpDims::*pt2UpFunc)(sg::base::DataVector&, sg::base::DataVector&, size_t), void (sg::pde::UpDownFourOpDims::*pt2DownFunc)(sg::base::DataVector&, sg::base::DataVector&, size_t), size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        // Doubles

        /**
         * If the current dimension is equal to the both special operation dimensions one and two.
         *
         * @param alpha the coefficients of the grid points
         * @param result the result of the operations
         * @param dim the current dimension in the recursion
         * @param op_dim_one the dimension in which to use the first gradient
         * @param op_dim_two the dimension in which to use the second gradient
         * @param op_dim_three the dimension in which to use the third gradient
         * @param op_dim_four the dimension in which to use the fourth gradient
         */
        void specialOpOneAndOpTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        /**
         * If the current dimension is equal to the both special operation dimensions one and three.
         * For an explanation of the parameters of this method, see the documentation for the method specialOpOneAndOpTwo in this class.
         */
        void specialOpOneAndOpThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        /**
         * If the current dimension is equal to the both special operation dimensions one and four.
         * For an explanation of the parameters of this method, see the documentation for the method specialOpOneAndOpTwo in this class.
         */
        void specialOpOneAndOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        /**
         * If the current dimension is equal to the both special operation dimensions two and three.
         * For an explanation of the parameters of this method, see the documentation for the method specialOpOneAndOpTwo in this class.
         */
        void specialOpTwoAndOpThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        /**
         * If the current dimension is equal to the both special operation dimensions two and four.
         * For an explanation of the parameters of this method, see the documentation for the method specialOpOneAndOpTwo in this class.
         */
        void specialOpTwoAndOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        /**
         * If the current dimension is equal to the both special operation dimensions three and four.
         * For an explanation of the parameters of this method, see the documentation for the method specialOpOneAndOpTwo in this class.
         */
        void specialOpThreeAndOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        /**
         * If the current dimension is equal to the all special operation dimensions one, two and three.
         * For an explanation of the parameters of this method, see the documentation for the method specialOpOneAndOpTwo in this class.
         */
        void specialOpOneAndOpTwoAndOpThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        /**
         * If the current dimension is equal to the all special operation dimensions one, two and four.
         * For an explanation of the parameters of this method, see the documentation for the method specialOpOneAndOpTwo in this class.
         */
        void specialOpOneAndOpTwoAndOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        /**
         * If the current dimension is equal to the all special operation dimensions one, three and four.
         * For an explanation of the parameters of this method, see the documentation for the method specialOpOneAndOpTwo in this class.
         */
        void specialOpOneAndOpThreeAndOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        /**
         * If the current dimension is equal to the all special operation dimensions two, three and four.
         * For an explanation of the parameters of this method, see the documentation for the method specialOpOneAndOpTwo in this class.
         */
        void specialOpTwoAndOpThreeAndOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

        /**
         * If the current dimension is equal to the all special operation dimensions one, two, three and four.
         * For an explanation of the parameters of this method, see the documentation for the method specialOpOneAndOpTwo in this class.
         */
        void specialOpOneAndOpTwoAndOpThreeAndOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);


        /**
         * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
         * Applies the up-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i < l_j} \alpha_j \phi_j(x) dx.\f]
         *
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         * @param dim dimension in which to apply the up-part
         */
        virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
         * Applies the down-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
         *
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         * @param dim dimension in which to apply the down-part
         */
        virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to i.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        virtual void downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to i.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to j.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to j.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to k.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to k.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to i and j.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to i and j.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to i and k.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to i and k.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to i and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to i and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to j and k.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to j and k.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to j and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to j and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to k and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to k and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to i and j and k.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to i and j and k.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to i and j and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to i and j and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to i and k and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to i and k and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D down if the current dim is equal to j and k and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to j and k and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;


        /**
         * 1D down if the current dim is equal to i and j and k and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * 1D up if the current dim is equal to i and j and k and l.
         * For an explanation of the parameters of this method, see the documentation for the method downOpDimOne in this class.
         */
        virtual void upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
    };

  }
}

#endif /* UPDOWNFOUROPDIMS_HPP */
