/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Dirk Pflüger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef GETAFFECTEDBASISFUNCTIONS_HPP
#define GETAFFECTEDBASISFUNCTIONS_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"
#include "basis/linear/boundary/linearboundaryBase.hpp"

#include <vector>
#include <utility>

namespace sg {

/**
 * Basic algorithm for getting all affected basis functions.
 * This implicitly assumes a tensor-product approach and local support.
 * No grid points on the border are supported.
 *
 * The main idea behind this algorithm is to spend as few function evaluations as possible.
 * Assume a regular sparse grid level 3 in two dimensions with the sparse grid basis
 * \f$\Phi:=\{\phi_i(x), i=1,\ldots,N\}\f$. Then the tableau of subspaces looks
 * as follows:
 *   \image html GetAffectedBasisFunctions_subspaces.png "Tableau of subspaces for a regular sparse grid level 3"
 * You could evaluate the function \f$ f_N(x) = \sum_{i=1}^N \alpha_i \phi_i(x)\f$ for all basis
 * functions \f$\phi_i(x)\f$, multiply them with the surplus and add them up.
 * In \f$d\f$ dimensions this would lead to \f$N\f$ evaluations of \f$d\f$ one-dimensional basis
 * functions each.
 *
 * A better way is to (recursively) look at each subspace, as only one basis function
 * per subspace can be non-zero (partially disjunct supports):
 *   \image html GetAffectedBasisFunctions_subspaces_affectedBasisFunctions.png "Traversal of subspaces for evaluation"
 * This can be done recursively in both the dimension and the level. In each subspace
 * the basis function concerned can be identified via a few index calculations and
 * evaluated at the given point in the domain.
 *
 * Even better would be to save further function evaluations and to reuse intermediate values obtained by
 * the evaluation of one-dimensional basis functions, see the following figure.
 *   \image html GetAffectedBasisFunctions_subspaces_affectedBasisFunctions_recursive.png "Minimize the number of evaluations" width=10cm
 * Descending recursively in the d-th dimension, one can propagate the value of the intermediate function
 * evaluation for the first d-1 dimensions that have already been looked at.
 *
 * @todo (heinecke) add bounding box support for every grid type
 */
template<class BASIS>
class GetAffectedBasisFunctions {
public:
    GetAffectedBasisFunctions(GridStorage* storage) :
        storage(storage) {
    }

    ~GetAffectedBasisFunctions() {
    }

    /**
     * Returns evaluations of all basis functions that are non-zero at a given evaluation point.
     * For a given evaluation point \f$x\f$, it stores tuples (std::pair) of
     * \f$(i,\phi_i(x))\f$ in the result vector for all basis functions that are non-zero.
     * If one wants to evaluate \f$f_N(x)\f$, one only has to compute
     * \f[ \sum_{r\in\mathbf{result}} \alpha[r\rightarrow\mathbf{first}] \cdot r\rightarrow\mathbf{second}. \f]
     *
     * @param basis a sparse grid basis
     * @param point evaluation point within the domain
     * @param result a vector to store the results in
     */
    void operator()(BASIS& basis, std::vector<double>& point, std::vector<
            std::pair<size_t, double> >& result) {
        GridStorage::grid_iterator working(storage);

        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;

        size_t bits = sizeof(index_type) * 8; // how many levels can we store in a index_type?

        size_t dim = storage->dim();

        index_type* source = new index_type[dim];

        for (size_t d = 0; d < dim; ++d) {
            // This does not really work on grids with borders.
            double temp = floor(point[d] * (1 << (bits - 2))) * 2;
            if (point[d] == 1.0) {
                source[d] = static_cast<index_type> (temp - 1);
            } else {
                source[d] = static_cast<index_type> (temp + 1);
            }

        }

        result.clear();
        rec(basis, point, 0, 1.0, working, source, result);

        delete[] source;

    }

protected:
    GridStorage* storage;

    /**
     * Recursive traversal of the "tree" of basis functions for evaluation, used in operator().
     * For a given evaluation point \f$x\f$, it stores tuples (std::pair) of
     * \f$(i,\phi_i(x))\f$ in the result vector for all basis functions that are non-zero.
     *
     * @param basis a sparse grid basis
     * @param point evaluation point within the domain
     * @param current_dim the dimension currently looked at (recursion parameter)
     * @param value the value of the evaluation of the current basis function up to (excluding) dimension current_dim (product of the evaluations of the one-dimensional ones)
     * @param working iterator working on the GridStorage of the basis
     * @param source array of indices for each dimension (identifying the indices of the current grid point)
     * @param result a vector to store the results in
     */
    void rec(BASIS& basis, std::vector<double>& point, size_t current_dim,
            double value, GridStorage::grid_iterator& working,
            GridStorage::index_type::index_type* source, std::vector<std::pair<
                    size_t, double> >& result) {
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;

        // @todo (blank) Remove 'magic' number
        level_type src_level = static_cast<level_type> (sizeof(index_type) * 8
                - 1);
        index_type src_index = source[current_dim];

        level_type work_level = 1;

        while (true) {
            size_t seq = working.seq();
            if (storage->end(seq)) {
                break;
            } else {
                index_type work_index;
                level_type temp;

                working.get(current_dim, temp, work_index);

                double new_value = basis.eval(work_level, work_index,
                        point[current_dim]);

                if (current_dim == storage->dim() - 1) {
                    result.push_back(std::make_pair(seq, value * new_value));
                } else {
                    rec(basis, point, current_dim + 1, value * new_value,
                            working, source, result);
                }
            }

            if (working.hint()) {
                break;
            }

            // this decides in which direction we should descend by evaluating the corresponding bit
            // the bits are coded from left to right starting with level 1 being in position src_level
            bool right = (src_index & (1 << (src_level - work_level))) > 0;
            ++work_level;

            if (right) {
                working.right_child(current_dim);
            } else {
                working.left_child(current_dim);
            }

        }

        working.top(current_dim);
    }

};


/**
 * Template Specialization for linearboundaryBase basis.
 */
template<>
class GetAffectedBasisFunctions<linearboundaryBase<unsigned int, unsigned int> >
{
	typedef linearboundaryBase<unsigned int, unsigned int> SLinearBoundaryBase;
public:
	GetAffectedBasisFunctions(GridStorage* storage) : storage(storage), BB(storage->getBoundingBox())
	{
	}

	~GetAffectedBasisFunctions() {}

	void operator()(SLinearBoundaryBase& basis, std::vector<double>& point, std::vector<std::pair<size_t, double> >& result)
	{
		bool useBB = false;

		// Check for special bounding box
		if (!this->BB->isTrivialCube())
		{
			useBB = true;
		}

		GridStorage::grid_iterator working(storage);

		working.resetToLevelZero();
		result.clear();
		if (useBB == false)
		{
			rec(basis, point, 0, 1.0, working, result);
		}
		else
		{
			recBB(basis, point, 0, 1.0, working, result);
		}
	}

protected:
	GridStorage* storage;
	BoundingBox* BB;

	void rec(SLinearBoundaryBase& basis, std::vector<double>& point, size_t current_dim, double value, GridStorage::grid_iterator& working, std::vector<std::pair<size_t, double> >& result)
	{
		typedef GridStorage::index_type::level_type level_type;
		typedef GridStorage::index_type::index_type index_type;

		level_type work_level = 0;

		while(true)
		{
			size_t seq = working.seq();
			index_type global_work_index = 0;


			if (storage->end(seq))
			{
				break;
			}
			else
			{
				index_type work_index;
				level_type temp;

				working.get(current_dim, temp, work_index);
				global_work_index = work_index;

				if (work_level > 0)
				{
					double new_value = basis.eval(work_level, work_index, point[current_dim]);

					if(current_dim == storage->dim()-1)
					{
						result.push_back(std::make_pair(seq, value*new_value));
					}
					else
					{
						rec(basis, point, current_dim + 1, value*new_value, working, result);
					}
				}
				// handle boundaries if we are on level 0
				else
				{
					// level 0, index 0
					working.left_levelzero(current_dim);
					size_t seq_lz_left = working.seq();
					double new_value_l_zero_left = basis.eval(0, 0, point[current_dim]);

					if(current_dim == storage->dim()-1)
					{
						result.push_back(std::make_pair(seq_lz_left, value*new_value_l_zero_left));
					}
					else
					{
						rec(basis, point, current_dim + 1, value*new_value_l_zero_left, working, result);
					}

					// level 0, index 1
					working.right_levelzero(current_dim);
					size_t seq_lz_right = working.seq();
					double new_value_l_zero_right = basis.eval(0, 1, point[current_dim]);

					if(current_dim == storage->dim()-1)
					{
						result.push_back(std::make_pair(seq_lz_right, value*new_value_l_zero_right));
					}
					else
					{
						rec(basis, point, current_dim + 1, value*new_value_l_zero_right, working, result);
					}
				}
			}

			// there are no levels left
			if(working.hint())
			{
				break;
			}

			// this decides in which direction we should descend by evaluating the corresponding bit
			// the bits are coded from left to right starting with level 1 being in position src_level
			if (work_level > 0)
			{
				double hat = 0.0;
				level_type h = 0;

				h = 1<<work_level;

				hat = (1.0/static_cast<double>(h))*static_cast<double>(global_work_index);

				if (point[current_dim] == hat)
					break;

				if(point[current_dim] < hat)
				{
					working.left_child(current_dim);
				}
				else
				{
					working.right_child(current_dim);
				}
			}
			else
			{
				if (point[current_dim] == 0.0 || point[current_dim] == 1.0)
					break;

				working.top(current_dim);
			}
			++work_level;
		}

		working.left_levelzero(current_dim);
	}


	void recBB(SLinearBoundaryBase& basis, std::vector<double>& point, size_t current_dim, double value, GridStorage::grid_iterator& working, std::vector<std::pair<size_t, double> >& result)
	{
		typedef GridStorage::index_type::level_type level_type;
		typedef GridStorage::index_type::index_type index_type;

		level_type work_level = 0;

		while(true)
		{
			size_t seq = working.seq();
			index_type global_work_index = 0;


			if (storage->end(seq))
			{
				break;
			}
			else
			{
				index_type work_index;
				level_type temp;

				working.get(current_dim, temp, work_index);
				global_work_index = work_index;

				if (work_level > 0)
				{
					double new_value = basis.eval(work_level, work_index, point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));

					if(current_dim == storage->dim()-1)
					{
						result.push_back(std::make_pair(seq, value*new_value));
					}
					else
					{
						recBB(basis, point, current_dim + 1, value*new_value, working, result);
					}
				}
				// handle boundaries if we are on level 0
				else
				{
					// level 0, index 0
					working.left_levelzero(current_dim);
					size_t seq_lz_left = working.seq();
					double new_value_l_zero_left = basis.eval(0, 0,  point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));

					if(current_dim == storage->dim()-1)
					{
						result.push_back(std::make_pair(seq_lz_left, value*new_value_l_zero_left));
					}
					else
					{
						recBB(basis, point, current_dim + 1, value*new_value_l_zero_left, working, result);
					}

					// level 0, index 1
					working.right_levelzero(current_dim);
					size_t seq_lz_right = working.seq();
					double new_value_l_zero_right = basis.eval(0, 1, point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));

					if(current_dim == storage->dim()-1)
					{
						result.push_back(std::make_pair(seq_lz_right, value*new_value_l_zero_right));
					}
					else
					{
						recBB(basis, point, current_dim + 1, value*new_value_l_zero_right, working, result);
					}
				}
			}

			// there are no levels left
			if(working.hint())
			{
				break;
			}

			// this decides in which direction we should descend by evaluating the corresponding bit
			// the bits are coded from left to right starting with level 1 being in position src_level
			if (work_level > 0)
			{
				double hat = 0.0;
				level_type h = 0;

				h = 1<<work_level;

				hat = (BB->getIntervalWidth(current_dim)*((1.0/static_cast<double>(h))*static_cast<double>(global_work_index))) + BB->getIntervalOffset(current_dim);

				if (point[current_dim] == hat)
					break;

				if(point[current_dim] < hat)
				{
					working.left_child(current_dim);
				}
				else
				{
					working.right_child(current_dim);
				}
			}
			else
			{
				if (point[current_dim] == (BB->getIntervalOffset(current_dim)) || point[current_dim] == ((BB->getIntervalWidth(current_dim))+BB->getIntervalOffset(current_dim)))
					break;

				working.top(current_dim);
			}
			++work_level;
		}

		working.left_levelzero(current_dim);
	}

};

}

#endif /* GETAFFECTEDBASISFUNCTIONS_HPP */
