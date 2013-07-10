/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Dirk Pflüger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef GETAFFECTEDBASISFUNCTIONS_HPP
#define GETAFFECTEDBASISFUNCTIONS_HPP

#include "base/grid/GridStorage.hpp"
#include "base/basis/modwavelet/ModifiedWaveletBasis.hpp"
#include "base/basis/modbspline/ModifiedBsplineBasis.hpp"
#include "base/basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "base/basis/linearstretched/boundary/LinearStretchedBoundaryBasis.hpp"
#include "base/basis/prewavelet/PrewaveletBasis.hpp"

#include <vector>
#include <utility>

namespace sg {
  namespace base {

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
        void operator()(BASIS& basis, std::vector<double>& point, std::vector <
                        std::pair<size_t, double> > & result) {
          GridStorage::grid_iterator working(storage);

          //typedef GridStorage::index_type::level_type level_type;
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
                 GridStorage::index_type::index_type* source, std::vector < std::pair <
                 size_t, double > > & result) {
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
     * Template Specialization for mod_Wavelet basis.
     */
    template<>
    class GetAffectedBasisFunctions<ModifiedWaveletBasis<unsigned int, unsigned int> > {
        typedef ModifiedWaveletBasis<unsigned int, unsigned int> SModWaveletBase;
      public:
        GetAffectedBasisFunctions(GridStorage* storage) :
          storage(storage) {
        }

        ~GetAffectedBasisFunctions() {
        }

        void operator()(SModWaveletBase& basis, std::vector<double>& point,
                        std::vector<std::pair<size_t, double> >& result) {
          GridStorage::grid_iterator working(storage);

          //typedef GridStorage::index_type::level_type level_type;
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

        void rec(SModWaveletBase& basis, std::vector<double>& point,
                 size_t current_dim, double value,
                 GridStorage::grid_iterator& working,
                 GridStorage::index_type::index_type* source, std::vector < std::pair <
                 size_t, double > > & result) {
          typedef GridStorage::index_type::level_type level_type;
          typedef GridStorage::index_type::index_type index_type;
          //size_t i;
          double tmpValue;
          size_t tmpSeq;

          // TODO: Remove 'magic' number
          level_type src_level = static_cast<level_type> (sizeof(index_type) * 8
                                 - 1);
          index_type src_index = source[current_dim];

          level_type work_level = 1;

          while (true) {
            size_t seq = working.seq();

            if (storage->end(seq)) {
              //std::cout << "Grid not found or dim exceed breaking..Grid: " <<seq<<" dim "<<current_dim<<std::endl;
              break;
            } else {
              index_type work_index;
              level_type temp;

              working.get(current_dim, temp, work_index);
              //std::cout <<"current dim: "<<current_dim<<" Point:"<<point[current_dim]<<" current_level: "<<work_level<<" current_Index: " << work_index << std::endl;

              double new_value = basis.eval(work_level, work_index,
                                            point[current_dim]);

              if (current_dim == storage->dim() - 1) {
                //for( int i=0;i < storage->dim();i++)
                //{
                //  working.get(i, temp, work_index);
                //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
                //}
                //std::cout << "Basis function " << seq <<" Value " <<value*new_value<<std::endl;
                result.push_back(std::make_pair(seq, value * new_value));

                //std::cout << "Checking Right Neighbour" << std::endl;
                working.step_right(current_dim);
                tmpSeq = working.seq();

                if (!(storage->end(tmpSeq))) {
                  //for( int i=0;i < storage->dim();i++)
                  //{
                  //  working.get(i, temp, work_index);
                  //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
                  //}
                  working.get(current_dim, temp, work_index);
                  tmpValue = basis.eval(work_level, work_index,
                                        point[current_dim]);
                  //std::cout << "RBasis function " << tmpSeq <<" Value " <<value*tmpValue<<std::endl;
                  result.push_back(std::make_pair(tmpSeq, value
                                                  * tmpValue));

                }

                working.step_left(current_dim);

                //std::cout << "Checking Left Neighbour" << std::endl;
                working.step_left(current_dim);
                tmpSeq = working.seq();

                if (!(storage->end(tmpSeq))) {
                  //for( int i=0;i < storage->dim();i++)
                  //{
                  //  working.get(i, temp, work_index);
                  //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
                  //}
                  working.get(current_dim, temp, work_index);
                  tmpValue = basis.eval(work_level, work_index,
                                        point[current_dim]);
                  //std::cout << "LBasis function " << tmpSeq <<" Value " <<value*tmpValue<<std::endl;
                  result.push_back(std::make_pair(tmpSeq, value
                                                  * tmpValue));

                }

                working.step_right(current_dim);

              } else {
                rec(basis, point, current_dim + 1, value * new_value,
                    working, source, result);

                //std::cout << "Checking Right Neighbour" << std::endl;
                working.step_right(current_dim);
                working.get(current_dim, temp, work_index);
                new_value = basis.eval(work_level, work_index,
                                       point[current_dim]);
                rec(basis, point, current_dim + 1, value * new_value,
                    working, source, result);
                working.step_left(current_dim);

                //std::cout << "Checking left Neighbour"<< std::endl;
                working.step_left(current_dim);
                working.get(current_dim, temp, work_index);
                new_value = basis.eval(work_level, work_index,
                                       point[current_dim]);
                rec(basis, point, current_dim + 1, value * new_value,
                    working, source, result);
                working.step_right(current_dim);
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
     * Template Specialization for B-Spline basis.
     **/
    template<>
    class GetAffectedBasisFunctions<ModifiedBsplineBasis<unsigned int, unsigned int> > {
        typedef ModifiedBsplineBasis<unsigned int, unsigned int> SModBsplineBase;
      public:
        GetAffectedBasisFunctions(GridStorage* storage) : storage(storage) {
        }

        ~GetAffectedBasisFunctions() {}

        void operator()(SModBsplineBase& base, std::vector<double>& point, std::vector<std::pair<size_t, double> >& result) {
          GridStorage::grid_iterator working(storage);


          //typedef GridStorage::index_type::level_type level_type;
          typedef GridStorage::index_type::index_type index_type;

          size_t bits = sizeof(index_type) * 8; // who many levels can we store in a index_type?

          size_t dim = storage->dim();

          index_type* source = new index_type[dim];

          for (size_t d = 0; d < dim; ++d) {
            // This does not really work on grids with borders.
            double temp = floor(point[d] * (1 << (bits - 2))) * 2;

            if (point[d] == 1.0) {
              source[d] = static_cast<index_type>(temp - 1);
            } else {
              source[d] = static_cast<index_type>(temp + 1);
            }

          }

          result.clear();
          rec(base, point, 0, 1.0, working, source, result);

          delete [] source;

        }

      protected:
        GridStorage* storage;

        void rec(SModBsplineBase& base, std::vector<double>& point, size_t current_dim, double value, GridStorage::grid_iterator& working, GridStorage::index_type::index_type* source, std::vector<std::pair<size_t, double> >& result) {
          typedef GridStorage::index_type::level_type level_type;
          typedef GridStorage::index_type::index_type index_type;
          //size_t i;
          double tmpValue;
          size_t tmpSeq;

          // TODO: Remove 'magic' number
          level_type src_level = static_cast<level_type>(sizeof(index_type) * 8 - 1);
          index_type src_index = source[current_dim];

          level_type work_level = 1;

          while (true) {
            size_t seq = working.seq();

            if ( storage->end(seq)  ) {
              //std::cout << "Grid not found or dim exceed breaking..Grid: " <<seq<<" dim "<<current_dim<<std::endl;
              break;
            } else {
              index_type work_index;
              level_type temp;

              working.get(current_dim, temp, work_index);
              //std::cout <<"current dim: "<<current_dim<<" Point:"<<point[current_dim]<<" current_level: "<<work_level<<" current_Index: " << work_index << std::endl;

              double new_value = base.eval(work_level, work_index, point[current_dim]);

              if (current_dim == storage->dim() - 1) {
                //for( int i=0;i < storage->dim();i++)
                //{
                //  working.get(i, temp, work_index);
                //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
                //}
                //std::cout << "Basis function " << seq <<" Value " <<value*new_value<<std::endl;
                result.push_back(std::make_pair(seq, value * new_value));

                //std::cout << "Checking Right Neighbour" << std::endl;
                working.step_right(current_dim);
                tmpSeq = working.seq();

                if ( !(storage->end(tmpSeq)) ) {
                  //for( int i=0;i < storage->dim();i++)
                  //{
                  //  working.get(i, temp, work_index);
                  //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
                  //}
                  working.get(current_dim, temp, work_index);
                  tmpValue = base.eval(work_level, work_index, point[current_dim]);
                  //std::cout << "RBasis function " << tmpSeq <<" Value " <<value*tmpValue<<std::endl;
                  result.push_back(std::make_pair(tmpSeq, value * tmpValue));

                }

                working.step_left(current_dim);

                //std::cout << "Checking Left Neighbour" << std::endl;
                working.step_left(current_dim);
                tmpSeq = working.seq();

                if (!(storage->end(tmpSeq)) ) {
                  //for( int i=0;i < storage->dim();i++)
                  //{
                  //  working.get(i, temp, work_index);
                  //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
                  //}
                  working.get(current_dim, temp, work_index);
                  tmpValue = base.eval(work_level, work_index, point[current_dim]);
                  //std::cout << "LBasis function " << tmpSeq <<" Value " <<value*tmpValue<<std::endl;
                  result.push_back(std::make_pair(tmpSeq, value * tmpValue));

                }

                working.step_right(current_dim);

              } else {
                rec(base, point, current_dim + 1, value * new_value, working, source, result);

                //std::cout << "Checking Right Neighbour" << std::endl;
                working.step_right(current_dim);
                working.get(current_dim, temp, work_index);
                new_value = base.eval(work_level, work_index, point[current_dim]);
                rec(base, point, current_dim + 1, value * new_value, working, source, result);
                working.step_left(current_dim);

                //std::cout << "Checking left Neighbour"<< std::endl;
                working.step_left(current_dim);
                working.get(current_dim, temp, work_index);
                new_value = base.eval(work_level, work_index, point[current_dim]);
                rec(base, point, current_dim + 1 , value * new_value, working, source, result);
                working.step_right(current_dim);
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
     * Template Specialization for LinearBoundaryBasis basis.
     */
    template<>
    class GetAffectedBasisFunctions<LinearBoundaryBasis<unsigned int, unsigned int> > {
        typedef LinearBoundaryBasis<unsigned int, unsigned int> SLinearBoundaryBase;
      public:
        GetAffectedBasisFunctions(GridStorage* storage) : storage(storage), BB(storage->getBoundingBox()) {
        }

        ~GetAffectedBasisFunctions() {}

        void operator()(SLinearBoundaryBase& basis, std::vector<double>& point, std::vector<std::pair<size_t, double> >& result) {
          bool useBB = false;

          // Check for special bounding box
          if (!this->BB->isTrivialCube()) {
            useBB = true;
          }

          GridStorage::grid_iterator working(storage);

          working.resetToLevelZero();
          result.clear();

          if (useBB == false) {
            rec(basis, point, 0, 1.0, working, result);
          } else {
            recBB(basis, point, 0, 1.0, working, result);
          }
        }

      protected:
        GridStorage* storage;
        BoundingBox* BB;

        void rec(SLinearBoundaryBase& basis, std::vector<double>& point, size_t current_dim, double value, GridStorage::grid_iterator& working, std::vector<std::pair<size_t, double> >& result) {
          typedef GridStorage::index_type::level_type level_type;
          typedef GridStorage::index_type::index_type index_type;

          level_type work_level = 0;

          while (true) {
            size_t seq = working.seq();
            index_type global_work_index = 0;


            if (storage->end(seq)) {
              break;
            } else {
              index_type work_index;
              level_type temp;

              working.get(current_dim, temp, work_index);
              global_work_index = work_index;

              if (work_level > 0) {
                double new_value = basis.eval(work_level, work_index, point[current_dim]);

                if (current_dim == storage->dim() - 1) {
                  result.push_back(std::make_pair(seq, value * new_value));
                } else {
                  rec(basis, point, current_dim + 1, value * new_value, working, result);
                }
              }
              // handle boundaries if we are on level 0
              else {
                // level 0, index 0
                working.left_levelzero(current_dim);
                size_t seq_lz_left = working.seq();
                double new_value_l_zero_left = basis.eval(0, 0, point[current_dim]);

                if (current_dim == storage->dim() - 1) {
                  result.push_back(std::make_pair(seq_lz_left, value * new_value_l_zero_left));
                } else {
                  rec(basis, point, current_dim + 1, value * new_value_l_zero_left, working, result);
                }

                // level 0, index 1
                working.right_levelzero(current_dim);
                size_t seq_lz_right = working.seq();
                double new_value_l_zero_right = basis.eval(0, 1, point[current_dim]);

                if (current_dim == storage->dim() - 1) {
                  result.push_back(std::make_pair(seq_lz_right, value * new_value_l_zero_right));
                } else {
                  rec(basis, point, current_dim + 1, value * new_value_l_zero_right, working, result);
                }
              }
            }

            // there are no levels left
            if (working.hint()) {
              break;
            }

            // this decides in which direction we should descend by evaluating the corresponding bit
            // the bits are coded from left to right starting with level 1 being in position src_level
            if (work_level > 0) {
              double hat = 0.0;
              level_type h = 0;

              h = 1 << work_level;

              hat = (1.0 / static_cast<double>(h)) * static_cast<double>(global_work_index);

              if (point[current_dim] == hat)
                break;

              if (point[current_dim] < hat) {
                working.left_child(current_dim);
              } else {
                working.right_child(current_dim);
              }
            } else {
              if (point[current_dim] == 0.0 || point[current_dim] == 1.0)
                break;

              working.top(current_dim);
            }

            ++work_level;
          }

          working.left_levelzero(current_dim);
        }


        void recBB(SLinearBoundaryBase& basis, std::vector<double>& point, size_t current_dim, double value, GridStorage::grid_iterator& working, std::vector<std::pair<size_t, double> >& result) {
          typedef GridStorage::index_type::level_type level_type;
          typedef GridStorage::index_type::index_type index_type;

          level_type work_level = 0;

          while (true) {
            size_t seq = working.seq();
            index_type global_work_index = 0;


            if (storage->end(seq)) {
              break;
            } else {
              index_type work_index;
              level_type temp;

              working.get(current_dim, temp, work_index);
              global_work_index = work_index;

              if (work_level > 0) {
                double new_value = basis.eval(work_level, work_index, point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));

                if (current_dim == storage->dim() - 1) {
                  result.push_back(std::make_pair(seq, value * new_value));
                } else {
                  recBB(basis, point, current_dim + 1, value * new_value, working, result);
                }
              }
              // handle boundaries if we are on level 0
              else {
                // level 0, index 0
                working.left_levelzero(current_dim);
                size_t seq_lz_left = working.seq();
                double new_value_l_zero_left = basis.eval(0, 0,  point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));

                if (current_dim == storage->dim() - 1) {
                  result.push_back(std::make_pair(seq_lz_left, value * new_value_l_zero_left));
                } else {
                  recBB(basis, point, current_dim + 1, value * new_value_l_zero_left, working, result);
                }

                // level 0, index 1
                working.right_levelzero(current_dim);
                size_t seq_lz_right = working.seq();
                double new_value_l_zero_right = basis.eval(0, 1, point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));

                if (current_dim == storage->dim() - 1) {
                  result.push_back(std::make_pair(seq_lz_right, value * new_value_l_zero_right));
                } else {
                  recBB(basis, point, current_dim + 1, value * new_value_l_zero_right, working, result);
                }
              }
            }

            // there are no levels left
            if (working.hint()) {
              break;
            }

            // this decides in which direction we should descend by evaluating the corresponding bit
            // the bits are coded from left to right starting with level 1 being in position src_level
            if (work_level > 0) {
              double hat = 0.0;
              level_type h = 0;

              h = 1 << work_level;

              hat = (BB->getIntervalWidth(current_dim) * ((1.0 / static_cast<double>(h)) * static_cast<double>(global_work_index))) + BB->getIntervalOffset(current_dim);

              if (point[current_dim] == hat)
                break;

              if (point[current_dim] < hat) {
                working.left_child(current_dim);
              } else {
                working.right_child(current_dim);
              }
            } else {
              if (point[current_dim] == (BB->getIntervalOffset(current_dim)) || point[current_dim] == ((BB->getIntervalWidth(current_dim)) + BB->getIntervalOffset(current_dim)))
                break;

              working.top(current_dim);
            }

            ++work_level;
          }

          working.left_levelzero(current_dim);
        }

    };

    /**
     * Template Specialization for LinearStretchedBoundaryBasis basis.
     */
    template<>
    class GetAffectedBasisFunctions<LinearStretchedBoundaryBasis<unsigned int, unsigned int> > {
        typedef LinearStretchedBoundaryBasis<unsigned int, unsigned int> SLinearStretchedBoundaryBase;
      public:
        GetAffectedBasisFunctions(GridStorage* storage) : storage(storage), stretch(storage->getStretching()) {
        }

        ~GetAffectedBasisFunctions() {}

        void operator()(SLinearStretchedBoundaryBase& basis, std::vector<double>& point, std::vector<std::pair<size_t, double> >& result) {
          //There's no BoundaryBox checking necessary

          GridStorage::grid_iterator working(storage);

          working.resetToLevelZero();
          result.clear();

          recBB(basis, point, 0, 1.0, working, result);


        }

      protected:
        GridStorage* storage;
        Stretching* stretch;

        void recBB(SLinearStretchedBoundaryBase& basis, std::vector<double>& point, size_t current_dim, double value, GridStorage::grid_iterator& working, std::vector<std::pair<size_t, double> >& result) {
          typedef GridStorage::index_type::level_type level_type;
          typedef GridStorage::index_type::index_type index_type;
          double posl = 0, posr = 0, posc = 0;

          level_type work_level = 0;

          while (true) {
            size_t seq = working.seq();
            index_type global_work_index = 0;


            if (storage->end(seq)) {
              break;
            } else {
              index_type work_index;
              level_type temp;

              working.get(current_dim, temp, work_index);
              global_work_index = work_index;




              if (work_level > 0) {
                //          stretch->getAdjacentPositions(static_cast<int>(temp), static_cast<int>(work_index), current_dim,posl,posr);
                //          posc=stretch->getCoordinates(static_cast<int>(temp), static_cast<int>(work_index), current_dim);
                stretch->getAdjacentPositions(static_cast<int>(temp), static_cast<int>(work_index), current_dim, posc, posl, posr);
                //double new_value = basis.eval(work_level, work_index, point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));
                double new_value;

                if (posc < point[current_dim]) {
                  /// calculate using the right-hand ramp
                  new_value = basis.eval(work_level, work_index, point[current_dim], posr, posc);
                } else {
                  ///calculate using the left-hand ramp
                  new_value = basis.eval(work_level, work_index, point[current_dim], posl, posc);
                }


                if (current_dim == storage->dim() - 1) {
                  result.push_back(std::make_pair(seq, value * new_value));
                } else {
                  recBB(basis, point, current_dim + 1, value * new_value, working, result);
                }
              }
              // handle boundaries if we are on level 0
              else {
                double left = stretch->getBoundary(current_dim).leftBoundary;
                double right = stretch->getBoundary(current_dim).rightBoundary;
                // level 0, index 0
                working.left_levelzero(current_dim);
                size_t seq_lz_left = working.seq();
                //double new_value_l_zero_left = basis.eval(0, 0,  point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));
                double new_value_l_zero_left = basis.eval(0, 0,  point[current_dim], right, left);

                if (current_dim == storage->dim() - 1) {
                  result.push_back(std::make_pair(seq_lz_left, value * new_value_l_zero_left));
                } else {
                  recBB(basis, point, current_dim + 1, value * new_value_l_zero_left, working, result);
                }

                // level 0, index 1
                working.right_levelzero(current_dim);
                size_t seq_lz_right = working.seq();
                //double new_value_l_zero_right = basis.eval(0, 1, point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));
                double new_value_l_zero_right = basis.eval(0, 1,  point[current_dim], left, right);

                if (current_dim == storage->dim() - 1) {
                  result.push_back(std::make_pair(seq_lz_right, value * new_value_l_zero_right));
                } else {
                  recBB(basis, point, current_dim + 1, value * new_value_l_zero_right, working, result);
                }
              }
            }

            // there are no levels left
            if (working.hint()) {
              break;
            }

            // this decides in which direction we should descend by evaluating the corresponding bit
            // the bits are coded from left to right starting with level 1 being in position src_level
            if (work_level > 0) {
              double hat = 0.0;
              //        level_type h = 0;
              //        h = 1<<work_level;
              //        hat = (stretch->getIntervalWidth(current_dim)*((1.0/static_cast<double>(h))*static_cast<double>(global_work_index))) + stretch->getIntervalOffset(current_dim);
              hat = stretch->getCoordinates(static_cast<int>(work_level), static_cast<int>(global_work_index), current_dim);


              if (point[current_dim] == hat)
                break;

              if (point[current_dim] < hat) {
                working.left_child(current_dim);
              } else {
                working.right_child(current_dim);
              }
            } else {
              //if (point[current_dim] == (BB->getIntervalOffset(current_dim)) || point[current_dim] == ((BB->getIntervalWidth(current_dim))+BB->getIntervalOffset(current_dim)))
              //break;

              if ((point[current_dim] == stretch->getCoordinates(0, 0, current_dim)) || (point[current_dim] == stretch->getCoordinates(0, 1, current_dim))) {
                break;
              }

              working.top(current_dim);
            }

            ++work_level;
          }

          working.left_levelzero(current_dim);
        }

    };



    /**
     * Template Specialization for prewavelet basis.
     */
    template<>
    class GetAffectedBasisFunctions<PrewaveletBasis<unsigned int, unsigned int> > {

        typedef PrewaveletBasis<unsigned int, unsigned int> SPrewaveletBase;
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
        void operator()(SPrewaveletBase& basis, std::vector<double>& point,
                        std::vector<std::pair<size_t, double> >& result) {

          GridStorage::grid_iterator iter(storage);
          result.clear();
          rec(basis, point, 0, iter, result);


        }


      protected:
        GridStorage* storage;
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;

        /**
         * Recursive traversal of the "tree" of basis functions for evaluation, used in operator().
         * For a given evaluation point \f$x\f$, it stores tuples (std::pair) of
         * \f$(i,\phi_i(x))\f$ in the result vector for all basis functions that are non-zero.
         *
         * @param basis a sparse grid basis
         * @param point evaluation point within the domain
         * @param current_dim the dimension currently looked at (recursion parameter)
         * @param iter iterator working on the GridStorage of the basis
         * @param result a vector to store the results in
         */

        void rec(SPrewaveletBase& basis, std::vector<double>& point,
                 size_t current_dim, GridStorage::grid_iterator& iter, std::vector <
                 std::pair<size_t, double> > & result) {

          //First, save this point
          double value = 1.0;

          for (size_t d = 0; d < storage->dim(); ++d) {
            index_type current_index;
            level_type current_level;
            iter.get(d, current_level, current_index);
            value *= basis.eval(current_level, current_index, point[d]);
          }

          result.push_back(std::make_pair(iter.seq(), value));

          for (size_t d = current_dim; d < storage->dim(); d++) {
            index_type save_index;
            level_type save_level;
            iter.get(d, save_level, save_index); //Save current index

            if (iter.hint_left(d)) {
              //Handle left Child in dimension d
              iter.left_child(d);
              index_type current_index;
              level_type current_level;
              iter.get(d, current_level, current_index);
              //Check if current index is still in the area of the point
              int int_index = current_index; // needed to avoid overrun with index-3

              if ((1.0 / (1 << current_level)) * (int_index - 3) < point[d]
                  && (1.0 / (1 << current_level)) * (int_index + 3)
                  > point[d]) {
                rec(basis, point, d, iter, result);
              }
            }

            iter.set(d, save_level, save_index); //reset index

            if (iter.hint_right(d)) {
              //Handle left Child in dimension d
              iter.right_child(d);
              index_type current_index;
              level_type current_level;
              iter.get(d, current_level, current_index);
              //Check if current index is still in the area of the point
              int int_index = current_index; // needed to avoid overrun with index-3

              if ((1.0 / (1 << current_level)) * (int_index - 3) < point[d]
                  && (1.0 / (1 << current_level)) * (int_index + 3)
                  > point[d]) {
                rec(basis, point, d, iter, result);
              }
            }

            iter.set(d, save_level, save_index); //reset index
          }

        }

    };
  }
}

#endif /* GETAFFECTEDBASISFUNCTIONS_HPP */
