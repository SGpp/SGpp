// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef ALGORITHMEVALUATIONTRANSPOSED_HPP
#define ALGORITHMEVALUATIONTRANSPOSED_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/basis/waveletModified/WaveletModifiedBasis.hpp>
#include <sgpp/base/basis/bsplineModified/BsplineModifiedBasis.hpp>
#include <sgpp/base/basis/linearBoundary/LinearBoundaryBasis.hpp>

#include <vector>
#include <utility>

#include <sgpp/globaldef.hpp>


namespace SGPP
{
namespace base
{

/**
 * Basic algorithm for getting all affected basis functions.
 * This implicitly assumes a tensor-product approach and local support.
 * No grid points on the border are supported.
 */
template<class BASIS>
class AlgorithmEvaluationTransposed
{
public:
    AlgorithmEvaluationTransposed(GridStorage* storage) :
            storage(storage)
    {}

    ~AlgorithmEvaluationTransposed()
    {}

    /**
     * Returns evaluations of all basis functions that are non-zero at a given evaluation point.
     * For a given evaluation point \f$x\f$, it stores tuples (std::pair) of
     * \f$(i,\phi_i(x))\f$ in the result vector for all basis functions that are non-zero.
     * If one wants to evaluate \f$f_N(x)\f$, one only has to compute
     * \f[ \sum_{r\in\mathbf{result}} \alpha[r\rightarrow\mathbf{first}] \cdot r\rightarrow\mathbf{second}. \f]
     *
     * @param basis a sparse grid basis
     * @param point evaluation point within the domain
     * @param alpha the coefficient of the regarded ansatzfunction
     * @param result vector that will contain the local support of the given ansatzfuction for all evaluations points
     */
    void operator()(BASIS& basis, std::vector<double>& point, double alpha, DataVector& result)
    {
        GridStorage::grid_iterator working(storage);

        //typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;

        size_t bits = sizeof(index_type) * 8; // how many levels can we store in a index_type?

        size_t dim = storage->dim();

        // Check for bounding box
        BoundingBox *bb = storage->getBoundingBox();
        if ( bb != NULL )
        {
            bool inside = true;
            for (size_t d = 0; d < dim; ++d)
            {
                DimensionBoundary dimbb = bb->getBoundary(d);

                if ( ! (dimbb.leftBoundary <= point[d] && point[d] <= dimbb.rightBoundary) )
                {
                    inside = false;
                    break;
                }
            }

            if ( !inside )
            {
            	//nothing to change
                return;
            }
            else
            {
                for (size_t d = 0; d < dim; ++d)
                {
                    DimensionBoundary dimbb = bb->getBoundary(d);

                    point[d] = (point[d] - dimbb.leftBoundary) / (dimbb.rightBoundary - dimbb.leftBoundary);
                }
            }
        }

        index_type* source = new index_type[dim];

        for (size_t d = 0; d < dim; ++d)
        {
            // This does not really work on grids with borders.
            double temp = floor(point[d] * (1 << (bits - 2))) * 2;

            if (point[d] == 1.0)
            {
                source[d] = static_cast<index_type> (temp - 1);
            }
            else
            {
                source[d] = static_cast<index_type> (temp + 1);
            }
        }

        rec(basis, point, 0, 1.0, working, source, alpha, result);
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
     * @param alpha the coefficient of current ansatzfunction
     * @param result vector that will contain the local support of the given ansatzfuction for all evaluations points
     */
    void rec(BASIS& basis, std::vector<double>& point, size_t current_dim,
             double value, GridStorage::grid_iterator& working,
             GridStorage::index_type::index_type* source, double alpha, DataVector& result)
    {
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;

        // @todo (blank) Remove 'magic' number
        level_type src_level = static_cast<level_type> (sizeof(index_type) * 8 - 1);
        index_type src_index = source[current_dim];

        level_type work_level = 1;

        while (true)
        {
            size_t seq = working.seq();

            if (storage->end(seq))
            {
                break;
            }
            else
            {
                index_type work_index;
                level_type temp;

                working.get(current_dim, temp, work_index);

                double new_value = basis.eval(work_level, work_index,
                                              point[current_dim]);
                new_value *= value;

                if (current_dim == storage->dim() - 1)
                {
                    result[seq] += (alpha * new_value);
                }
                else
                {
                    rec(basis, point, current_dim + 1, new_value,
                        working, source, alpha, result);
                }
            }

            if (working.hint())
            {
                break;
            }

            // this decides in which direction we should descend by evaluating the corresponding bit
            // the bits are coded from left to right starting with level 1 being in position src_level
            bool right = (src_index & (1 << (src_level - work_level))) > 0;
            ++work_level;

            if (right)
            {
                working.right_child(current_dim);
            }
            else
            {
                working.left_child(current_dim);
            }

        }

        working.top(current_dim);
    }

};

///**
// * Template Specialization for mod_Wavelet basis.
// */
//template<>
//class GetAffectedBasisFunctions<WaveletModifiedBasis<unsigned int, unsigned int> > {
//    typedef WaveletModifiedBasis<unsigned int, unsigned int> SWaveletModifiedBase;
//public:
//    GetAffectedBasisFunctions(GridStorage* storage) :
//        storage(storage) {
//    }
//
//    ~GetAffectedBasisFunctions() {
//    }
//
//    void operator()(SWaveletModifiedBase& basis, std::vector<double>& point,
//            std::vector<std::pair<size_t, double> >& result) {
//        GridStorage::grid_iterator working(storage);
//
//        typedef GridStorage::index_type::level_type level_type;
//        typedef GridStorage::index_type::index_type index_type;
//
//        size_t bits = sizeof(index_type) * 8; // how many levels can we store in a index_type?
//
//        size_t dim = storage->dim();
//
//        index_type* source = new index_type[dim];
//
//        for (size_t d = 0; d < dim; ++d) {
//            // This does not really work on grids with borders.
//            double temp = floor(point[d] * (1 << (bits - 2))) * 2;
//            if (point[d] == 1.0) {
//                source[d] = static_cast<index_type> (temp - 1);
//            } else {
//                source[d] = static_cast<index_type> (temp + 1);
//            }
//
//        }
//
//        result.clear();
//        rec(basis, point, 0, 1.0, working, source, result);
//
//        delete[] source;
//
//    }
//
//protected:
//    GridStorage* storage;
//
//    void rec(SWaveletModifiedBase& basis, std::vector<double>& point,
//            size_t current_dim, double value,
//            GridStorage::grid_iterator& working,
//            GridStorage::index_type::index_type* source, std::vector<std::pair<
//                    size_t, double> >& result) {
//        typedef GridStorage::index_type::level_type level_type;
//        typedef GridStorage::index_type::index_type index_type;
//        //size_t i;
//        double tmpValue;
//        size_t tmpSeq;
//
//        // TODO: Remove 'magic' number
//        level_type src_level = static_cast<level_type> (sizeof(index_type) * 8
//                - 1);
//        index_type src_index = source[current_dim];
//
//        level_type work_level = 1;
//
//        while (true) {
//            size_t seq = working.seq();
//            if (storage->end(seq)) {
//                //std::cout << "Grid not found or dim exceed breaking..Grid: " <<seq<<" dim "<<current_dim<<std::endl;
//                break;
//            } else {
//                index_type work_index;
//                level_type temp;
//
//                working.get(current_dim, temp, work_index);
//                //std::cout <<"current dim: "<<current_dim<<" Point:"<<point[current_dim]<<" current_level: "<<work_level<<" current_Index: " << work_index << std::endl;
//
//                double new_value = basis.eval(work_level, work_index,
//                        point[current_dim]);
//
//                if (current_dim == storage->dim() - 1) {
//                    //for( int i=0;i < storage->dim();i++)
//                    //{
//                    //  working.get(i, temp, work_index);
//                    //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
//                    //}
//                    //std::cout << "Basis function " << seq <<" Value " <<value*new_value<<std::endl;
//                    result.push_back(std::make_pair(seq, value * new_value));
//
//                    //std::cout << "Checking Right Neighbour" << std::endl;
//                    working.step_right(current_dim);
//                    tmpSeq = working.seq();
//                    if (!(storage->end(tmpSeq))) {
//                        //for( int i=0;i < storage->dim();i++)
//                        //{
//                        //  working.get(i, temp, work_index);
//                        //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
//                        //}
//                        working.get(current_dim, temp, work_index);
//                        tmpValue = basis.eval(work_level, work_index,
//                                point[current_dim]);
//                        //std::cout << "RBasis function " << tmpSeq <<" Value " <<value*tmpValue<<std::endl;
//                        result.push_back(std::make_pair(tmpSeq, value
//                                * tmpValue));
//
//                    }
//                    working.step_left(current_dim);
//
//                    //std::cout << "Checking Left Neighbour" << std::endl;
//                    working.step_left(current_dim);
//                    tmpSeq = working.seq();
//                    if (!(storage->end(tmpSeq))) {
//                        //for( int i=0;i < storage->dim();i++)
//                        //{
//                        //  working.get(i, temp, work_index);
//                        //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
//                        //}
//                        working.get(current_dim, temp, work_index);
//                        tmpValue = basis.eval(work_level, work_index,
//                                point[current_dim]);
//                        //std::cout << "LBasis function " << tmpSeq <<" Value " <<value*tmpValue<<std::endl;
//                        result.push_back(std::make_pair(tmpSeq, value
//                                * tmpValue));
//
//                    }
//                    working.step_right(current_dim);
//
//                } else {
//                    rec(basis, point, current_dim + 1, value * new_value,
//                            working, source, result);
//
//                    //std::cout << "Checking Right Neighbour" << std::endl;
//                    working.step_right(current_dim);
//                    working.get(current_dim, temp, work_index);
//                    new_value = basis.eval(work_level, work_index,
//                            point[current_dim]);
//                    rec(basis, point, current_dim + 1, value * new_value,
//                            working, source, result);
//                    working.step_left(current_dim);
//
//                    //std::cout << "Checking left Neighbour"<< std::endl;
//                    working.step_left(current_dim);
//                    working.get(current_dim, temp, work_index);
//                    new_value = basis.eval(work_level, work_index,
//                            point[current_dim]);
//                    rec(basis, point, current_dim + 1, value * new_value,
//                            working, source, result);
//                    working.step_right(current_dim);
//                }
//
//            }
//
//            if (working.hint()) {
//                break;
//            }
//
//            // this decides in which direction we should descend by evaluating the corresponding bit
//            // the bits are coded from left to right starting with level 1 being in position src_level
//            bool right = (src_index & (1 << (src_level - work_level))) > 0;
//            ++work_level;
//
//            if (right) {
//                working.right_child(current_dim);
//            } else {
//                working.left_child(current_dim);
//            }
//
//        }
//
//        working.top(current_dim);
//    }
//
//};
//
//
///**
// * Template Specialization for B-Spline basis.
// **/
//template<>
//class GetAffectedBasisFunctions<BsplineModifiedBasis<unsigned int, unsigned int> >
//{
//    typedef BsplineModifiedBasis<unsigned int, unsigned int> SBsplineModifiedBase;
//public:
//    GetAffectedBasisFunctions(GridStorage* storage) : storage(storage)
//    {
//    }
//
//    ~GetAffectedBasisFunctions() {}
//
//    void operator()(SBsplineModifiedBase& base, std::vector<double>& point, std::vector<std::pair<size_t, double> >& result)
//    {
//        GridStorage::grid_iterator working(storage);
//
//
//        typedef GridStorage::index_type::level_type level_type;
//        typedef GridStorage::index_type::index_type index_type;
//
//        size_t bits = sizeof(index_type) * 8; // who many levels can we store in a index_type?
//
//        size_t dim = storage->dim();
//
//        index_type* source = new index_type[dim];
//
//        for(size_t d = 0; d < dim; ++d)
//        {
//            // This does not really work on grids with borders.
//            double temp = floor(point[d]*(1<<(bits-2)))*2;
//            if(point[d] == 1.0)
//            {
//                source[d] = static_cast<index_type>(temp-1);
//            }
//            else
//            {
//                source[d] = static_cast<index_type>(temp+1);
//            }
//
//        }
//
//        result.clear();
//        rec(base, point, 0, 1.0, working, source, result);
//
//        delete [] source;
//
//    }
//
//protected:
//    GridStorage* storage;
//
//    /**
//     * Example implementation of storage agnostic algorithm.
//     * Returns all affected base functions
//     */
//    void rec(SBsplineModifiedBase& base, std::vector<double>& point, size_t current_dim, double value, GridStorage::grid_iterator& working, GridStorage::index_type::index_type* source, std::vector<std::pair<size_t, double> >& result)
//    {
//        typedef GridStorage::index_type::level_type level_type;
//        typedef GridStorage::index_type::index_type index_type;
//        //size_t i;
//        double tmpValue;
//        size_t tmpSeq;
//
//        // TODO: Remove 'magic' number
//        level_type src_level = static_cast<level_type>(sizeof(index_type) * 8 - 1);
//        index_type src_index = source[current_dim];
//
//        level_type work_level = 1;
//
//        while(true)
//        {
//            size_t seq = working.seq();
//            if( storage->end(seq)  )
//            {
//                //std::cout << "Grid not found or dim exceed breaking..Grid: " <<seq<<" dim "<<current_dim<<std::endl;
//                break;
//            }
//            else
//            {
//                index_type work_index;
//                level_type temp;
//
//                working.get(current_dim, temp, work_index);
//                //std::cout <<"current dim: "<<current_dim<<" Point:"<<point[current_dim]<<" current_level: "<<work_level<<" current_Index: " << work_index << std::endl;
//
//                double new_value = base.eval(work_level, work_index, point[current_dim]);
//
//                if(current_dim == storage->dim()-1)
//                {
//                    //for( int i=0;i < storage->dim();i++)
//                    //{
//                    //  working.get(i, temp, work_index);
//                    //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
//                    //}
//                    //std::cout << "Basis function " << seq <<" Value " <<value*new_value<<std::endl;
//                    result.push_back(std::make_pair(seq, value*new_value));
//
//                    //std::cout << "Checking Right Neighbour" << std::endl;
//                    working.step_right(current_dim);
//                    tmpSeq = working.seq();
//                    if( !(storage->end(tmpSeq)) )
//                    {
//                        //for( int i=0;i < storage->dim();i++)
//                        //{
//                        //  working.get(i, temp, work_index);
//                        //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
//                        //}
//                        working.get(current_dim, temp, work_index);
//                        tmpValue=base.eval(work_level, work_index, point[current_dim]);
//                        //std::cout << "RBasis function " << tmpSeq <<" Value " <<value*tmpValue<<std::endl;
//                        result.push_back(std::make_pair(tmpSeq, value*tmpValue));
//
//                    }
//                    working.step_left(current_dim);
//
//                    //std::cout << "Checking Left Neighbour" << std::endl;
//                    working.step_left(current_dim);
//                    tmpSeq = working.seq();
//                    if(!(storage->end(tmpSeq)) )
//                    {
//                        //for( int i=0;i < storage->dim();i++)
//                        //{
//                        //  working.get(i, temp, work_index);
//                        //  std::cout <<" dim "<<i <<" level "<<temp<<" Index "<<work_index<<std::endl;
//                        //}
//                        working.get(current_dim, temp, work_index);
//                        tmpValue=base.eval(work_level, work_index, point[current_dim]);
//                        //std::cout << "LBasis function " << tmpSeq <<" Value " <<value*tmpValue<<std::endl;
//                        result.push_back(std::make_pair(tmpSeq, value*tmpValue));
//
//                    }
//                    working.step_right(current_dim);
//
//                }
//                else
//                {
//                    rec(base, point, current_dim + 1, value*new_value, working, source, result);
//
//                    //std::cout << "Checking Right Neighbour" << std::endl;
//                    working.step_right(current_dim);
//                    working.get(current_dim, temp, work_index);
//                    new_value = base.eval(work_level, work_index, point[current_dim]);
//                    rec(base, point, current_dim + 1, value*new_value, working, source, result);
//                    working.step_left(current_dim);
//
//                    //std::cout << "Checking left Neighbour"<< std::endl;
//                    working.step_left(current_dim);
//                    working.get(current_dim, temp, work_index);
//                    new_value = base.eval(work_level, work_index, point[current_dim]);
//                    rec(base, point, current_dim +1 , value*new_value, working, source, result);
//                    working.step_right(current_dim);
//                }
//
//            }
//
//            if(working.hint())
//            {
//                break;
//            }
//
//            // this decides in which direction we should descend by evaluating the corresponding bit
//            // the bits are coded from left to right starting with level 1 being in position src_level
//            bool right = (src_index & (1 << (src_level - work_level))) > 0;
//            ++work_level;
//
//            if(right)
//            {
//                working.right_child(current_dim);
//            }
//            else
//            {
//                working.left_child(current_dim);
//            }
//
//        }
//
//        working.top(current_dim);
//    }
//
//};
//
//template<>
//class GetAffectedBasisFunctions<LinearBoundaryBasis<unsigned int, unsigned int> >
//{
//  typedef LinearBoundaryBasis<unsigned int, unsigned int> SLinearBoundaryBase;
//public:
//  GetAffectedBasisFunctions(GridStorage* storage) : storage(storage), BB(storage->getBoundingBox())
//  {
//  }
//
//  ~GetAffectedBasisFunctions() {}
//
//  void operator()(SLinearBoundaryBase& basis, std::vector<double>& point, std::vector<std::pair<size_t, double> >& result)
//  {
//    bool useBB = false;
//
//    // Check for special bounding box
//    if (!this->BB->isTrivialCube())
//    {
//      useBB = true;
//    }
//
//    GridStorage::grid_iterator working(storage);
//
//    working.resetToLevelZero();
//    result.clear();
//    if (useBB == false)
//    {
//      rec(basis, point, 0, 1.0, working, result);
//    }
//    else
//    {
//      recBB(basis, point, 0, 1.0, working, result);
//    }
//  }
//
//protected:
//  GridStorage* storage;
//  BoundingBox* BB;
//
//  void rec(SLinearBoundaryBase& basis, std::vector<double>& point, size_t current_dim, double value, GridStorage::grid_iterator& working, std::vector<std::pair<size_t, double> >& result)
//  {
//    typedef GridStorage::index_type::level_type level_type;
//    typedef GridStorage::index_type::index_type index_type;
//
//    level_type work_level = 0;
//
//    while(true)
//    {
//      size_t seq = working.seq();
//      index_type global_work_index = 0;
//
//
//      if (storage->end(seq))
//      {
//        break;
//      }
//      else
//      {
//        index_type work_index;
//        level_type temp;
//
//        working.get(current_dim, temp, work_index);
//        global_work_index = work_index;
//
//        if (work_level > 0)
//        {
//          double new_value = basis.eval(work_level, work_index, point[current_dim]);
//
//          if(current_dim == storage->dim()-1)
//          {
//            result.push_back(std::make_pair(seq, value*new_value));
//          }
//          else
//          {
//            rec(basis, point, current_dim + 1, value*new_value, working, result);
//          }
//        }
//        // handle boundaries if we are on level 0
//        else
//        {
//          // level 0, index 0
//          working.left_levelzero(current_dim);
//          size_t seq_lz_left = working.seq();
//          double new_value_l_zero_left = basis.eval(0, 0, point[current_dim]);
//
//          if(current_dim == storage->dim()-1)
//          {
//            result.push_back(std::make_pair(seq_lz_left, value*new_value_l_zero_left));
//          }
//          else
//          {
//            rec(basis, point, current_dim + 1, value*new_value_l_zero_left, working, result);
//          }
//
//          // level 0, index 1
//          working.right_levelzero(current_dim);
//          size_t seq_lz_right = working.seq();
//          double new_value_l_zero_right = basis.eval(0, 1, point[current_dim]);
//
//          if(current_dim == storage->dim()-1)
//          {
//            result.push_back(std::make_pair(seq_lz_right, value*new_value_l_zero_right));
//          }
//          else
//          {
//            rec(basis, point, current_dim + 1, value*new_value_l_zero_right, working, result);
//          }
//        }
//      }
//
//      // there are no levels left
//      if(working.hint())
//      {
//        break;
//      }
//
//      // this decides in which direction we should descend by evaluating the corresponding bit
//      // the bits are coded from left to right starting with level 1 being in position src_level
//      if (work_level > 0)
//      {
//        double hat = 0.0;
//        level_type h = 0;
//
//        h = 1<<work_level;
//
//        hat = (1.0/static_cast<double>(h))*static_cast<double>(global_work_index);
//
//        if (point[current_dim] == hat)
//          break;
//
//        if(point[current_dim] < hat)
//        {
//          working.left_child(current_dim);
//        }
//        else
//        {
//          working.right_child(current_dim);
//        }
//      }
//      else
//      {
//        if (point[current_dim] == 0.0 || point[current_dim] == 1.0)
//          break;
//
//        working.top(current_dim);
//      }
//      ++work_level;
//    }
//
//    working.left_levelzero(current_dim);
//  }
//
//
//  void recBB(SLinearBoundaryBase& basis, std::vector<double>& point, size_t current_dim, double value, GridStorage::grid_iterator& working, std::vector<std::pair<size_t, double> >& result)
//  {
//    typedef GridStorage::index_type::level_type level_type;
//    typedef GridStorage::index_type::index_type index_type;
//
//    level_type work_level = 0;
//
//    while(true)
//    {
//      size_t seq = working.seq();
//      index_type global_work_index = 0;
//
//
//      if (storage->end(seq))
//      {
//        break;
//      }
//      else
//      {
//        index_type work_index;
//        level_type temp;
//
//        working.get(current_dim, temp, work_index);
//        global_work_index = work_index;
//
//        if (work_level > 0)
//        {
//          double new_value = basis.eval(work_level, work_index, point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));
//
//          if(current_dim == storage->dim()-1)
//          {
//            result.push_back(std::make_pair(seq, value*new_value));
//          }
//          else
//          {
//            recBB(basis, point, current_dim + 1, value*new_value, working, result);
//          }
//        }
//        // handle boundaries if we are on level 0
//        else
//        {
//          // level 0, index 0
//          working.left_levelzero(current_dim);
//          size_t seq_lz_left = working.seq();
//          double new_value_l_zero_left = basis.eval(0, 0,  point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));
//
//          if(current_dim == storage->dim()-1)
//          {
//            result.push_back(std::make_pair(seq_lz_left, value*new_value_l_zero_left));
//          }
//          else
//          {
//            recBB(basis, point, current_dim + 1, value*new_value_l_zero_left, working, result);
//          }
//
//          // level 0, index 1
//          working.right_levelzero(current_dim);
//          size_t seq_lz_right = working.seq();
//          double new_value_l_zero_right = basis.eval(0, 1, point[current_dim], BB->getIntervalWidth(current_dim), BB->getIntervalOffset(current_dim));
//
//          if(current_dim == storage->dim()-1)
//          {
//            result.push_back(std::make_pair(seq_lz_right, value*new_value_l_zero_right));
//          }
//          else
//          {
//            recBB(basis, point, current_dim + 1, value*new_value_l_zero_right, working, result);
//          }
//        }
//      }
//
//      // there are no levels left
//      if(working.hint())
//      {
//        break;
//      }
//
//      // this decides in which direction we should descend by evaluating the corresponding bit
//      // the bits are coded from left to right starting with level 1 being in position src_level
//      if (work_level > 0)
//      {
//        double hat = 0.0;
//        level_type h = 0;
//
//        h = 1<<work_level;
//
//        hat = (BB->getIntervalWidth(current_dim)*((1.0/static_cast<double>(h))*static_cast<double>(global_work_index))) + BB->getIntervalOffset(current_dim);
//
//        if (point[current_dim] == hat)
//          break;
//
//        if(point[current_dim] < hat)
//        {
//          working.left_child(current_dim);
//        }
//        else
//        {
//          working.right_child(current_dim);
//        }
//      }
//      else
//      {
//        if (point[current_dim] == (BB->getIntervalOffset(current_dim)) || point[current_dim] == ((BB->getIntervalWidth(current_dim))+BB->getIntervalOffset(current_dim)))
//          break;
//
//        working.top(current_dim);
//      }
//      ++work_level;
//    }
//
//    working.left_levelzero(current_dim);
//  }
//
//};

}
}

#endif /* ALGORITHMEVALUATIONTRANSPOSED_HPP */