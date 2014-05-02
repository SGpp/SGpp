/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef GETAFFECTEDBASISFUNCTIONSGRADIENT_HPP
#define GETAFFECTEDBASISFUNCTIONSGRADIENT_HPP

#include "base/grid/GridStorage.hpp"
#include "base/basis/wavelet/noboundary/WaveletBasis.hpp"
#include "base/basis/wavelet/boundary/WaveletBoundaryBasis.hpp"
#include "base/basis/modwavelet/ModifiedWaveletBasis.hpp"
#include "base/basis/bspline/noboundary/BsplineBasis.hpp"
#include "base/basis/bspline/boundary/BsplineBoundaryBasis.hpp"
#include "base/basis/bspline/clenshawcurtis/BsplineClenshawCurtisBasis.hpp"
#include "base/basis/modbspline/ModifiedBsplineBasis.hpp"

#include <vector>
#include <tuple>

namespace sg
{
namespace base
{

template<class BASIS>
class GetAffectedBasisFunctionsGradient
{
public:
    GetAffectedBasisFunctionsGradient(GridStorage *storage) : storage(storage) {}
    
    ~GetAffectedBasisFunctionsGradient() {}
    
    void operator()(BASIS &basis, const std::vector<double> &point,
                    std::vector<std::tuple<size_t, double, DataVector> > &result)
    {
        GridStorage::grid_iterator working(storage);
        
        typedef GridStorage::index_type::index_type index_type;
        
        size_t bits = sizeof(index_type) * 8;
        size_t dim = storage->dim();
        index_type *source = new index_type[dim];
        
        for (size_t d = 0; d < dim; d++)
        {
            double temp = floor(point[d] * (1 << (bits - 2))) * 2;
            
            if (point[d] == 1.0)
            {
                source[d] = static_cast<index_type>(temp - 1);
            } else
            {
                source[d] = static_cast<index_type>(temp + 1);
            }
        }
        
        DataVector gradient(dim);
        gradient.setAll(1.0);
        
        result.clear();
        rec(basis, point, 0, 1.0, gradient, working, source, result);
        
        delete[] source;
    }
    
protected:
    GridStorage* storage;
    
    inline void updateGradient(DataVector &gradient, size_t current_dim, double new_value,
                               double new_dx, DataVector &new_gradient)
    {
        for (size_t i = 0; i < storage->dim(); i++)
        {
            if (i == current_dim)
            {
                new_gradient.set(i, new_dx * gradient.get(i));
            } else
            {
                new_gradient.set(i, new_value * gradient.get(i));
            }
        }
    }
    
    void rec(BASIS &basis, const std::vector<double> &point, size_t current_dim,
             double value, DataVector &gradient,
             GridStorage::grid_iterator &working,
             const GridStorage::index_type::index_type *source,
             std::vector<std::tuple<size_t, double, DataVector> > &result)
    {
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;
        
        level_type src_level = static_cast<level_type>(sizeof(index_type) * 8 - 1);
        index_type src_index = source[current_dim];
        
        level_type work_level = 1;
        DataVector new_gradient(storage->dim());
        
        while (true)
        {
            size_t seq = working.seq();
            
            if (storage->end(seq))
            {
                break;
            } else
            {
                index_type work_index;
                level_type temp;
                
                working.get(current_dim, temp, work_index);
                
                double new_value = basis.eval(work_level, work_index, point[current_dim]);
                double new_dx = basis.evalDx(work_level, work_index, point[current_dim]);
                updateGradient(gradient, current_dim, new_value, new_dx, new_gradient);
                
                if (current_dim == storage->dim() - 1)
                {
                    result.push_back(std::make_tuple(seq, value * new_value, new_gradient));
                } else
                {
                    rec(basis, point, current_dim + 1, value * new_value, new_gradient,
                        working, source, result);
                }
            }
            
            if (working.hint())
            {
                break;
            }
            
            bool right = ((src_index & (1 << (src_level - work_level))) > 0);
            work_level++;
            
            if (right)
            {
                working.right_child(current_dim);
            } else
            {
                working.left_child(current_dim);
            }
        }
        
        working.top(current_dim);
    }
};

template<class BASIS>
class GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders
{
public:
    GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders(GridStorage *storage) :
        storage(storage) {}
    
    ~GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders() {}
    
    void operator()(BASIS &base, const std::vector<double> &point,
                    std::vector<std::tuple<size_t, double, DataVector> > &result)
    {
        GridStorage::grid_iterator working(storage);
        
        typedef GridStorage::index_type::index_type index_type;
        
        size_t bits = sizeof(index_type) * 8;
        size_t dim = storage->dim();
        index_type *source = new index_type[dim];
        
        for (size_t d = 0; d < dim; d++)
        {
            double temp = floor(point[d] * (1 << (bits - 2))) * 2;
            
            if (point[d] == 1.0)
            {
                source[d] = static_cast<index_type>(temp - 1);
            } else
            {
                source[d] = static_cast<index_type>(temp + 1);
            }
        }
        
        DataVector gradient(dim);
        gradient.setAll(1.0);
        
        result.clear();
        rec(base, point, 0, 1.0, gradient, working, source, result);
        
        delete[] source;
    }
    
protected:
    GridStorage* storage;
    
    inline void updateGradient(DataVector &gradient, size_t current_dim, double new_value,
                               double new_dx, DataVector &new_gradient)
    {
        for (size_t i = 0; i < storage->dim(); i++)
        {
            if (i == current_dim)
            {
                new_gradient.set(i, new_dx * gradient.get(i));
            } else
            {
                new_gradient.set(i, new_value * gradient.get(i));
            }
        }
    }
    
    void rec(BASIS &basis, const std::vector<double> &point, size_t current_dim,
             double value, DataVector &gradient,
             GridStorage::grid_iterator &working,
             const GridStorage::index_type::index_type *source,
             std::vector<std::tuple<size_t, double, DataVector> > &result)
    {
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;
        double tmpValue, tmpDx;
        size_t tmpSeq;
        
        level_type src_level = static_cast<level_type>(sizeof(index_type) * 8 - 1);
        index_type src_index = source[current_dim];
        
        level_type work_level = 1;
        DataVector new_gradient(storage->dim());
        
        while (true)
        {
            size_t seq = working.seq();
            
            if (storage->end(seq))
            {
                break;
            } else
            {
                index_type work_index;
                level_type temp;
                
                working.get(current_dim, temp, work_index);
                
                double new_value = basis.eval(work_level, work_index, point[current_dim]);
                double new_dx = basis.evalDx(work_level, work_index, point[current_dim]);
                updateGradient(gradient, current_dim, new_value, new_dx, new_gradient);
                
                if (current_dim == storage->dim() - 1)
                {
                    result.push_back(std::make_tuple(seq, value * new_value, new_gradient));
                    
                    working.step_right(current_dim);
                    tmpSeq = working.seq();
                    
                    if (!storage->end(tmpSeq))
                    {
                        working.get(current_dim, temp, work_index);
                        tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                        tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                        updateGradient(gradient, current_dim, tmpValue, tmpDx, new_gradient);
                        result.push_back(std::make_tuple(tmpSeq, value * tmpValue, new_gradient));
                    }
                    
                    working.step_left(current_dim);
                    
                    working.step_left(current_dim);
                    tmpSeq = working.seq();
                    
                    if (!storage->end(tmpSeq))
                    {
                        working.get(current_dim, temp, work_index);
                        tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                        tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                        updateGradient(gradient, current_dim, tmpValue, tmpDx, new_gradient);
                        result.push_back(std::make_tuple(tmpSeq, value * tmpValue, new_gradient));
                    }
                    
                    working.step_right(current_dim);
                } else
                {
                    rec(basis, point, current_dim + 1, value * new_value, new_gradient,
                        working, source, result);
                    
                    working.step_right(current_dim);
                    working.get(current_dim, temp, work_index);
                    new_value = basis.eval(work_level, work_index, point[current_dim]);
                    new_dx = basis.evalDx(work_level, work_index, point[current_dim]);
                    updateGradient(gradient, current_dim, new_value, new_dx, new_gradient);
                    rec(basis, point, current_dim + 1, value * new_value, new_gradient,
                        working, source, result);
                    working.step_left(current_dim);
                    
                    working.step_left(current_dim);
                    working.get(current_dim, temp, work_index);
                    new_value = basis.eval(work_level, work_index, point[current_dim]);
                    new_dx = basis.evalDx(work_level, work_index, point[current_dim]);
                    updateGradient(gradient, current_dim, new_value, new_dx, new_gradient);
                    rec(basis, point, current_dim + 1, value * new_value, new_gradient,
                        working, source, result);
                    working.step_right(current_dim);
                }
            }
            
            if (working.hint())
            {
                break;
            }
            
            bool right = (src_index & (1 << (src_level - work_level))) > 0;
            ++work_level;
            
            if (right)
            {
                working.right_child(current_dim);
            } else
            {
                working.left_child(current_dim);
            }
        }
        
        working.top(current_dim);
    }
};

template<>
class GetAffectedBasisFunctionsGradient<SBsplineBase> :
        public GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders<SBsplineBase>
{
public:
    GetAffectedBasisFunctionsGradient(GridStorage* storage) :
        GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders(storage) {}
};

template<>
class GetAffectedBasisFunctionsGradient<SModBsplineBase> :
        public GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders<SModBsplineBase>
{
public:
    GetAffectedBasisFunctionsGradient(GridStorage* storage) :
        GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders(storage) {}
};

template<>
class GetAffectedBasisFunctionsGradient<SWaveletBase> :
        public GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders<SWaveletBase>
{
public:
    GetAffectedBasisFunctionsGradient(GridStorage* storage) :
        GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders(storage) {}
};

template<>
class GetAffectedBasisFunctionsGradient<SModWaveletBase> :
        public GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders<SModWaveletBase>
{
public:
    GetAffectedBasisFunctionsGradient(GridStorage* storage) :
        GetAffectedBasisFunctionsGradientWaveletBsplineNoBorders(storage) {}
};

template<class BASIS>
class GetAffectedBasisFunctionsGradientWaveletBsplineWithBorders
{
public:
    GetAffectedBasisFunctionsGradientWaveletBsplineWithBorders(GridStorage *storage) :
        storage(storage) {}
    
    ~GetAffectedBasisFunctionsGradientWaveletBsplineWithBorders() {}
    
    void operator()( BASIS &basis, const std::vector<double> &point,
                    std::vector<std::tuple<size_t, double, DataVector> > &result)
    {
        GridStorage::grid_iterator working(storage);
        
        working.resetToLevelZero();
        result.clear();
        
        size_t dim = storage->dim();
        DataVector gradient(dim);
        gradient.setAll(1.0);
        
        rec(basis, point, 0, 1.0, gradient, working, result);
    }
    
protected:
    GridStorage* storage;
    
    virtual double getGridPoint(GridStorage::index_type::level_type level,
                                GridStorage::index_type::index_type index)
    {
        return (1.0 / static_cast<double>(1 << level)) * static_cast<double>(index);
    }
    
    inline void updateGradient(DataVector &gradient, size_t current_dim, double new_value,
                               double new_dx, DataVector &new_gradient)
    {
        for (size_t i = 0; i < storage->dim(); i++)
        {
            if (i == current_dim)
            {
                new_gradient.set(i, new_dx * gradient.get(i));
            } else
            {
                new_gradient.set(i, new_value * gradient.get(i));
            }
        }
    }
    
    void rec(BASIS &basis, const std::vector<double> &point, size_t current_dim,
             double value, DataVector &gradient,
             GridStorage::grid_iterator &working,
             std::vector<std::tuple<size_t, double, DataVector> > &result)
    {
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;
        double tmpValue, tmpDx;
        size_t tmpSeq;
        
        level_type work_level = 0;
        DataVector new_gradient(storage->dim());
        
        //std::cout << "\nrec called with current_dim = " << current_dim << "\n";
        
        while (true)
        {
            //std::cout << "rec loop, current_dim = " << current_dim << "\n";
            size_t seq = working.seq();
            index_type global_work_index = 0;
            
            if (storage->end(seq))
            {
                break;
            } else
            {
                index_type work_index;
                level_type temp;

                working.get(current_dim, temp, work_index);
                global_work_index = work_index;
                
                if (work_level > 0)
                {
                    double new_value = basis.eval(work_level, work_index, point[current_dim]);
                    double new_dx = basis.evalDx(work_level, work_index, point[current_dim]);
                    updateGradient(gradient, current_dim, new_value, new_dx, new_gradient);
                    
                    if (current_dim == storage->dim() - 1)
                    {
                        result.push_back(std::make_tuple(seq, value * new_value, new_gradient));
                        
                        working.step_right(current_dim);
                        tmpSeq = working.seq();
                        
                        if (!storage->end(tmpSeq))
                        {
                            working.get(current_dim, temp, work_index);
                            tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                            tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                            updateGradient(gradient, current_dim, tmpValue, tmpDx, new_gradient);
                            result.push_back(std::make_tuple(tmpSeq, value * tmpValue,
                                                             new_gradient));
                        }
                        
                        working.step_left(current_dim);
                        
                        working.step_left(current_dim);
                        tmpSeq = working.seq();
                        
                        if (!storage->end(tmpSeq))
                        {
                            working.get(current_dim, temp, work_index);
                            tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                            tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                            updateGradient(gradient, current_dim, tmpValue, tmpDx, new_gradient);
                            result.push_back(std::make_tuple(tmpSeq, value * tmpValue,
                                                             new_gradient));
                        }
                        
                        working.step_right(current_dim);
                    } else
                    {
                        rec(basis, point, current_dim + 1, value * new_value, new_gradient,
                            working, result);
                        
                        working.step_right(current_dim);
                        tmpSeq = working.seq();
                        
                        if (!storage->end(tmpSeq))
                        {
                            working.get(current_dim, temp, work_index);
                            tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                            tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                            updateGradient(gradient, current_dim, tmpValue, tmpDx, new_gradient);
                            rec(basis, point, current_dim + 1, value * tmpValue, new_gradient,
                                working, result);
                        }
                        
                        working.step_left(current_dim);
                        
                        working.step_left(current_dim);
                        tmpSeq = working.seq();
                        
                        if (!storage->end(tmpSeq))
                        {
                            working.get(current_dim, temp, work_index);
                            tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                            tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                            updateGradient(gradient, current_dim, tmpValue, tmpDx, new_gradient);
                            rec(basis, point, current_dim + 1, value * tmpValue, new_gradient,
                                working, result);
                        }
                        
                        working.step_right(current_dim);
                    }
                } else
                {
                    working.left_levelzero(current_dim);
                    size_t seq_lz_left = working.seq();
                    double new_value_l_zero_left = basis.eval(0, 0, point[current_dim]);
                    double new_dx_l_zero_left = basis.evalDx(0, 0, point[current_dim]);
                    updateGradient(gradient, current_dim, new_value_l_zero_left,
                                   new_dx_l_zero_left, new_gradient);
                    
                    if (current_dim == storage->dim() - 1)
                    {
                        result.push_back(std::make_tuple(
                                seq_lz_left, value * new_value_l_zero_left, new_gradient));
                    } else
                    {
                        rec(basis, point, current_dim + 1, value * new_value_l_zero_left,
                            new_gradient, working, result);
                    }
                    
                    working.right_levelzero(current_dim);
                    size_t seq_lz_right = working.seq();
                    double new_value_l_zero_right = basis.eval(0, 1, point[current_dim]);
                    double new_dx_l_zero_right = basis.evalDx(0, 1, point[current_dim]);
                    updateGradient(gradient, current_dim, new_value_l_zero_right,
                                   new_dx_l_zero_right, new_gradient);
                    
                    if (current_dim == storage->dim() - 1)
                    {
                        result.push_back(std::make_tuple(
                                seq_lz_right, value * new_value_l_zero_right, new_gradient));
                    } else
                    {
                        rec(basis, point, current_dim + 1, value * new_value_l_zero_right,
                            new_gradient, working, result);
                    }
                }
            }
            
            if (working.hint())
            {
                break;
            }
            
            if (work_level > 0)
            {
                double hat = getGridPoint(work_level, global_work_index);
                //std::cout << "hat = " << hat << "\n";
                //std::cout << "point[current_dim] = " << point[current_dim] << "\n";
                
                /*if (point[current_dim] == hat)
                {
                    break;
                }*/
                
                if (point[current_dim] < hat)
                {
                    working.left_child(current_dim);
                } else
                {
                    working.right_child(current_dim);
                }
            } else
            {
                /*if (point[current_dim] == 0.0 || point[current_dim] == 1.0)
                {
                    break;
                }*/
                
                working.top(current_dim);
            }
            
            work_level++;
        }
        
        working.left_levelzero(current_dim);
        //std::cout << "rec exiting current_dim = " << current_dim << "\n\n";
    }
};

template<>
class GetAffectedBasisFunctionsGradient<SBsplineBoundaryBase> :
        public GetAffectedBasisFunctionsGradientWaveletBsplineWithBorders<SBsplineBoundaryBase>
{
public:
    GetAffectedBasisFunctionsGradient(GridStorage* storage) :
        GetAffectedBasisFunctionsGradientWaveletBsplineWithBorders(storage) {}
};

template<>
class GetAffectedBasisFunctionsGradient<SBsplineClenshawCurtisBase> :
        public GetAffectedBasisFunctionsGradientWaveletBsplineWithBorders<
                SBsplineClenshawCurtisBase>
{
public:
    GetAffectedBasisFunctionsGradient(GridStorage* storage) :
        GetAffectedBasisFunctionsGradientWaveletBsplineWithBorders(storage) {}
    
    void operator()(SBsplineClenshawCurtisBase &basis, const std::vector<double> &point,
                    std::vector<std::tuple<size_t, double, DataVector> > &result)
    {
        //std::cout << "\nEntering...\n";
        GetAffectedBasisFunctionsGradientWaveletBsplineWithBorders<
                SBsplineClenshawCurtisBase>::operator()(basis, point, result);
        //std::cout << "\nExiting...\n";
    }
protected:
    virtual double getGridPoint(GridStorage::index_type::level_type level,
                                GridStorage::index_type::index_type index)
    {
        //std::cout << "getGridPoint(level, index): level = " << level << ", index = "
        //          << index << "\n";
        return (cos(M_PI * (1.0 - static_cast<double>(index) /
                            static_cast<double>(1 << level))) + 1.0) / 2.0;
    }
};

template<>
class GetAffectedBasisFunctionsGradient<SWaveletBoundaryBase> :
        public GetAffectedBasisFunctionsGradientWaveletBsplineWithBorders<SWaveletBoundaryBase>
{
public:
    GetAffectedBasisFunctionsGradient(GridStorage* storage) :
        GetAffectedBasisFunctionsGradientWaveletBsplineWithBorders(storage) {}
};

}
}

#endif /* GETAFFECTEDBASISFUNCTIONSGRADIENT_HPP */
