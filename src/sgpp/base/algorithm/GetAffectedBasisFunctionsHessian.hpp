/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef GETAFFECTEDBASISFUNCTIONSHESSIAN_HPP
#define GETAFFECTEDBASISFUNCTIONSHESSIAN_HPP

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
class GetAffectedBasisFunctionsHessian
{
public:
    GetAffectedBasisFunctionsHessian(GridStorage *storage) : storage(storage) {}
    
    ~GetAffectedBasisFunctionsHessian() {}
    
    void operator()(BASIS& basis, const std::vector<double> &point,
                    std::vector<std::tuple<size_t, double, DataVector, DataMatrix> > &result)
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
        
        DataMatrix hessian(dim, dim);
        hessian.setAll(1.0);
        
        result.clear();
        rec(basis, point, 0, 1.0, gradient, hessian, working, source, result);
        
        delete[] source;
    }
    
protected:
    GridStorage* storage;
    
    inline void updateHessian(DataVector &gradient, DataMatrix &hessian,
                              size_t current_dim, double new_value, double new_dx, double new_dxdx,
                              DataVector &new_gradient, DataMatrix &new_hessian)
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
            
            for (size_t j = 0; j < storage->dim(); j++)
            {
                if ((i == current_dim) || (j == current_dim))
                {
                    if (i == j)
                    {
                        new_hessian.set(i, j, new_dxdx * hessian.get(i, j));
                    } else
                    {
                        new_hessian.set(i, j, new_dx * hessian.get(i, j));
                    }
                } else
                {
                    new_hessian.set(i, j, new_value * hessian.get(i, j));
                }
            }
        }
    }
    
    void rec(BASIS &basis, const std::vector<double> &point, size_t current_dim,
             double value, DataVector &gradient, DataMatrix &hessian,
             GridStorage::grid_iterator &working,
             const GridStorage::index_type::index_type *source,
             std::vector<std::tuple<size_t, double, DataVector, DataMatrix> > &result)
    {
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;
        
        level_type src_level = static_cast<level_type>(sizeof(index_type) * 8 - 1);
        index_type src_index = source[current_dim];
        
        level_type work_level = 1;
        DataVector new_gradient(storage->dim());
        DataMatrix new_hessian(storage->dim(), storage->dim());
        
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
                double new_dxdx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                updateHessian(gradient, hessian, current_dim, new_value, new_dx, new_dxdx,
                              new_gradient, new_hessian);
                
                if (current_dim == storage->dim() - 1)
                {
                    result.push_back(std::make_tuple(seq, value * new_value,
                                                     new_gradient, new_hessian));
                } else
                {
                    rec(basis, point, current_dim + 1,
                        value * new_value, new_gradient, new_hessian,
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
class GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders
{
public:
    GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders(GridStorage *storage) :
        storage(storage) {}
    
    ~GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders() {}
    
    void operator()(BASIS &base, const std::vector<double> &point,
                    std::vector<std::tuple<size_t, double, DataVector, DataMatrix> > &result)
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
        
        DataMatrix hessian(dim, dim);
        hessian.setAll(1.0);
        
        result.clear();
        rec(base, point, 0, 1.0, gradient, hessian, working, source, result);
        
        delete[] source;
    }
    
protected:
    GridStorage* storage;
    
    inline void updateHessian(DataVector &gradient, DataMatrix &hessian,
                              size_t current_dim, double new_value, double new_dx, double new_dxdx,
                              DataVector &new_gradient, DataMatrix &new_hessian)
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
            
            for (size_t j = 0; j < storage->dim(); j++)
            {
                if ((i == current_dim) || (j == current_dim))
                {
                    if (i == j)
                    {
                        new_hessian.set(i, j, new_dxdx * hessian.get(i, j));
                    } else
                    {
                        new_hessian.set(i, j, new_dx * hessian.get(i, j));
                    }
                } else
                {
                    new_hessian.set(i, j, new_value * hessian.get(i, j));
                }
            }
        }
    }
    
    void rec(BASIS &basis, const std::vector<double> &point, size_t current_dim,
             double value, DataVector &gradient, DataMatrix &hessian,
             GridStorage::grid_iterator &working,
             const GridStorage::index_type::index_type *source,
             std::vector<std::tuple<size_t, double, DataVector, DataMatrix> > &result)
    {
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;
        double tmpValue, tmpDx, tmpDxDx;
        size_t tmpSeq;
        
        level_type src_level = static_cast<level_type>(sizeof(index_type) * 8 - 1);
        index_type src_index = source[current_dim];
        
        level_type work_level = 1;
        DataVector new_gradient(storage->dim());
        DataMatrix new_hessian(storage->dim(), storage->dim());
        
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
                double new_dxdx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                updateHessian(gradient, hessian, current_dim, new_value, new_dx, new_dxdx,
                              new_gradient, new_hessian);
                
                if (current_dim == storage->dim() - 1)
                {
                    result.push_back(std::make_tuple(seq, value * new_value,
                                                     new_gradient, new_hessian));
                    
                    working.step_right(current_dim);
                    tmpSeq = working.seq();
                    
                    if (!storage->end(tmpSeq))
                    {
                        working.get(current_dim, temp, work_index);
                        tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                        tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                        tmpDxDx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                        updateHessian(gradient, hessian, current_dim, tmpValue, tmpDx, tmpDxDx,
                                      new_gradient, new_hessian);
                        result.push_back(std::make_tuple(tmpSeq, value * tmpValue,
                                                         new_gradient, new_hessian));
                    }
                    
                    working.step_left(current_dim);
                    
                    working.step_left(current_dim);
                    tmpSeq = working.seq();
                    
                    if (!storage->end(tmpSeq))
                    {
                        working.get(current_dim, temp, work_index);
                        tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                        tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                        tmpDxDx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                        updateHessian(gradient, hessian, current_dim, tmpValue, tmpDx, tmpDxDx,
                                      new_gradient, new_hessian);
                        result.push_back(std::make_tuple(tmpSeq, value * tmpValue,
                                                         new_gradient, new_hessian));
                    }
                    
                    working.step_right(current_dim);
                } else
                {
                    rec(basis, point, current_dim + 1, value * new_value,
                        new_gradient, new_hessian, working, source, result);
                    
                    working.step_right(current_dim);
                    working.get(current_dim, temp, work_index);
                    new_value = basis.eval(work_level, work_index, point[current_dim]);
                    new_dx = basis.evalDx(work_level, work_index, point[current_dim]);
                    new_dxdx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                    updateHessian(gradient, hessian, current_dim, new_value, new_dx, new_dxdx,
                                  new_gradient, new_hessian);
                    rec(basis, point, current_dim + 1, value * new_value,
                        new_gradient, new_hessian, working, source, result);
                    working.step_left(current_dim);
                    
                    working.step_left(current_dim);
                    working.get(current_dim, temp, work_index);
                    new_value = basis.eval(work_level, work_index, point[current_dim]);
                    new_dx = basis.evalDx(work_level, work_index, point[current_dim]);
                    new_dxdx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                    updateHessian(gradient, hessian, current_dim, new_value, new_dx, new_dxdx,
                                  new_gradient, new_hessian);
                    rec(basis, point, current_dim + 1, value * new_value,
                        new_gradient, new_hessian, working, source, result);
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
class GetAffectedBasisFunctionsHessian<SBsplineBase> :
        public GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders<SBsplineBase>
{
public:
    GetAffectedBasisFunctionsHessian(GridStorage* storage) :
        GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders(storage) {}
};

template<>
class GetAffectedBasisFunctionsHessian<SModBsplineBase> :
        public GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders<SModBsplineBase>
{
public:
    GetAffectedBasisFunctionsHessian(GridStorage* storage) :
        GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders(storage) {}
};

template<>
class GetAffectedBasisFunctionsHessian<SWaveletBase> :
        public GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders<SWaveletBase>
{
public:
    GetAffectedBasisFunctionsHessian(GridStorage* storage) :
        GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders(storage) {}
};

template<>
class GetAffectedBasisFunctionsHessian<SModWaveletBase> :
        public GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders<SModWaveletBase>
{
public:
    GetAffectedBasisFunctionsHessian(GridStorage* storage) :
        GetAffectedBasisFunctionsHessianWaveletBsplineNoBorders(storage) {}
};

template<class BASIS>
class GetAffectedBasisFunctionsHessianWaveletBsplineWithBorders
{
public:
    GetAffectedBasisFunctionsHessianWaveletBsplineWithBorders(GridStorage *storage) :
        storage(storage) {}
    
    ~GetAffectedBasisFunctionsHessianWaveletBsplineWithBorders() {}
    
    void operator()(BASIS &basis, const std::vector<double> &point,
                    std::vector<std::tuple<size_t, double, DataVector, DataMatrix> > &result)
    {
        GridStorage::grid_iterator working(storage);
        
        working.resetToLevelZero();
        result.clear();
        
        size_t dim = storage->dim();
        DataVector gradient(dim);
        gradient.setAll(1.0);
        
        DataMatrix hessian(dim, dim);
        hessian.setAll(1.0);
        
        rec(basis, point, 0, 1.0, gradient, hessian, working, result);
    }
    
protected:
    GridStorage* storage;
    
    virtual double getGridPoint(GridStorage::index_type::level_type level,
                                GridStorage::index_type::index_type index)
    {
        return (1.0 / static_cast<double>(1 << level)) * static_cast<double>(index);
    }
    
    inline void updateHessian(DataVector &gradient, DataMatrix &hessian,
                              size_t current_dim, double new_value, double new_dx, double new_dxdx,
                              DataVector &new_gradient, DataMatrix &new_hessian)
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
            
            for (size_t j = 0; j < storage->dim(); j++)
            {
                if ((i == current_dim) || (j == current_dim))
                {
                    if (i == j)
                    {
                        new_hessian.set(i, j, new_dxdx * hessian.get(i, j));
                    } else
                    {
                        new_hessian.set(i, j, new_dx * hessian.get(i, j));
                    }
                } else
                {
                    new_hessian.set(i, j, new_value * hessian.get(i, j));
                }
            }
        }
    }
    
    void rec(BASIS &basis, const std::vector<double> &point, size_t current_dim,
             double value, DataVector &gradient, DataMatrix &hessian,
             GridStorage::grid_iterator &working,
             std::vector<std::tuple<size_t, double, DataVector, DataMatrix> > &result)
    {
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;
        double tmpValue, tmpDx, tmpDxDx;
        size_t tmpSeq;
        
        level_type work_level = 0;
        DataVector new_gradient(storage->dim());
        DataMatrix new_hessian(storage->dim(), storage->dim());
        
        while (true)
        {
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
                    double new_dxdx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                    updateHessian(gradient, hessian, current_dim, new_value, new_dx, new_dxdx,
                                  new_gradient, new_hessian);
                    
                    if (current_dim == storage->dim() - 1)
                    {
                        result.push_back(std::make_tuple(seq, value * new_value,
                                                         new_gradient, new_hessian));
                        
                        working.step_right(current_dim);
                        tmpSeq = working.seq();
                        
                        if (!storage->end(tmpSeq))
                        {
                            working.get(current_dim, temp, work_index);
                            tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                            tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                            tmpDxDx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                            updateHessian(gradient, hessian, current_dim, tmpValue, tmpDx, tmpDxDx,
                                          new_gradient, new_hessian);
                            result.push_back(std::make_tuple(tmpSeq, value * tmpValue,
                                                             new_gradient, new_hessian));
                        }
                        
                        working.step_left(current_dim);
                        
                        working.step_left(current_dim);
                        tmpSeq = working.seq();
                        
                        if (!storage->end(tmpSeq))
                        {
                            working.get(current_dim, temp, work_index);
                            tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                            tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                            tmpDxDx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                            updateHessian(gradient, hessian, current_dim, tmpValue, tmpDx, tmpDxDx,
                                          new_gradient, new_hessian);
                            result.push_back(std::make_tuple(tmpSeq, value * tmpValue,
                                                             new_gradient, new_hessian));
                        }
                        
                        working.step_right(current_dim);
                    } else
                    {
                        rec(basis, point, current_dim + 1, value * new_value,
                            new_gradient, new_hessian, working, result);
                        
                        working.step_right(current_dim);
                        tmpSeq = working.seq();
                        
                        if (!storage->end(tmpSeq))
                        {
                            working.get(current_dim, temp, work_index);
                            tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                            tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                            tmpDxDx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                            updateHessian(gradient, hessian, current_dim, tmpValue, tmpDx, tmpDxDx,
                                          new_gradient, new_hessian);
                            rec(basis, point, current_dim + 1, value * tmpValue,
                                new_gradient, new_hessian, working, result);
                        }
                        
                        working.step_left(current_dim);
                        
                        working.step_left(current_dim);
                        tmpSeq = working.seq();
                        
                        if (!storage->end(tmpSeq))
                        {
                            working.get(current_dim, temp, work_index);
                            tmpValue = basis.eval(work_level, work_index, point[current_dim]);
                            tmpDx = basis.evalDx(work_level, work_index, point[current_dim]);
                            tmpDxDx = basis.evalDxDx(work_level, work_index, point[current_dim]);
                            updateHessian(gradient, hessian, current_dim, tmpValue, tmpDx, tmpDxDx,
                                          new_gradient, new_hessian);
                            rec(basis, point, current_dim + 1, value * tmpValue,
                                new_gradient, new_hessian, working, result);
                        }
                        
                        working.step_right(current_dim);
                    }
                } else
                {
                    working.left_levelzero(current_dim);
                    size_t seq_lz_left = working.seq();
                    double new_value_l_zero_left = basis.eval(0, 0, point[current_dim]);
                    double new_dx_l_zero_left = basis.evalDx(0, 0, point[current_dim]);
                    double new_dxdx_l_zero_left = basis.evalDxDx(0, 0, point[current_dim]);
                    updateHessian(gradient, hessian, current_dim, new_value_l_zero_left,
                                  new_dx_l_zero_left, new_dxdx_l_zero_left,
                                  new_gradient, new_hessian);
                    
                    if (current_dim == storage->dim() - 1)
                    {
                        result.push_back(std::make_tuple(
                                seq_lz_left, value * new_value_l_zero_left,
                                new_gradient, new_hessian));
                    } else
                    {
                        rec(basis, point, current_dim + 1, value * new_value_l_zero_left,
                            new_gradient, new_hessian, working, result);
                    }
                    
                    working.right_levelzero(current_dim);
                    size_t seq_lz_right = working.seq();
                    double new_value_l_zero_right = basis.eval(0, 1, point[current_dim]);
                    double new_dx_l_zero_right = basis.evalDx(0, 1, point[current_dim]);
                    double new_dxdx_l_zero_right = basis.evalDxDx(0, 1, point[current_dim]);
                    updateHessian(gradient, hessian, current_dim, new_value_l_zero_right,
                                  new_dx_l_zero_right, new_dxdx_l_zero_right,
                                  new_gradient, new_hessian);
                    
                    if (current_dim == storage->dim() - 1)
                    {
                        result.push_back(std::make_tuple(
                                seq_lz_right, value * new_value_l_zero_right,
                                new_gradient, new_hessian));
                    } else
                    {
                        rec(basis, point, current_dim + 1, value * new_value_l_zero_right,
                            new_gradient, new_hessian, working, result);
                    }
                }
            }
            
            if (working.hint())
            {
                break;
            }
            
            if (work_level > 0)
            {
                double hat = 0.0;
                level_type h = 0;
                
                h = 1 << work_level;
                hat = (1.0 / static_cast<double>(h)) * static_cast<double>(global_work_index);
                
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
    }
};

template<>
class GetAffectedBasisFunctionsHessian<SBsplineBoundaryBase> :
        public GetAffectedBasisFunctionsHessianWaveletBsplineWithBorders<SBsplineBoundaryBase>
{
public:
    GetAffectedBasisFunctionsHessian(GridStorage* storage) :
        GetAffectedBasisFunctionsHessianWaveletBsplineWithBorders(storage) {}
};

template<>
class GetAffectedBasisFunctionsHessian<SBsplineClenshawCurtisBase> :
        public GetAffectedBasisFunctionsHessianWaveletBsplineWithBorders<
                SBsplineClenshawCurtisBase>
{
public:
    GetAffectedBasisFunctionsHessian(GridStorage* storage) :
        GetAffectedBasisFunctionsHessianWaveletBsplineWithBorders(storage) {}
protected:
    virtual double getGridPoint(GridStorage::index_type::level_type level,
                                GridStorage::index_type::index_type index)
    {
        return (cos(M_PI * (1.0 - static_cast<double>(index) /
                            static_cast<double>(1 << level))) + 1.0) / 2.0;
    }
};

}
}

#endif /* GETAFFECTEDBASISFUNCTIONSHESSIAN_HPP */
