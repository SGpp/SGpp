/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_SLE_SYSTEM_HIERARCHISATION_HPP
#define SGPP_OPT_SLE_SYSTEM_HIERARCHISATION_HPP

#include "opt/sle/system/Cloneable.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"

#include "opt/basis/bspline/noboundary/BsplineBasis.hpp"
#include "opt/basis/bspline/boundary/BsplineBoundaryBasis.hpp"
#include "opt/basis/bspline/clenshawcurtis/BsplineClenshawCurtisBasis.hpp"
#include "opt/basis/bspline/modified/ModBsplineBasis.hpp"
#include "opt/basis/linear/noboundary/LinearBasis.hpp"
#include "opt/basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "opt/basis/linear/clenshawcurtis/LinearClenshawCurtisBasis.hpp"
#include "opt/basis/linear/modified/ModLinearBasis.hpp"
#include "opt/basis/wavelet/noboundary/WaveletBasis.hpp"
#include "opt/basis/wavelet/boundary/WaveletBoundaryBasis.hpp"
#include "opt/basis/wavelet/modified/ModWaveletBasis.hpp"

#include "opt/grid/LinearClenshawCurtisGrid.hpp"
#include "opt/grid/BsplineGrid.hpp"
#include "opt/grid/BsplineTrapezoidBoundaryGrid.hpp"
#include "opt/grid/BsplineClenshawCurtisGrid.hpp"
#include "opt/grid/ModBsplineGrid.hpp"

#include <vector>
#include <cstddef>
#include <cstring>
#include <stdexcept>

namespace sg
{
namespace opt
{
namespace sle
{
namespace system
{

/**
 * Linear system of the hierarchisation in a sparse grid.
 */
class Hierarchisation : public Cloneable
{
public:
    /**
     * Constructor.
     * Do not destruct the grid before this object!
     * 
     * @param grid              sparse grid
     */
    Hierarchisation(base::Grid &grid) :
        Cloneable(),
        grid(grid),
        grid_storage(grid.getStorage()),
        basis_type(INVALID)
    {
        initialize();
    }
    
    /**
     * Constructor.
     * Do not destruct the grid before this object!
     * 
     * @param grid              sparse grid
     * @param grid_storage      custom grid storage (use basis function according to grid,
     *                          but use another set of grid points according to grid_storage)
     */
    Hierarchisation(base::Grid &grid, base::GridStorage *grid_storage) :
        Cloneable(),
        grid(grid),
        grid_storage(grid_storage),
        basis_type(INVALID)
    {
        initialize();
    }
    
    /**
     * @param i     row index
     * @param j     column index
     * @return      whether the i-th grid point lies in the support of the j-th basis function
     */
    inline bool isMatrixEntryNonZero(size_t i, size_t j)
    {
        return (evalBasisFunctionAtGridPoint(j, i) != 0.0);
    }
    
    /**
     * @param i     row index
     * @param j     column index
     * @return      value of the j-th basis function at the i-th grid point
     */
    inline double getMatrixEntry(size_t i, size_t j)
    {
        return evalBasisFunctionAtGridPoint(j, i);
    }
    
    /**
     * @return          sparse grid
     */
    base::Grid &getGrid()
    {
        return grid;
    }
    
    /**
     * @return grid     sparse grid
     */
    void setGrid(base::Grid &grid)
    {
        this->grid = grid;
    }
    
    /**
     * @return              grid storage
     */
    base::GridStorage *getGridStorage()
    {
        return grid_storage;
    }
    
    /**
     * @param grid_storage  grid storage (do not destruct before this object!)
     */
    void setGridStorage(base::GridStorage *grid_storage)
    {
        this->grid_storage = grid_storage;
    }
    
    size_t getDimension() const
    {
        return grid_storage->size();
    }
    
    /**
     * @return smart pointer to cloned object
     */
    virtual tools::SmartPointer<Cloneable> clone()
    {
        return tools::SmartPointer<Cloneable>(new Hierarchisation(grid, grid_storage));
    }
    
protected:
    /// sparse grid
    base::Grid &grid;
    /// grid storage
    base::GridStorage *grid_storage;
    
    /// B-spline basis
    tools::SmartPointer<SBsplineBase> bspline_basis;
    /// B-spline boundary basis
    tools::SmartPointer<SBsplineBoundaryBase> bspline_boundary_basis;
    /// B-spline Clenshaw-Curtis basis
    tools::SmartPointer<SBsplineClenshawCurtisBase> bspline_clenshaw_curtis_basis;
    /// modified B-spline basis
    tools::SmartPointer<SModBsplineBase> mod_bspline_basis;
    /// linear basis
    tools::SmartPointer<SLinearBase> linear_basis;
    /// linear boundary basis
    tools::SmartPointer<SLinearBoundaryBase> linear_boundary_basis;
    /// linear Clenshaw-Curtis basis
    tools::SmartPointer<SLinearClenshawCurtisBase> linear_clenshaw_curtis_basis;
    /// modified linear basis
    tools::SmartPointer<SModLinearBase> mod_linear_basis;
    /// wavelet basis
    tools::SmartPointer<SWaveletBase> wavelet_basis;
    /// wavelet boundary basis
    tools::SmartPointer<SWaveletBoundaryBase> wavelet_boundary_basis;
    /// modified wavelet basis
    tools::SmartPointer<SModWaveletBase> mod_wavelet_basis;
    
    /// type of grid/basis functions
    enum
    {
        INVALID,
        BSPLINE, BSPLINE_BOUNDARY, BSPLINE_CLENSHAW_CURTIS, MOD_BSPLINE,
        LINEAR, LINEAR_BOUNDARY, LINEAR_CLENSHAW_CURTIS, MOD_LINEAR,
        WAVELET, WAVELET_BOUNDARY, MOD_WAVELET
    } basis_type;
    
    /**
     * Initialize the correct basis (according to the grid).
     */
    void initialize()
    {
        if (strcmp(grid.getType(), "Bspline") == 0)
        {
            bspline_basis = tools::SmartPointer<SBsplineBase>(
                    new SBsplineBase(
                    dynamic_cast<BsplineGrid &>(grid).getDegree()));
            basis_type = BSPLINE;
        } else if (strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0)
        {
            bspline_boundary_basis = tools::SmartPointer<SBsplineBoundaryBase>(
                    new SBsplineBoundaryBase(
                    dynamic_cast<BsplineTrapezoidBoundaryGrid &>(grid).getDegree()));
            basis_type = BSPLINE_BOUNDARY;
        } else if (strcmp(grid.getType(), "BsplineClenshawCurtis") == 0)
        {
            bspline_clenshaw_curtis_basis = tools::SmartPointer<SBsplineClenshawCurtisBase>(
                    new SBsplineClenshawCurtisBase(
                    dynamic_cast<BsplineClenshawCurtisGrid &>(grid).getDegree(),
                    dynamic_cast<BsplineClenshawCurtisGrid &>(grid).getCosineTable()));
            basis_type = BSPLINE_CLENSHAW_CURTIS;
        } else if (strcmp(grid.getType(), "modBspline") == 0)
        {
            mod_bspline_basis = tools::SmartPointer<SModBsplineBase>(
                    new SModBsplineBase(
                    dynamic_cast<ModBsplineGrid &>(grid).getDegree()));
            basis_type = MOD_BSPLINE;
        } else if (strcmp(grid.getType(), "linear") == 0)
        {
            linear_basis = tools::SmartPointer<SLinearBase>(
                    new SLinearBase());
            basis_type = LINEAR;
        } else if (strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
            linear_boundary_basis = tools::SmartPointer<SLinearBoundaryBase>(
                    new SLinearBoundaryBase());
            basis_type = LINEAR_BOUNDARY;
        } else if (strcmp(grid.getType(), "linearClenshawCurtis") == 0)
        {
            linear_clenshaw_curtis_basis = tools::SmartPointer<SLinearClenshawCurtisBase>(
                    new SLinearClenshawCurtisBase(
                    dynamic_cast<LinearClenshawCurtisGrid &>(grid).getCosineTable()));
            basis_type = LINEAR_CLENSHAW_CURTIS;
        } else if (strcmp(grid.getType(), "modLinear") == 0)
        {
            mod_linear_basis = tools::SmartPointer<SModLinearBase>(
                    new SModLinearBase());
            basis_type = MOD_LINEAR;
        } else if (strcmp(grid.getType(), "Wavelet") == 0)
        {
            wavelet_basis = tools::SmartPointer<SWaveletBase>(
                    new SWaveletBase());
            basis_type = WAVELET;
        } else if (strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0)
        {
            wavelet_boundary_basis = tools::SmartPointer<SWaveletBoundaryBase>(
                    new SWaveletBoundaryBase());
            basis_type = WAVELET_BOUNDARY;
        } else if (strcmp(grid.getType(), "modWavelet") == 0)
        {
            mod_wavelet_basis = tools::SmartPointer<SModWaveletBase>(
                    new SModWaveletBase());
            basis_type = MOD_WAVELET;
        } else
        {
            throw std::invalid_argument("Grid type not supported.");
        }
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th basis function at the point_j-th grid point
     */
    inline double evalBasisFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        if (basis_type == BSPLINE)
        {
            return evalBsplineFunctionAtGridPoint(basis_i, point_j);
        } else if (basis_type == BSPLINE_BOUNDARY)
        {
            return evalBsplineBoundaryFunctionAtGridPoint(basis_i, point_j);
        } else if (basis_type == BSPLINE_CLENSHAW_CURTIS)
        {
            return evalBsplineClenshawCurtisFunctionAtGridPoint(basis_i, point_j);
        } else if (basis_type == MOD_BSPLINE)
        {
            return evalModBsplineFunctionAtGridPoint(basis_i, point_j);
        } else if (basis_type == LINEAR)
        {
            return evalLinearFunctionAtGridPoint(basis_i, point_j);
        } else if (basis_type == LINEAR_BOUNDARY)
        {
            return evalLinearBoundaryFunctionAtGridPoint(basis_i, point_j);
        } else if (basis_type == LINEAR_CLENSHAW_CURTIS)
        {
            return evalLinearClenshawCurtisFunctionAtGridPoint(basis_i, point_j);
        } else if (basis_type == MOD_LINEAR)
        {
            return evalModLinearFunctionAtGridPoint(basis_i, point_j);
        } else if (basis_type == WAVELET)
        {
            return evalWaveletFunctionAtGridPoint(basis_i, point_j);
        } else if (basis_type == WAVELET_BOUNDARY)
        {
            return evalWaveletBoundaryFunctionAtGridPoint(basis_i, point_j);
        } else if (basis_type == MOD_WAVELET)
        {
            return evalModWaveletFunctionAtGridPoint(basis_i, point_j);
        } else
        {
            return 0.0;
        }
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th B-spline basis function
     *                  at the point_j-th grid point
     */
    inline double evalBsplineFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = bspline_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th B-spline boundary basis function
     *                  at the point_j-th grid point
     */
    inline double evalBsplineBoundaryFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = bspline_boundary_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th B-spline Clenshaw-Curtis basis function
     *                  at the point_j-th grid point
     */
    inline double evalBsplineClenshawCurtisFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = bspline_clenshaw_curtis_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th modified B-spline basis function
     *                  at the point_j-th grid point
     */
    inline double evalModBsplineFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = mod_bspline_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th linear basis function
     *                  at the point_j-th grid point
     */
    inline double evalLinearFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = linear_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th linear boundary basis function
     *                  at the point_j-th grid point
     */
    inline double evalLinearBoundaryFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = linear_boundary_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th linear Clenshaw-Curtis basis function
     *                  at the point_j-th grid point
     */
    inline double evalLinearClenshawCurtisFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = linear_clenshaw_curtis_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th modified linear basis function
     *                  at the point_j-th grid point
     */
    inline double evalModLinearFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = mod_linear_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th wavelet basis function
     *                  at the point_j-th grid point
     */
    inline double evalWaveletFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = wavelet_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th wavelet boundary basis function
     *                  at the point_j-th grid point
     */
    inline double evalWaveletBoundaryFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = wavelet_boundary_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
    
    /**
     * @param basis_i   basis function index
     * @param point_j   grid point index
     * @return          value of the basis_i-th modified wavelet basis function
     *                  at the point_j-th grid point
     */
    inline double evalModWaveletFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = mod_wavelet_basis->eval(
                    gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
};

}
}
}
}

#endif
