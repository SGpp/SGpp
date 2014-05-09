#include "opt/gridgen/IterativeGridGeneratorLinearSurplus.hpp"
#include "opt/tools/Printer.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "base/operation/OperationHierarchisation.hpp"
#include "base/grid/generation/hashmap/HashRefinement.hpp"
#include "base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"

#include "base/basis/linear/noboundary/LinearBasis.hpp"
#include "base/basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "base/basis/linear/clenshawcurtis/LinearClenshawCurtisBasis.hpp"
#include "base/basis/modlinear/ModifiedLinearBasis.hpp"

//#include <iostream>
//#include <cmath>
#include <cstring>

namespace sg
{
namespace opt
{
namespace gridgen
{

IterativeGridGeneratorLinearSurplus::IterativeGridGeneratorLinearSurplus(
        function::Objective &f, base::Grid &grid, size_t N, double alpha,
        const base::CosineTable *cosine_table) :
    IterativeGridGenerator(f, grid, N),
    linear_base(NULL),
    linear_grid(NULL),
    alpha(alpha)
{
    if ((strcmp(grid.getType(), "Bspline") == 0) ||
        (strcmp(grid.getType(), "Wavelet") == 0))
    {
        linear_base = new base::SLinearBase();
        linear_grid = new base::LinearGrid(f.getDimension());
    } else if ((strcmp(grid.getType(), "BsplineBoundary") == 0) ||
        (strcmp(grid.getType(), "WaveletBoundary") == 0))
    {
        linear_base = new base::SLinearBoundaryBase();
        linear_grid = new base::LinearBoundaryGrid(f.getDimension());
    } else if (strcmp(grid.getType(), "BsplineClenshawCurtis") == 0)
    {
        linear_base = new base::SLinearClenshawCurtisBase(cosine_table);
        linear_grid = new base::LinearClenshawCurtisGrid(f.getDimension());
    } else if ((strcmp(grid.getType(), "modBspline") == 0) ||
        (strcmp(grid.getType(), "modWavelet") == 0))
    {
        linear_base = new base::SModLinearBase();
        linear_grid = new base::ModLinearGrid(f.getDimension());
    } else
    {
        throw std::invalid_argument("Grid type not supported.");
    }
}

IterativeGridGeneratorLinearSurplus::~IterativeGridGeneratorLinearSurplus()
{
    if (linear_base != NULL)
    {
        delete linear_base;
    }
    
    if (linear_grid != NULL)
    {
        delete linear_grid;
    }
}

double IterativeGridGeneratorLinearSurplus::getAlpha() const
{
    return alpha;
}

void IterativeGridGeneratorLinearSurplus::setAlpha(double alpha)
{
    this->alpha = alpha;
}

inline double IterativeGridGeneratorLinearSurplus::evalBasisFunctionAtGridPoint(
        base::GridStorage *grid_storage, size_t basis_i, size_t point_j)
{
    base::GridIndex *gp_basis = grid_storage->get(basis_i);
    base::GridIndex *gp_point = grid_storage->get(point_j);
    size_t d = grid_storage->dim();
    double result = 1.0;
    
    for (size_t t = 0; t < d; t++)
    {
        double result1d = linear_base->eval(
                gp_basis->getLevel(t), gp_basis->getIndex(t), gp_point->abs(t));
        
        if (result1d == 0.0)
        {
            return 0.0;
        }
        
        result *= result1d;
    }
    
    return result;
}

void IterativeGridGeneratorLinearSurplus::generate()
{
    tools::printer.printStatusBegin("Adaptive grid generation (linear surplus)...");
    
    base::GridIndex::PointDistribution distr = base::GridIndex::PointDistribution::Normal;
    
    if (strcmp(grid.getType(), "BsplineClenshawCurtis") == 0)
    {
        distr = base::GridIndex::PointDistribution::ClenshawCurtis;
    }
    
    base::HashRefinement hash_refinement;
    base::HashRefinementBoundaries hash_refinement_boundaries;
    
    base::AbstractRefinement *abstract_refinement;
    
    if ((strcmp(grid.getType(), "BsplineBoundary") == 0) ||
        (strcmp(grid.getType(), "WaveletBoundary") == 0) ||
        (strcmp(grid.getType(), "BsplineClenshawCurtis") == 0))
    {
        abstract_refinement = &hash_refinement_boundaries;
    } else
    {
        abstract_refinement = &hash_refinement;
    }
    
    base::GridStorage *grid_storage = grid.getStorage();
    
    base::GridGenerator *grid_gen = grid.createGridGenerator();
    grid_gen->regular(static_cast<int>(INITIAL_LEVEL));
    delete grid_gen;
    
    base::GridGenerator *linear_grid_gen = linear_grid->createGridGenerator();
    linear_grid_gen->regular(static_cast<int>(INITIAL_LEVEL));
    delete linear_grid_gen;
    
    size_t d = grid_storage->dim();
    size_t current_N = grid_storage->size();
    base::DataVector coeffs(current_N);
    
    std::vector<double> x(d, 0.0);
    function_values.assign(N, 0.0);
    
    if (current_N > N)
    {
        function_values.resize(grid_storage->size(), 0.0);
    }
    
    for (size_t i = 0; i < current_N; i++)
    {
        base::GridIndex *gp = grid_storage->get(i);
        gp->setPointDistribution(distr);
        
        linear_grid->getStorage()->get(i)->setPointDistribution(distr);
        
        for (size_t t = 0; t < d; t++)
        {
            x[t] = gp->abs(t);
        }
        
        function_values[i] = f.eval(x);
        coeffs[i] = function_values[i];
    }
    
    base::OperationHierarchisation *linear_hier_op =
            op_factory::createOperationHierarchisation(*linear_grid);
    linear_hier_op->doHierarchisation(coeffs);
    delete linear_hier_op;
    
    size_t k = 0;
    
    size_t refinable_pts_count;
    size_t pts_to_be_refined_count;
    double refine_factor = 1.0;
    
    while (true)
    {
        {
            char str[10];
            snprintf(str, 10, "%.1f%%",
                     static_cast<double>(current_N) / static_cast<double>(N) * 100.0);
            tools::printer.printStatusUpdate(std::string(str) +
                                             " (N = " + std::to_string(current_N) + ")");
            //tools::printer.printStatusNewLine();
            //std::cout << "pts_to_be_refined_count = " << pts_to_be_refined_count << "\n";
        }
        
        refinable_pts_count = abstract_refinement->getNumberOfRefinablePoints(grid_storage);
        pts_to_be_refined_count = static_cast<int>(
                1.0 + refine_factor * alpha * static_cast<double>(refinable_pts_count));
        
        base::SurplusRefinementFunctor refine_func(&coeffs, pts_to_be_refined_count);
        abstract_refinement->free_refine(grid_storage, &refine_func);
        
        /*base::SurplusRefinementFunctor refine_func(&refinement_alpha, pts_to_be_refined_count);
        //grid_gen->refine(&refine_func);
        base::HashRefinementSGOpt hash_refinement;
        hash_refinement.free_refine(grid_storage, &refine_func);*/
        
        if (grid_storage->size() == current_N)
        {
            // size unchanged ==> nothing refined
            std::cout << "IterativeGridGeneratorLinearSurplus: size unchanged\n";
            break;
        }
        
        if (grid_storage->size() > N)
        {
            //function_values.resize(grid_storage->size(), 0.0);
            std::list<size_t> indices_to_remove;
            
            for (size_t i = current_N; i < grid_storage->size(); i++)
            {
                indices_to_remove.push_back(i);
            }
            
            grid_storage->deletePoints(indices_to_remove);
            
            if (pts_to_be_refined_count == 1)
            {
                break;
            } else
            {
                refine_factor /= 2.0;
                k++;
                continue;
            }
        }
        
        coeffs.resize(grid_storage->size());
        
        for (size_t i = current_N; i < grid_storage->size(); i++)
        {
            base::GridIndex *gp = grid_storage->get(i);
            gp->setPointDistribution(distr);
            
            for (size_t t = 0; t < d; t++)
            {
                x[t] = gp->abs(t);
            }
            
            function_values[i] = f.eval(x);
            coeffs[i] = function_values[i];
            
            for (uint j = 0; j < i; j++)
            {
                coeffs[i] -= evalBasisFunctionAtGridPoint(grid_storage, j, i) * coeffs[j];
            }
        }
        
        current_N = grid_storage->size();
        k++;
    }
    
    tools::printer.printStatusUpdate("100.0% (N = " + std::to_string(current_N) + ")");
    
    function_values.erase(function_values.begin() + current_N, function_values.end());
    
    tools::printer.printStatusEnd();
}

/*void IterativeGridGeneratorLinearSurplus::generate()
{
    base::GridStorage *grid_storage = grid->getStorage();
    base::GridGenerator *grid_gen = grid->createGridGenerator();
    grid_gen->regular(static_cast<int>(INITIAL_LEVEL));
    
    size_t current_N = grid_storage->size();
    base::DataVector coeffs(current_N);
    base::DataVector coeffs_refine(current_N);
    size_t refinable_pts_count = 0;
    
    std::vector<double> x(d, 0.0);
    function_values.assign(N, 0.0);
    
    if (current_N > N)
    {
        function_values.resize(grid_storage->size(), 0.0);
    }
    
    for (size_t i = 0; i < current_N; i++)
    {
        base::GridIndex *gp = grid_storage->get(i);
        
        for (size_t t = 0; t < d; t++)
        {
            x[t] = gp->abs(t);
        }
        
        function_values[i] = f.eval(x);
        coeffs[i] = function_values[i];
    }
    
    base::OperationHierarchisation *hier_op = op_factory::createOperationHierarchisation(*grid);
    hier_op->doHierarchisation(coeffs);
    
    for (size_t i = 0; i < current_N; i++)
    {
        base::GridIndex *gp = grid_storage->get(i);
        
        if (gp->getLevelSum() == INITIAL_LEVEL + d - 1)
        {
            refinable_pts_count++;
            coeffs_refine[i] = coeffs[i];
        } else
        {
            coeffs_refine[i] = 0.0;
        }
    }
    
    base::HashRefinement hash_refinement;
    size_t k = 0;
    
    while (current_N < N)
    {
        size_t pts_to_be_refined_count =
                static_cast<int>(1.0 + alpha * static_cast<double>(refinable_pts_count));
        base::SurplusRefinementFunctor refine_func(&coeffs_refine, pts_to_be_refined_count);
        hash_refinement.free_refine(grid_storage, &refine_func);
        
        //base::SurplusRefinementFunctor refine_func(&refinement_alpha, pts_to_be_refined_count);
        ////grid_gen->refine(&refine_func);
        //base::HashRefinementSGOpt hash_refinement;
        //hash_refinement.free_refine(grid_storage, &refine_func);
        
        if (grid_storage->size() == current_N)
        {
            // size unchanged ==> nothing refined
            break;
        }
        
        if (grid_storage->size() > N)
        {
            function_values.resize(grid_storage->size(), 0.0);
        }
        
        coeffs.resizeZero(grid_storage->size());
        coeffs_refine.resize(grid_storage->size());
        coeffs_refine.setAll(0.0);
        
        for (size_t i = current_N; i < grid_storage->size(); i++)
        {
            base::GridIndex *gp = grid_storage->get(i);
            //grid_index_to_point(gp, X[i]);
            //fX[i] = f(X[i].coords, data);
            
            for (size_t t = 0; t < d; t++)
            {
                x[t] = gp->abs(t);
            }
            
            function_values[i] = f.eval(x);
            coeffs[i] = function_values[i];
            
            for (uint j = 0; j < i; j++)
            {
                coeffs[i] -= evalBasisFunctionAtGridPoint(grid_storage, j, i) * coeffs[j];
            }
            
            coeffs_refine[i] = coeffs[i];
        }
        
        refinable_pts_count = grid_storage->size() - current_N;
        current_N = grid_storage->size();
        
        k++;
    }
    
    function_values.erase(function_values.begin() + current_N, function_values.end());
    delete hier_op;
    delete grid_gen;
}*/

}
}
}
