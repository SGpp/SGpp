#include "opt/gridgen/IterativeGridGeneratorFerenczi.hpp"
#include "opt/tools/Printer.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "base/grid/generation/hashmap/HashRefinementMultiple.hpp"
#include "base/grid/generation/hashmap/HashRefinementMultipleBoundaries.hpp"

#include <iostream>
#include <cmath>
#include <cstring>

namespace sg
{
namespace opt
{
namespace gridgen
{

// http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
inline double fastPow(double a, double b)
{
    union
    {
        double d;
        int x[2];
    } u = {a};
    
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    
    return u.d;
}

IterativeGridGeneratorFerenczi::IterativeGridGeneratorFerenczi(
        function::Objective &f, base::Grid &grid, size_t N, double alpha) :
    IterativeGridGenerator(f, grid, N),
    alpha(alpha)
{
}

double IterativeGridGeneratorFerenczi::getAlpha() const
{
    return alpha;
}

void IterativeGridGeneratorFerenczi::setAlpha(double alpha)
{
    this->alpha = alpha;
}

/*void IterativeGridGeneratorFerenczi::rankVector(const std::vector<double> &vec,
                                                std::vector<size_t> &rank) const
{
    size_t n = vec.size();
    std::vector<size_t> index(n, 0);
    rank.assign(n, 0);
    
    for (size_t i = 0; i < n; i++)
    {
        index[i] = i;
        rank[i] = i;
    }
    
    std::sort(index.begin(), index.end(),
            [&](const size_t &a, const size_t &b)
            {
                return (vec[a] < vec[b]);
            }
    );
    
    std::sort(rank.begin(), rank.end(),
            [&](const size_t &a, const size_t &b)
            {
                return (index[a] < index[b]);
            }
    );
}*/

bool IterativeGridGeneratorFerenczi::generate()
{
    tools::printer.printStatusBegin("Adaptive grid generation (Ferenczi)...");
    
    base::GridIndex::PointDistribution distr = base::GridIndex::PointDistribution::Normal;
    
    if (strcmp(grid.getType(), "BsplineClenshawCurtis") == 0)
    {
        distr = base::GridIndex::PointDistribution::ClenshawCurtis;
    }
    
    bool is_boundary_grid = ((strcmp(grid.getType(), "BsplineBoundary") == 0) ||
                             (strcmp(grid.getType(), "WaveletBoundary") == 0) ||
                             (strcmp(grid.getType(), "BsplineClenshawCurtis") == 0));
    int initial_level = (is_boundary_grid ? 1 : static_cast<int>(INITIAL_LEVEL_WO_BOUNDARIES));
    
    base::GridStorage *grid_storage = grid.getStorage();
    
    {
        std::unique_ptr<base::GridGenerator> grid_gen(grid.createGridGenerator());
        grid_gen->regular(initial_level);
    }
    //std::cout << "number of grid points: " << grid_storage->size() << "\n";
    
    size_t d = grid_storage->dim();
    size_t current_N = grid_storage->size();
    base::DataVector refinement_alpha(current_N);
    
    //X = std::vector<Point>(N, Point(d));
    //fX = std::vector<double>(N, 0.0);
    
    std::vector<double> &fX = function_values;
    std::vector<double> fX_sorted(current_N, 0);
    std::vector<size_t> fX_order(current_N, 0);
    // fX_sorted[i] = fX[fX_order[i]] (fX sorted in ascending order)
    
    std::vector<double> x(d, 0.0);
    
    std::vector<size_t> degree(N, 0);
    std::vector<size_t> level_sum(N, 0);
    std::vector<size_t> level_max(N, 0);
    std::vector<size_t> rank(N, 0);
    fX.assign(N, 0.0);
    
    if (current_N > N)
    {
        fX.resize(grid_storage->size(), 0.0);
        degree.resize(grid_storage->size(), 0);
        level_sum.resize(grid_storage->size(), 0);
        level_max.resize(grid_storage->size(), 0);
        rank.resize(grid_storage->size(), 0);
    }
    
    for (size_t i = 0; i < current_N; i++)
    {
        base::GridIndex *gp = grid_storage->get(i);
        gp->setPointDistribution(distr);
        //std::cout << "grid point " << i << ": [" << gp->abs(0) << ", " << gp->abs(1) << "] (" <<
        //             gp->toString() << ")\n";
        //grid_index_to_point(gp, X[i]);
        //fX[i] = f(X[i].coords, data);
        refinement_alpha[i] = 0.0;
        
        for (size_t t = 0; t < d; t++)
        {
            x[t] = gp->abs(t);
            
            size_t level = gp->getLevel(t);
            level_sum[i] += level;
            
            if (level_max[i] < level)
            {
                level_max[i] = level;
            }
        }
        
        fX[i] = f.eval(x);
        fX_order[i] = i;
    }
    
    std::sort(fX_order.begin(), fX_order.end(),
            [&](const size_t &a, const size_t &b)
            {
                return (fX[a] < fX[b]);
            }
    );
    
    for (size_t i = 0; i < current_N; i++)
    {
        fX_sorted[i] = fX[fX_order[i]];
    }
    
    size_t k = 0;
    size_t xhat = 0;
    bool result = true;
    
    while (current_N < N)
    {
        //if (k % 10 == 0)
        {
            char str[10];
            snprintf(str, 10, "%.1f%%",
                     static_cast<double>(current_N) / static_cast<double>(N) * 100.0);
            tools::printer.printStatusUpdate(std::string(str) +
                                             " (N = " + std::to_string(current_N) + ")");
        }
        
        degree[xhat]++;
        refinement_alpha[xhat] = 1.0;
        
        /*std::cout << "\nk = " << k << ", current_N = " << current_N << "\n";
        std::cout << "refining X[" << xhat << "] = " << grid_storage->get(xhat)
                  << ", level = [" << grid_storage->get(xhat)->getLevel(0)
                  << ", " << grid_storage->get(xhat)->getLevel(1)
                  << "], degree = " << degree[xhat] << "\n";*/
        
        base::SurplusRefinementFunctor refine_func(&refinement_alpha, 1);
        //grid_gen->refine(&refine_func);
        //base::HashRefinementSGOpt hash_refinement;
        //hash_refinement.free_refine(grid_storage, &refine_func);
        
        if (is_boundary_grid)
        {
            base::HashRefinementMultipleBoundaries hash_refinement;
            hash_refinement.free_refine(grid_storage, &refine_func);
        } else
        {
            base::HashRefinementMultiple hash_refinement;
            hash_refinement.free_refine(grid_storage, &refine_func);
        }
        
        if (grid_storage->size() == current_N)
        {
            // size unchanged ==> point not refined (should not happen)
            tools::printer.printStatusEnd(
                    "error: size unchanged in IterativeGridGeneratorFerenczi");
            result = false;
            break;
        }
        
        if (grid_storage->size() > N)
        {
            std::list<size_t> indices_to_remove;
            
            for (size_t i = current_N; i < grid_storage->size(); i++)
            {
                indices_to_remove.push_back(i);
            }
            
            grid_storage->deletePoints(indices_to_remove);
            break;
        }
        
        refinement_alpha.resize(grid_storage->size());
        refinement_alpha[xhat] = 0.0;
        
        for (size_t i = current_N; i < grid_storage->size(); i++)
        {
            refinement_alpha[i] = 0.0;
            
            base::GridIndex *gp = grid_storage->get(i);
            gp->setPointDistribution(distr);
            //grid_index_to_point(gp, X[i]);
            //fX[i] = f(X[i].coords, data);
            
            for (size_t t = 0; t < d; t++)
            {
                x[t] = gp->abs(t);
                
                size_t level = gp->getLevel(t);
                level_sum[i] += level;
                
                if (level_max[i] < level)
                {
                    level_max[i] = level;
                }
            }
            
            /*std::cout << "new point X[" << i << "] = " << gp
                      << ", level = [" << gp->getLevel(0) << ", " << gp->getLevel(1) << "]\n";*/
            
            fX[i] = f.eval(x);
            
            // update rank and fX_order by insertion sort
            // ==> go through fX from biggest entry to lowest
            for (size_t j = i - 1; j-- > 0; )
            {
                if (fX_sorted[j] < fX[i])
                {
                    // new function value is bigger as current one ==> insert here
                    fX_order.insert(fX_order.begin() + (j+1), i);
                    fX_sorted.insert(fX_sorted.begin() + (j+1), fX[i]);
                    rank[i] = j + 1;
                    break;
                } else
                {
                    // new function value value not bigger yet ==> increase rank
                    rank[fX_order[j]]++;
                }
            }
            
            if (fX_order.size() == i)
            {
                // this happens if the new function value is smaller than all of the previous ones
                // ==> insert at the beginning of fX_order, rank = 0
                fX_order.insert(fX_order.begin(), i);
                fX_sorted.insert(fX_sorted.begin(), fX[i]);
                rank[i] = 0;
            }
        }
        
        current_N = grid_storage->size();
        
        /*std::vector<double> fX_part(fX.begin(), fX.begin() + current_N);
        rankVector(fX_part, rank);*/
        
        double beta_hat = INFINITY;
        xhat = 0;
        
        for (size_t i = 0; i < current_N; i++)
        {
            // refinement criterion
            /*double beta = std::pow((double)level_sum[i] + (double)degree[i], alpha) *
                          std::pow((double)rank[i] + 1.0, 1.0 - alpha);*/
            // TODO: do static_cast instead of C-style cast
            double beta = fastPow((double)level_max[i] + (double)degree[i], alpha) *
                          fastPow((double)rank[i] + 1.0, 1.0 - alpha);
            
            // determine the smallest beta, but make sure MAX_LEVEL isn't exceeded
            //if ((beta < beta_hat) && (level_max[i] + degree[i] < MAX_LEVEL))
            if (beta < beta_hat)
            {
                base::GridIndex *gp = grid_storage->get(i);
                base::GridIndex::index_type source_index, child_index;
                base::GridIndex::level_type source_level, child_level;
                bool skip_point = false;
                
                for (size_t t = 0; t < d; t++)
                {
                    gp->get(t, source_level, source_index);
                    
                    child_index = source_index;
                    child_level = source_level;
                    
                    while (grid_storage->has_key(gp))
                    {
                        child_index *= 2;
                        child_level++;
                        gp->set(t, child_level, child_index - 1);
                    }
                    
                    gp->set(t, source_level, source_index);
                    
                    if (child_level > MAX_LEVEL)
                    {
                        skip_point = true;
                        break;
                    }
                    
                    child_index = source_index;
                    child_level = source_level;
                    
                    while (grid_storage->has_key(gp))
                    {
                        child_index *= 2;
                        child_level++;
                        gp->set(t, child_level, child_index + 1);
                    }
                    
                    gp->set(t, source_level, source_index);
                    
                    if (child_level > MAX_LEVEL)
                    {
                        skip_point = true;
                        break;
                    }
                }
                
                if (!skip_point)
                {
                    beta_hat = beta;
                    xhat = i;
                }
            }
        }
        
        k++;
    }
    
    function_values.erase(function_values.begin() + current_N, function_values.end());
    
    if (result)
    {
        tools::printer.printStatusUpdate("100.0% (N = " + std::to_string(current_N) + ")");
        tools::printer.printStatusEnd();
        return true;
    } else
    {
        return false;
    }
}

}
}
}
