#include "opt/gridgen/IterativeGridGeneratorFerenczi.hpp"
#include "opt/tools/Printer.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "base/grid/generation/hashmap/HashRefinementMultiple.hpp"
#include "base/grid/generation/hashmap/HashRefinementMultipleBoundaries.hpp"

#include <iostream>
#include <cmath>

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
        function::ObjectiveFunction &f, GridType grid_type, size_t N, double alpha) :
    IterativeGridGenerator(f, grid_type, N),
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

void IterativeGridGeneratorFerenczi::rankVector(const std::vector<double> &vec,
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
}

void IterativeGridGeneratorFerenczi::generate()
{
    tools::printer.printStatusBegin("Adaptive grid generation (Ferenczi)...");
    
    base::GridIndex::PointDistribution distr = ((grid_type == GridType::ClenshawCurtis) ?
            base::GridIndex::PointDistribution::ClenshawCurtis :
            base::GridIndex::PointDistribution::Normal);
    
    base::GridStorage *grid_storage = grid->getStorage();
    base::GridGenerator *grid_gen = grid->createGridGenerator();
    grid_gen->regular(static_cast<int>(INITIAL_LEVEL));
    //std::cout << "number of grid points: " << grid_storage->size() << "\n";
    
    size_t d = grid_storage->dim();
    size_t current_N = grid_storage->size();
    base::DataVector refinement_alpha(current_N);
    
    //X = std::vector<Point>(N, Point(d));
    //fX = std::vector<double>(N, 0.0);
    
    std::vector<double> x(d, 0.0);
    
    std::vector<size_t> degree(N, 0);
    std::vector<size_t> level_sum(N, 0);
    std::vector<size_t> level_max(N, 0);
    std::vector<size_t> rank(N, 0);
    function_values.assign(N, 0.0);
    
    if (current_N > N)
    {
        function_values.resize(grid_storage->size(), 0.0);
        degree.resize(grid_storage->size(), 0.0);
        level_sum.resize(grid_storage->size(), 0.0);
        level_max.resize(grid_storage->size(), 0.0);
        rank.resize(grid_storage->size(), 0.0);
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
        
        function_values[i] = f.eval(x);
    }
    
    size_t k = 0;
    size_t xhat = 0;
    
    while (current_N < N)
    {
        if (k % 10 == 0)
        {
            char str[10];
            snprintf(str, 10, "%.1f%%",
                     static_cast<double>(current_N) / static_cast<double>(N) * 100.0);
            tools::printer.printStatusUpdate(std::string(str) +
                                             " (N = " + std::to_string(N) + ")");
        }
        
        degree[xhat]++;
        refinement_alpha[xhat] = 1.0;
        
        //std::cout << "k = " << k << ", current_N = " << current_N << "\n";
        //std::cout << "should refine " << X[xhat] << " with surplus " << alpha[xhat] << "\n";
        
        base::SurplusRefinementFunctor refine_func(&refinement_alpha, 1);
        //grid_gen->refine(&refine_func);
        //base::HashRefinementSGOpt hash_refinement;
        //hash_refinement.free_refine(grid_storage, &refine_func);
        
        if ((grid_type == GridType::Boundary) || (grid_type == GridType::ClenshawCurtis))
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
            // size unchanged ==> point not refined
            break;
        }
        
        if (grid_storage->size() > N)
        {
            function_values.resize(grid_storage->size(), 0.0);
            degree.resize(grid_storage->size(), 0.0);
            level_sum.resize(grid_storage->size(), 0.0);
            level_max.resize(grid_storage->size(), 0.0);
            rank.resize(grid_storage->size(), 0.0);
        }
        
        refinement_alpha.resize(grid_storage->size());
        
        for (size_t i = 0; i < grid_storage->size(); i++)
        {
            refinement_alpha[i] = 0.0;
            
            if (i >= current_N)
            {
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
                
                function_values[i] = f.eval(x);
            }
        }
        
        current_N = grid_storage->size();
        
        std::vector<double> function_values_part(function_values.begin(),
                                                 function_values.begin() + current_N);
        rankVector(function_values_part, rank);
        
        double beta_hat = INFINITY;
        xhat = 0;
        
        for (size_t i = 0; i < current_N; i++)
        {
            // refinement criterion
            /*double beta = std::pow((double)level_sum[i] + (double)degree[i], alpha) *
                          std::pow((double)rank[i] + 1.0, 1.0 - alpha);*/
            double beta = fastPow((double)level_sum[i] + (double)degree[i], alpha) *
                          fastPow((double)rank[i] + 1.0, 1.0 - alpha);
            
            // determine the smallest beta, but make sure MAX_LEVEL isn't exceeded
            if ((beta < beta_hat) && (level_max[i] + degree[i] < MAX_LEVEL))
            {
                beta_hat = beta;
                xhat = i;
            }
        }
        
        // log((l1 + deg)^a * (r + 1)^(1-a))
        // = a * log(l1 + deg) + (1-a) * log(r + 1)
        
        k++;
    }
    
    tools::printer.printStatusUpdate("100.0% (N = " + std::to_string(current_N) + ")");
    
    function_values.erase(function_values.begin() + current_N, function_values.end());
    delete grid_gen;
    
    tools::printer.printStatusEnd();
}

}
}
}
