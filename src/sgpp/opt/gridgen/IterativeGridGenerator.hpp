#ifndef SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP
#define SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP

#include "opt/function/ObjectiveFunction.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/type/LinearBoundaryGrid.hpp"
#include "base/grid/type/LinearClenshawCurtisGrid.hpp"
#include "base/grid/type/ModLinearGrid.hpp"

#include <vector>
#include <cstddef>

namespace sg
{
namespace opt
{
namespace gridgen
{

class IterativeGridGenerator
{
public:
    enum GridType
    {
        Noboundary,
        Boundary,
        Modified,
        ClenshawCurtis
    };
    
    IterativeGridGenerator(function::ObjectiveFunction &f, GridType grid_type, size_t N) :
            f(f),
            N(N),
        grid_type(grid_type),
        linear_grid(base::LinearGrid(f.getDimension())),
        linear_boundary_grid(base::LinearBoundaryGrid(f.getDimension())),
        linear_clenshawcurtis_grid(base::LinearClenshawCurtisGrid(f.getDimension())),
        mod_linear_grid(base::ModLinearGrid(f.getDimension()))
    {
        if (grid_type == GridType::Noboundary)
        {
            grid = &linear_grid;
        } else if (grid_type == GridType::Boundary)
        {
            grid = &linear_boundary_grid;
        } else if (grid_type == GridType::ClenshawCurtis)
        {
            grid = &linear_clenshawcurtis_grid;
        } else
        {
            grid = &mod_linear_grid;
        }
    }
    
    virtual ~IterativeGridGenerator() {}
    
    virtual void generate() = 0;
    
    base::Grid *getGrid() const { return grid; }
    GridType getGridType() const { return grid_type; }
    
    const std::vector<double> &getFunctionValues() const { return function_values; }
    
protected:
    function::ObjectiveFunction &f;
    size_t N;
    GridType grid_type;
    
    base::LinearGrid linear_grid;
    base::LinearBoundaryGrid linear_boundary_grid;
    base::LinearClenshawCurtisGrid linear_clenshawcurtis_grid;
    base::ModLinearGrid mod_linear_grid;
    
    base::Grid *grid;
    std::vector<double> function_values;
};

}
}
}

#endif
