#ifndef SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP
#define SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP

#include "opt/function/Objective.hpp"
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
    IterativeGridGenerator(function::Objective &f, base::Grid &grid, size_t N) :
        f(f), grid(grid), N(N) {}
    virtual ~IterativeGridGenerator() {}
    
    virtual bool generate() = 0;
    
    base::Grid &getGrid() const { return grid; }
    const std::vector<double> &getFunctionValues() const { return function_values; }
    
protected:
    function::Objective &f;
    base::Grid &grid;
    size_t N;
    
    std::vector<double> function_values;
};

}
}
}

#endif
