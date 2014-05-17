#ifndef SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP
#define SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP

#include "opt/gridgen/IterativeGridGenerator.hpp"
#include "base/basis/Basis.hpp"

#include <cstddef>

namespace sg
{
namespace opt
{
namespace gridgen
{

class IterativeGridGeneratorLinearSurplus : public IterativeGridGenerator
{
public:
    static const size_t INITIAL_LEVEL_WO_BOUNDARIES = 3;
    
    IterativeGridGeneratorLinearSurplus(function::Objective &f, base::Grid &grid,
                                        size_t N, double alpha);
    
    bool generate();
    
    double getAlpha() const;
    void setAlpha(double alpha);
    
protected:
    std::unique_ptr<base::SBase> linear_base;
    std::unique_ptr<base::Grid> linear_grid;
    
    /*base::SLinearBase linear_base;
    base::SLinearBoundaryBase linear_boundary_base;
    base::SLinearClenshawCurtisBase linear_clenshawcurtis_base;
    base::SModLinearBase mod_linear_base;*/
    
    double alpha;
    
    double evalBasisFunctionAtGridPoint(base::GridStorage *grid_storage, size_t i, size_t point_j);
};

}
}
}

#endif
