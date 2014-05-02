#ifndef SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP
#define SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP

#include "opt/gridgen/IterativeGridGenerator.hpp"
#include "base/basis/linear/noboundary/LinearBasis.hpp"
#include "base/basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "base/basis/linear/clenshawcurtis/LinearClenshawCurtisBasis.hpp"
#include "base/basis/modlinear/ModifiedLinearBasis.hpp"

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
    static const size_t INITIAL_LEVEL = 1;
    
    IterativeGridGeneratorLinearSurplus(function::ObjectiveFunction &f, GridType grid_type,
                                        size_t N, double alpha);
    IterativeGridGeneratorLinearSurplus(function::ObjectiveFunction &f, GridType grid_type,
            size_t N, double alpha, const base::CosineTable *cosine_table);
    
    void generate();
    
    double getAlpha() const;
    void setAlpha(double alpha);
    
protected:
    base::SLinearBase linear_base;
    base::SLinearBoundaryBase linear_boundary_base;
    base::SLinearClenshawCurtisBase linear_clenshawcurtis_base;
    base::SModLinearBase mod_linear_base;
    
    double alpha;
    
    double evalBasisFunctionAtGridPoint(base::GridStorage *grid_storage, size_t i, size_t point_j);
};

}
}
}

#endif
