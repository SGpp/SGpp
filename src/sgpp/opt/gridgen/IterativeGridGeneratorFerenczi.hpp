#ifndef SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATORFERENCZI_HPP
#define SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATORFERENCZI_HPP

#include "opt/gridgen/IterativeGridGenerator.hpp"

#include <vector>
#include <cstddef>

namespace sg
{
namespace opt
{
namespace gridgen
{

class IterativeGridGeneratorFerenczi : public IterativeGridGenerator
{
public:
    static const size_t INITIAL_LEVEL = 1;
    static const size_t MAX_LEVEL = 25;
    
    IterativeGridGeneratorFerenczi(function::ObjectiveFunction &f, GridType grid_type,
                                   size_t N, double alpha);
    
    void generate();
    
    double getAlpha() const;
    void setAlpha(double alpha);
    
protected:
    double alpha;
    
    void rankVector(const std::vector<double> &vec, std::vector<size_t> &rank) const;
};

}
}
}

#endif
