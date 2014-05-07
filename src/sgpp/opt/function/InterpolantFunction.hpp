#ifndef SGPP_OPT_FUNCTION_INTERPOLANTFUNCTION_HPP
#define SGPP_OPT_FUNCTION_INTERPOLANTFUNCTION_HPP

#include "opt/function/ObjectiveFunction.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/operation/OperationEval.hpp"
//#include "base/grid/Grid.hpp"

#include "base/basis/basis.hpp"

/*#include "base/grid/type/BsplineGrid.hpp"
#include "base/grid/type/BsplineBoundaryGrid.hpp"
#include "base/grid/type/BsplineClenshawCurtisGrid.hpp"
#include "base/grid/type/ModBsplineGrid.hpp"*/

#include <vector>
#include <cstring>

namespace sg
{
namespace opt
{
namespace function
{

class InterpolantFunction : public ObjectiveFunction
{
public:
    InterpolantFunction(size_t d, base::OperationEval *op_eval, base::DataVector &alpha) :
        ObjectiveFunction(d),
        op_eval(op_eval),
        /*grid(NULL),
        grid_storage(NULL),
        basis(NULL),*/
        alpha(alpha)
    {
    }
    
    /*InterpolantFunction(size_t d, base::Grid *grid, base::DataVector &alpha) :
        ObjectiveFunction(d),
        op_eval(NULL),
        grid(grid),
        grid_storage(grid->getStorage()),
        basis(NULL),
        alpha(alpha)
    {
        if (strcmp(grid->getType(), "Bspline") == 0)
        {
            basis = new sg::base::SBsplineBase(
                    ((base::BsplineGrid *)grid)->getDegree());
        } else if (strcmp(grid->getType(), "BsplineBoundary") == 0)
        {
            basis = new sg::base::SBsplineBoundaryBase(
                    ((base::BsplineBoundaryGrid *)grid)->getDegree());
        } else if (strcmp(grid->getType(), "BsplineClenshawCurtis") == 0)
        {
            basis = new sg::base::SBsplineClenshawCurtisBase(
                    ((base::BsplineClenshawCurtisGrid *)grid)->getDegree(),
                    ((base::BsplineClenshawCurtisGrid *)grid)->getCosineTable());
        } else if (strcmp(grid->getType(), "modBspline") == 0)
        {
            basis = new sg::base::SModBsplineBase(
                    ((base::ModBsplineGrid *)grid)->getDegree());
        } else if (strcmp(grid->getType(), "Wavelet") == 0)
        {
            basis = new sg::base::SWaveletBase();
        } else if (strcmp(grid->getType(), "WaveletBoundary") == 0)
        {
            basis = new sg::base::SWaveletBoundaryBase();
        } else if (strcmp(grid->getType(), "modWavelet") == 0)
        {
            basis = new sg::base::SModWaveletBase();
        } else
        {
            throw std::invalid_argument("Grid type not supported.");
        }
    }*/
    
    /*~InterpolantFunction()
    {
        if (basis != NULL)
        {
            delete basis;
        }
    }*/
    
    inline double eval(const std::vector<double> &x)
    {
        /*if (grid == NULL)
        {*/
            std::vector<double> y = x;
            return op_eval->eval(alpha, y);
        /*} else
        {
            double result = 0.0;
            
            for (size_t i = 0; i < grid_storage->size(); i++)
            {
                const base::GridIndex *gp = grid_storage->get(i);
                double cur_result = alpha[i];
                
                for (size_t t = 0; t < x.size(); t++)
                {
                    double cur_result1d = basis->eval(gp->getLevel(t), gp->getIndex(t), x[t]);
                    
                    if (cur_result1d == 0.0)
                    {
                        cur_result = 0.0;
                        break;
                    }
                    
                    cur_result *= cur_result1d;
                }
                
                result += cur_result;
            }
            
            return result;
        }*/
    }
    
protected:
    base::OperationEval *op_eval;
    
    /*base::Grid *grid;
    base::GridStorage *grid_storage;
    base::SBase *basis;*/
    
    base::DataVector &alpha;
};

}
}
}

#endif
