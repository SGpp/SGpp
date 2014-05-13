#ifndef SGPP_OPT_FUNCTION_INTERPOLANT_HPP
#define SGPP_OPT_FUNCTION_INTERPOLANT_HPP

#include "opt/function/Objective.hpp"
#include "base/datatypes/DataVector.hpp"
//include "base/operation/OperationEval.hpp"
#include "base/grid/Grid.hpp"

/*#include "base/basis/Basis.hpp"
#include "base/grid/type/BsplineGrid.hpp"
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

class Interpolant : public Objective
{
public:
    /*Interpolant(size_t d, base::OperationEval *op_eval, base::DataVector &alpha) :
        Objective(d),
        op_eval(op_eval),
        alpha(alpha)
    {
    }*/
    
    Interpolant(size_t d, base::Grid &grid, base::DataVector &alpha) :
        Objective(d),
        grid(grid),
        op_eval(op_factory::createOperationEval(grid)),
        alpha(alpha)
    {
    }
    
    /*Interpolant(size_t d, base::OperationEval *op_eval, base::DataVector &alpha) :
        Objective(d),
        op_eval(op_eval),
        grid(nullptr),
        grid_storage(nullptr),
        basis(nullptr),
        alpha(alpha)
    {
    }
    
    Interpolant(size_t d, base::Grid *grid, base::DataVector &alpha) :
        ObjectiveFunction(d),
        op_eval(nullptr),
        grid(grid),
        grid_storage(grid->getStorage()),
        basis(nullptr),
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
    
    inline double eval(const std::vector<double> &x)
    {
        /*if (grid == nullptr)
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
    
    virtual std::unique_ptr<Objective> clone()
    {
        return std::unique_ptr<Objective>(new Interpolant(d, grid, alpha));
    }
    
protected:
    base::Grid &grid;
    std::unique_ptr<base::OperationEval> op_eval;
    
    /*base::Grid *grid;
    base::GridStorage *grid_storage;
    base::SBase *basis;*/
    
    base::DataVector &alpha;
};

}
}
}

#endif
