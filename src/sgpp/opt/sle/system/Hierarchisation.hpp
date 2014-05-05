#ifndef SGPP_OPT_SLE_SYSTEM_HIERARCHISATION_HPP
#define SGPP_OPT_SLE_SYSTEM_HIERARCHISATION_HPP

#include "opt/sle/system/Parallelizable.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/algorithm/GetAffectedBasisFunctions.hpp"

#include "base/basis/basis.hpp"
/*#include "base/basis/bspline/noboundary/BsplineBasis.hpp"
#include "base/basis/bspline/boundary/BsplineBoundaryBasis.hpp"
#include "base/basis/bspline/clenshawcurtis/BsplineClenshawCurtisBasis.hpp"
#include "base/basis/modbspline/ModifiedBsplineBasis.hpp"
#include "base/basis/wavelet/noboundary/WaveletBasis.hpp"
#include "base/basis/wavelet/boundary/WaveletBoundaryBasis.hpp"
#include "base/basis/modwavelet/ModifiedWaveletBasis.hpp"*/

#include "base/grid/type/BsplineGrid.hpp"
#include "base/grid/type/BsplineBoundaryGrid.hpp"
#include "base/grid/type/BsplineClenshawCurtisGrid.hpp"
#include "base/grid/type/ModBsplineGrid.hpp"

#include <vector>
#include <cstddef>
#include <cstring>
#include <stdexcept>
//#include <boost/functional/hash.hpp>
//#include <unordered_map>

namespace sg
{
namespace opt
{
namespace sle
{
namespace system
{

class Hierarchisation : public Parallelizable
{
public:
    Hierarchisation(base::Grid *grid, const std::vector<double> &function_values) :
        Parallelizable(function_values),
        grid(grid),
        grid_storage(grid->getStorage()),
        basis(NULL)
        /*cached_row_index(0),
        cached_row({}),
        row_cached(false)*/
        //hash_map(std::unordered_map<size_t, double>())
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
    }
    
    virtual ~Hierarchisation() { if (basis != NULL) delete basis; }
    
    inline bool isMatrixEntryNonZero(size_t i, size_t j)
            { return (evalBasisFunctionAtGridPoint(j, i) != 0.0); }
    inline double getMatrixEntry(size_t i, size_t j)
            { return evalBasisFunctionAtGridPoint(j, i); }
    
    base::Grid *getGrid() { return grid; }
    void setGrid(base::Grid *grid) { this->grid = grid; /*row_cached = false;*/ }
    
    virtual Parallelizable *clone() { return new Hierarchisation(grid, b); }
    
protected:
    base::Grid *grid;
    base::GridStorage *grid_storage;
    base::SBase *basis;
    
    /*size_t cached_row_index;
    std::vector<double> cached_row;
    bool row_cached;*/
    
    //std::unordered_map<size_t, double> hash_map;
    
    // TODO: currently unused
    /*inline void calculateRow(size_t i)
    {
        typedef std::vector<std::pair<size_t, double> > IndexValVector;
        IndexValVector vec;
        base::GetAffectedBasisFunctions<BASIS> ga(grid_storage);
        base::GridIndex *gp = grid_storage->get(i);
        size_t d = grid_storage->dim();
        size_t N = grid_storage->size();
        std::vector<double> point(d, 0.0);
        
        for (size_t t = 0; t < d; t++)
        {
            point[t] = gp->abs(t);
        }
        
        ga(basis, point, vec);
        cached_row.assign(N, 0.0);
        
        for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
        {
            cached_row[iter->first] = iter->second;
            //std::cout << "cached_row[" << iter->first << "] = " << iter->second << "\n";
        }
    }*/
    
    inline double evalBasisFunctionAtGridPoint(size_t basis_i, size_t point_j)
    {
        // ridiculously slow due to rehashing
        /*if (!row_cached || (cached_row_index != point_j))
        {
            calculateRow(point_j);
            cached_row_index = point_j;
            row_cached = true;
        }
        
        return cached_row[basis_i];*/
        
        // ridiculously slow due to hashing
        /*const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            size_t hash = 0;
            boost::hash_combine(hash, gp_basis->getLevel(t));
            boost::hash_combine(hash, gp_basis->getIndex(t));
            boost::hash_combine(hash, gp_point->getLevel(t));
            boost::hash_combine(hash, gp_point->getIndex(t));
            
            double result1d;
            
            if (hash_map.find(hash) == hash_map.end())
            {
                result1d = basis->eval(gp_basis->getLevel(t), gp_basis->getIndex(t),
                                      gp_point->abs(t));
                hash_map[hash] = result1d;
            } else
            {
                result1d = hash_map[hash];
            }
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;*/
        
        const base::GridIndex *gp_basis = grid_storage->get(basis_i);
        const base::GridIndex *gp_point = grid_storage->get(point_j);
        double result = 1.0;
        
        for (size_t t = 0; t < grid_storage->dim(); t++)
        {
            double result1d = basis->eval(gp_basis->getLevel(t), gp_basis->getIndex(t),
                                         gp_point->abs(t));
            
            if (result1d == 0.0)
            {
                return 0.0;
            }
            
            result *= result1d;
        }
        
        return result;
    }
};

/*template class HierarchisationSystem<base::SBsplineBase>;
template class HierarchisationSystem<base::SBsplineBoundaryBase>;
template class HierarchisationSystem<base::SBsplineClenshawCurtisBase>;
template class HierarchisationSystem<base::SModBsplineBase>;
template class HierarchisationSystem<base::SWaveletBase>;
template class HierarchisationSystem<base::SWaveletBoundaryBase>;
template class HierarchisationSystem<base::SModWaveletBase>;*/

}
}
}
}

#endif
