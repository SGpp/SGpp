#ifndef SGPP_OPT_SLE_HIERARCHISATIONSYSTEM_HPP
#define SGPP_OPT_SLE_HIERARCHISATIONSYSTEM_HPP

#include "opt/sle/System.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/algorithm/GetAffectedBasisFunctions.hpp"

#include <vector>
#include <cstddef>
//#include <boost/functional/hash.hpp>
//#include <unordered_map>

namespace sg
{
namespace opt
{
namespace sle
{

template <class BASIS>
class HierarchisationSystem : public System
{
public:
    HierarchisationSystem(base::Grid *grid, BASIS &basis,
                          const std::vector<double> &function_values);
    
    inline bool isMatrixEntryNonZero(size_t i, size_t j)
    {
        return (evalBasisFunctionAtGridPoint(j, i) != 0.0);
    }
    
    inline double getMatrixEntry(size_t i, size_t j)
    {
        return evalBasisFunctionAtGridPoint(j, i);
    }
    
    base::Grid *getGrid();
    void setGrid(base::Grid *grid);
    
protected:
    base::Grid *grid;
    base::GridStorage *grid_storage;
    BASIS &basis;
    
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
                result1d = basis.eval(gp_basis->getLevel(t), gp_basis->getIndex(t),
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
            double result1d = basis.eval(gp_basis->getLevel(t), gp_basis->getIndex(t),
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

}
}
}

#endif
