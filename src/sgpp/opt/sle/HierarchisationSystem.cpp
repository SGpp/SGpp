#include "opt/sle/HierarchisationSystem.hpp"
#include "base/basis/bspline/noboundary/BsplineBasis.hpp"
#include "base/basis/bspline/boundary/BsplineBoundaryBasis.hpp"
#include "base/basis/bspline/clenshawcurtis/BsplineClenshawCurtisBasis.hpp"
#include "base/basis/modbspline/ModifiedBsplineBasis.hpp"
#include "base/basis/wavelet/noboundary/WaveletBasis.hpp"
#include "base/basis/wavelet/boundary/WaveletBoundaryBasis.hpp"
#include "base/basis/modwavelet/ModifiedWaveletBasis.hpp"

namespace sg
{
namespace opt
{
namespace sle
{

template <class BASIS>
HierarchisationSystem<BASIS>::HierarchisationSystem(
        base::Grid *grid, BASIS &basis,
        const std::vector<double> &function_values) :
    System(grid->getStorage()->size(), function_values),
    grid(grid),
    grid_storage(grid->getStorage()),
    basis(basis)
    /*cached_row_index(0),
    cached_row({}),
    row_cached(false)*/
    //hash_map(std::unordered_map<size_t, double>())
{
}

template <class BASIS>
base::Grid *HierarchisationSystem<BASIS>::getGrid()
{
    return grid;
}

template <class BASIS>
void HierarchisationSystem<BASIS>::setGrid(base::Grid *grid)
{
    this->grid = grid;
    //row_cached = false;
}

template class HierarchisationSystem<base::SBsplineBase>;
template class HierarchisationSystem<base::SBsplineBoundaryBase>;
template class HierarchisationSystem<base::SBsplineClenshawCurtisBase>;
template class HierarchisationSystem<base::SModBsplineBase>;
template class HierarchisationSystem<base::SWaveletBase>;
template class HierarchisationSystem<base::SWaveletBoundaryBase>;
template class HierarchisationSystem<base::SModWaveletBase>;

}
}
}
