#ifndef SGPP_OPT_SLE_FULLSYSTEM_HPP
#define SGPP_OPT_SLE_FULLSYSTEM_HPP

#include "opt/sle/System.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/operation/OperationMatrix.hpp"

#include <vector>
#include <cstddef>

namespace sg
{
namespace opt
{
namespace sle
{

class FullSystem : public System
{
public:
    FullSystem(size_t n);
    FullSystem(size_t n, base::DataMatrix &A, const std::vector<double> &b);
    FullSystem(size_t n, base::OperationMatrix &operation_matrix,
               const std::vector<double> &b);
    
    inline bool isMatrixEntryNonZero(size_t i, size_t j) { return (A.get(i, j) != 0.0); }
    inline double getMatrixEntry(size_t i, size_t j) { return A.get(i, j); }
    
    base::DataMatrix &getA();
    void setA(base::DataMatrix &A);
    
protected:
    base::DataMatrix A;
};

}
}
}

#endif
