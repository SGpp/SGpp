#ifndef SGPP_OPT_SLE_SYSTEM_FULL_HPP
#define SGPP_OPT_SLE_SYSTEM_FULL_HPP

#include "opt/sle/system/Cloneable.hpp"
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
namespace system
{

class Full : public Cloneable
{
public:
    Full(base::DataMatrix &A, const std::vector<double> &b) : Cloneable(b), A(A) {}
    virtual ~Full() {}
    
    inline bool isMatrixEntryNonZero(size_t i, size_t j) { return (A.get(i, j) != 0.0); }
    inline double getMatrixEntry(size_t i, size_t j) { return A.get(i, j); }
    
    base::DataMatrix &getA() { return A; }
    void setA(base::DataMatrix &A) { this->A = A; }
    
    virtual std::unique_ptr<Cloneable> clone() {
        return std::unique_ptr<Cloneable>(new Full(A, b));
    }
    
protected:
    base::DataMatrix &A;
};

}
}
}
}

#endif
