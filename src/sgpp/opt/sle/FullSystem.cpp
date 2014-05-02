#include "opt/sle/FullSystem.hpp"

namespace sg
{
namespace opt
{
namespace sle
{

FullSystem::FullSystem(size_t n) :
    System(n),
    A(base::DataMatrix(n, n))
{
}

FullSystem::FullSystem(size_t n, base::DataMatrix &A, const std::vector<double> &b) :
    System(n, b),
    A(A)
{
}

FullSystem::FullSystem(size_t n, base::OperationMatrix &operation_matrix,
                       const std::vector<double> &b) :
    System(n, b),
    A(base::DataMatrix(n, n))
{
    base::DataVector x(n);
    base::DataVector y(n);
    
    x.setAll(0.0);
    
    for (size_t j = 0; j < n; j++)
    {
        if (j > 0)
        {
            x[j-1] = 0.0;
        }
        
        x[j] = 1.0;
        operation_matrix.mult(x, y);
        
        for (size_t i = 0; i < n; i++)
        {
            A.set(i, j, y[i]);
        }
    }
}

base::DataMatrix &FullSystem::getA()
{
    return A;
}

void FullSystem::setA(base::DataMatrix &A)
{
    this->A = A;
}

}
}
}
