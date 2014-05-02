#ifndef SGPP_OPT_SLE_SYSTEM_HPP
#define SGPP_OPT_SLE_SYSTEM_HPP

#include <vector>
#include <cstddef>

namespace sg
{
namespace opt
{
namespace sle
{

class System
{
public:
    System(size_t n) : System(n, std::vector<double>(n, 0.0)) {}
    System(size_t n, const std::vector<double> &b) : n(n), b(b) {}
    virtual ~System() {}
    
    virtual bool isMatrixEntryNonZero(size_t i, size_t j) = 0;
    virtual double getMatrixEntry(size_t i, size_t j) = 0;
    
    virtual void matrixVectorMultiplication(const std::vector<double> &x,
                                            std::vector<double> &y)
    {
        y = std::vector<double>(n, 0.0);
        
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                y[i] += getMatrixEntry(i, j) * x[j];
            }
        }
    }
    
    virtual size_t countNNZ()
    {
        size_t nnz = 0;
        
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                if (isMatrixEntryNonZero(i, j))
                {
                    nnz++;
                }
            }
        }
        
        return nnz;
    }
    
    const std::vector<double> &getRHS() const { return b; }
    void setRHS(const std::vector<double> &b) { this->b = b; }
    /*void setRHS(const std::vector<double> &b) {
        if (b.size() != n)
        {
            std::cerr << "sg::opt::sle::System::setRHS: Wrong size of b!\n";
            return;
        }
        
        this->b = b;
    }*/
    
    size_t getDimension() const { return n; }
    
protected:
    size_t n;
    std::vector<double> b;
};

}
}
}

#endif
