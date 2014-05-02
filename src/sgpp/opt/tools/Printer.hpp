#ifndef SGPP_OPT_TOOLS_PRINTER_HPP
#define SGPP_OPT_TOOLS_PRINTER_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "opt/gridgen/IterativeGridGenerator.hpp"

#include <cstddef>
#include <chrono>
#include <stack>
#include <type_traits>

namespace sg
{
namespace opt
{
namespace tools
{

class Printer
{
public:
    static const size_t DEFAULT_VERBOSE = 0;
    
    Printer();
    
    void printStatusBegin(const std::string &msg);
    void printStatusUpdate(const std::string &msg);
    void printStatusNewLine();
    void printStatusEnd(const std::string &msg = "");
    
    void printGridToFile(const std::string &filename,
                         const gridgen::IterativeGridGenerator &grid_gen) const;
    
    void printVectorToFile(const std::string &filename, base::DataVector &x) const;
    template <class T>
    void printVectorToFile(const std::string &filename, std::vector<T> &x) const
    {
        printMatrixToFile(filename, x, x.size(), 1);
    }
    
    void printMatrixToFile(const std::string &filename, base::DataMatrix &A) const;
    template <class T>
    void printMatrixToFile(const std::string &filename, std::vector<T> &A,
                           size_t m, size_t n) const
    {
        std::ofstream f(filename, std::ios::out | std::ios::binary);
        
        std::vector<const char *> types =
                {"uint8           ", "uint16          ", "uint32          ", "uint64          ",
                 "double          ", "other           "};
        size_t index;
        
        if (std::is_same<T, uint8_t>::value)
        {
            index = 0;
        } else if (std::is_same<T, uint16_t>::value)
        {
            index = 1;
        } else if (std::is_same<T, uint32_t>::value)
        {
            index = 2;
        } else if (std::is_same<T, uint64_t>::value)
        {
            index = 3;
        } else if (std::is_same<T, double>::value)
        {
            index = 4;
        } else
        {
            index = 5;
        }
        
        f.write(reinterpret_cast<const char *>(&m), sizeof(m));
        f.write(reinterpret_cast<const char *>(&n), sizeof(n));
        f.write(types[index], 16);
        
        for (size_t i = 0; i < m*n; i++)
        {
            f.write(reinterpret_cast<const char *>(&A[i]), sizeof(T));
        }

        f.close();
    }
    
    inline size_t getVerbosity() const { return verbose; }
    inline void setVerbosity(size_t level) { verbose = level; }
    
    double getLastDurationSecs() const;
    
protected:
    size_t verbose;
    
    size_t current_level;
    bool cursor_in_clear_line;
    size_t last_msg_length;
    std::stack<std::chrono::time_point<std::chrono::system_clock> > start_times;
    std::chrono::duration<double> last_duration;
};

extern Printer printer;

}
}
}

#endif
