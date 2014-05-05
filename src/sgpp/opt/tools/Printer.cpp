#include "opt/tools/Printer.hpp"
#include "base/grid/GridStorage.hpp"

#include <iostream>

namespace sg
{
namespace opt
{
namespace tools
{

Printer printer;

Printer::Printer() :
    verbose(DEFAULT_VERBOSE),
    current_level(0),
    cursor_in_clear_line(true),
    last_msg_length(0)
{
}

void Printer::printStatusBegin(const std::string &msg)
{
    if (current_level > verbose)
    {
        current_level++;
        return;
    }
    
    if (!cursor_in_clear_line)
    {
        printStatusNewLine();
    }
    
    std::cout << msg;
    printStatusNewLine();
    std::cout << "    ";
    current_level++;
    
    start_times.push(std::chrono::system_clock::now());
    last_msg_length = 0;
    cursor_in_clear_line = true;
}

void Printer::printStatusUpdate(const std::string &msg)
{
    if (current_level > verbose)
    {
        return;
    }
    
    std::cout << std::string(last_msg_length, '\b') << msg;
    
    if (last_msg_length > msg.length())
    {
        std::cout << std::string(last_msg_length - msg.length(), ' ') <<
                     std::string(last_msg_length - msg.length(), '\b');
    }
    
    std::cout << std::flush;
    last_msg_length = msg.length();
    cursor_in_clear_line = false;
}

void Printer::printStatusNewLine()
{
    if (current_level > verbose)
    {
        return;
    }
    
    std::cout << std::endl;
    
    for (size_t i = 0; i < current_level; i++)
    {
        std::cout << "    ";
    }
    
    last_msg_length = 0;
    cursor_in_clear_line = true;
}

void Printer::printStatusEnd(const std::string &msg)
{
    current_level--;
    last_duration = std::chrono::system_clock::now() - start_times.top();
    
    if (current_level > verbose)
    {
        return;
    }
    
    if (cursor_in_clear_line)
    {
        std::cout << "\b\b\b\b";
    } else
    {
        printStatusNewLine();
    }
    
    std::string time_msg =
            "Done in " + std::to_string(static_cast<int>(1000.0 * last_duration.count())) + "ms";
    
    if (msg == "")
    {
        std::cout << time_msg << ".";
    } else
    {
        std::cout << time_msg << ", " << msg;
    }
    
    printStatusNewLine();
    last_msg_length = 0;
    start_times.pop();
    cursor_in_clear_line = true;
}

void Printer::printGridToFile(const std::string &filename,
                              const gridgen::IterativeGridGenerator &grid_gen) const
{
    base::GridStorage *grid_storage = grid_gen.getGrid().getStorage();
    const std::vector<double> &function_values = grid_gen.getFunctionValues();
    size_t N = grid_storage->size();
    size_t d = grid_storage->dim();
    
    std::ofstream f(filename, std::ios::out | std::ios::binary);
    
    f.write(reinterpret_cast<const char *>(&N), sizeof(N));
    f.write(reinterpret_cast<const char *>(&d), sizeof(d));
    
    for (size_t j = 0; j < N; j++)
    {
        base::GridIndex *gp = grid_storage->get(j);
        
        for (size_t t = 0; t < d; t++)
        {
            double x = gp->abs(t);
            unsigned int l = gp->getLevel(t);
            unsigned int i = gp->getIndex(t);
            
            f.write(reinterpret_cast<const char *>(&x), sizeof(x));
            f.write(reinterpret_cast<const char *>(&l), sizeof(l));
            f.write(reinterpret_cast<const char *>(&i), sizeof(i));
        }
        
        f.write(reinterpret_cast<const char *>(&function_values[j]), sizeof(double));
    }
    
    f.close();
}

void Printer::printVectorToFile(const std::string &filename, base::DataVector &x) const
{
    std::vector<double> x_vector(x.getPointer(), x.getPointer() + x.getSize());
    printMatrixToFile(filename, x_vector, x.getSize(), 1);
}

void Printer::printMatrixToFile(const std::string &filename, base::DataMatrix &A) const
{
    std::vector<double> A_vector(A.getPointer(), A.getPointer() + A.getNrows() * A.getNcols());
    printMatrixToFile(filename, A_vector, A.getNrows(), A.getNcols());
}

double Printer::getLastDurationSecs() const
{
    return last_duration.count();
}

}
}
}
