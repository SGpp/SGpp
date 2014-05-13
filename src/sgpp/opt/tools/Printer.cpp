#include "opt/tools/Printer.hpp"
#include "base/grid/GridStorage.hpp"
#include "opt/tools/ScopedLock.hpp"

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
    status_printing_enabled(true),
    current_level(0),
    cursor_in_clear_line(true),
    last_msg_length(0)
{
}

void Printer::printStatusBegin(const std::string &msg)
{
    ScopedLock lock(mutex);
    
    if (!status_printing_enabled || (current_level > verbose))
    {
        return;
    }
    
    if (!cursor_in_clear_line)
    {
        printStatusNewLine();
    }
    
    printStatusIdentation();
    current_level++;
    
    std::cout << msg;
    printStatusNewLine();
    
    start_times.push(std::chrono::system_clock::now());
    last_msg_length = 0;
    cursor_in_clear_line = true;
}

void Printer::printStatusUpdate(const std::string &msg)
{
    ScopedLock lock(mutex);
    
    if (!status_printing_enabled || (current_level > verbose))
    {
        return;
    }
    
    if (cursor_in_clear_line)
    {
        printStatusIdentation();
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
    if (!status_printing_enabled || (current_level > verbose))
    {
        return;
    }
    
    std::cout << std::endl;
    last_msg_length = 0;
    cursor_in_clear_line = true;
}

void Printer::printStatusIdentation()
{
    if (!status_printing_enabled || (current_level > verbose))
    {
        return;
    }
    
    for (size_t i = 0; i < current_level; i++)
    {
        std::cout << "    ";
    }
}

void Printer::printStatusEnd(const std::string &msg)
{
    ScopedLock lock(mutex);
    
    if (!status_printing_enabled || (current_level > verbose))
    {
        return;
    }
    
    if (!cursor_in_clear_line)
    {
        printStatusNewLine();
    }
    
    current_level--;
    printStatusIdentation();
    
    last_duration = std::chrono::system_clock::now() - start_times.top();
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

void Printer::enableStatusPrinting()
{
    status_printing_enabled = true;
}

void Printer::disableStatusPrinting()
{
    status_printing_enabled = false;
}

/*void Printer::increaseCurrentLevel(size_t inc)
{
    current_level += inc;
}

void Printer::decreaseCurrentLevel(size_t dec)
{
    current_level -= dec;
}*/

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

void Printer::printIterativeGridGenerator(const gridgen::IterativeGridGenerator &grid_gen) const
{
    base::GridStorage *grid_storage = grid_gen.getGrid().getStorage();
    const std::vector<double> &function_values = grid_gen.getFunctionValues();
    
    for (size_t i = 0; i < grid_storage->size(); i++)
    {
        if (i > 0)
        {
            std::cout << "\n";
        }
        
        base::GridIndex *gp = grid_storage->get(i);
        std::cout << gp << ", " << function_values[i];
    }
}

void Printer::printSLE(sle::system::System &system) const
{
    size_t n = system.getDimension();
    const std::vector<double> &b = system.getRHS();
    
    std::cout << "A = [";
    
    for (size_t i = 0; i < n; i++)
    {
        if (i > 0)
        {
            std::cout << ";\n";
        }
        
        for (size_t j = 0; j < n; j++)
        {
            if (j > 0)
            {
                std::cout << ", ";
            }
            
            std::cout << system.getMatrixEntry(i, j);
        }
    }
    
    std::cout << "],\nb = " << b << "]";
}

double Printer::getLastDurationSecs() const
{
    return last_duration.count();
}

MutexType &Printer::getMutex()
{
    return mutex;
}

}
}
}
