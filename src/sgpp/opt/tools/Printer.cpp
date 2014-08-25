/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

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
    verbose(DEFAULT_VERBOSITY),
    status_printing_enabled(true),
    status_level(0),
    indentation_level(0),
    cursor_in_clear_line(true),
    last_msg_length(0)
{
}

void Printer::printStatusBegin(const std::string &msg)
{
    ScopedLock lock(mutex);
    
    if (!status_printing_enabled || (status_level > verbose))
    {
        // status printing disabled or verbose level too low
        status_level++;
        return;
    }
    
    // go to new line
    if (!cursor_in_clear_line)
    {
        printStatusNewLine();
    }
    
    // print indentation
    printStatusIdentation();
    
    // print message
    std::cout << msg;
    printStatusNewLine();
    
    // timing
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    start_times.push(ts);
    
    last_msg_length = 0;
    cursor_in_clear_line = true;
    status_level++;
    indentation_level++;
}

void Printer::printStatusUpdate(const std::string &msg)
{
    ScopedLock lock(mutex);
    
    if (!status_printing_enabled || (status_level > verbose))
    {
        // status printing disabled or verbose level too low
        return;
    }
    
    // print indentation
    if (cursor_in_clear_line)
    {
        printStatusIdentation();
    }
    
    // go back to start of last message and overwrite it with the new message
    std::cout << std::string(last_msg_length, '\b') << msg;
    
    // new message is too short ==> print spaces to hide the old one and go back again
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
    if (!status_printing_enabled || (status_level > verbose))
    {
        // status printing disabled or verbose level too low
        return;
    }
    
    std::cout << std::endl;
    last_msg_length = 0;
    cursor_in_clear_line = true;
}

void Printer::printStatusIdentation()
{
    if (!status_printing_enabled || (status_level > verbose))
    {
        // status printing disabled or verbose level too low
        return;
    }
    
    // print indentation
    for (int i = 0; i < indentation_level; i++)
    {
        std::cout << "    ";
    }
}

void Printer::printStatusEnd(const std::string &msg)
{
    ScopedLock lock(mutex);
    
    status_level--;
    
    if (!status_printing_enabled || (status_level > verbose))
    {
        // status printing disabled or verbose level too low
        return;
    }
    
    // go to new line
    if (!cursor_in_clear_line)
    {
        printStatusNewLine();
    }
    
    // print indentation
    indentation_level--;
    printStatusIdentation();
    
    // timing
    timespec ts1 = start_times.top(), ts2;
    clock_gettime(CLOCK_MONOTONIC, &ts2);
    last_duration = static_cast<double>(ts2.tv_sec - ts1.tv_sec) +
            static_cast<double>(ts2.tv_nsec - ts1.tv_nsec) / 1e9;
    
    // print message
    std::string time_msg =
            "Done in " + toString(static_cast<size_t>(1000.0 * last_duration)) + "ms";
    
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

void Printer::printGridToFile(const std::string &filename,
                              const gridgen::IterativeGridGenerator &grid_gen) const
{
    base::GridStorage *grid_storage = grid_gen.getGrid().getStorage();
    const std::vector<double> &function_values = grid_gen.getFunctionValues();
    size_t N = grid_storage->size();
    size_t d = grid_storage->dim();
    
    std::ofstream f(filename.c_str(), std::ios::out | std::ios::binary);
    
    // header (number of points and dimensions)
    f.write(reinterpret_cast<const char *>(&N), sizeof(N));
    f.write(reinterpret_cast<const char *>(&d), sizeof(d));
    
    for (size_t j = 0; j < N; j++)
    {
        base::GridIndex *gp = grid_storage->get(j);
        
        for (size_t t = 0; t < d; t++)
        {
            double x = gp->abs(t);
            base::GridIndex::level_type l = gp->getLevel(t);
            base::GridIndex::index_type i = gp->getIndex(t);
            
            // coordinate, level and index of current grid point
            f.write(reinterpret_cast<const char *>(&x), sizeof(x));
            f.write(reinterpret_cast<const char *>(&l), sizeof(l));
            f.write(reinterpret_cast<const char *>(&i), sizeof(i));
        }
        
        // function value at the current grid point
        f.write(reinterpret_cast<const char *>(&function_values[j]), sizeof(double));
    }
    
    f.close();
}

void Printer::printVectorToFile(const std::string &filename, base::DataVector &x) const
{
    // convert DataVector to std::vector
    std::vector<double> x_vector(x.getPointer(), x.getPointer() + x.getSize());
    printMatrixToFile(filename, x_vector, x.getSize(), 1);
}

void Printer::printMatrixToFile(const std::string &filename, base::DataMatrix &A) const
{
    // convert DataMatrix to std::vector
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
        
        // print grid point and function value
        std::cout << grid_storage->get(i) << ", " << function_values[i];
    }
}

void Printer::printSLE(sle::system::System &system) const
{
    const size_t n = system.getDimension();
    
    std::cout << "A = [";
    
    // print matrix
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
    
    std::cout << "]";
}

double Printer::getLastDurationSecs() const
{
    return last_duration;
}

MutexType &Printer::getMutex()
{
    return mutex;
}

template <>
const char *Printer::getTypeString(const std::vector<uint8_t> &A) const
{
    (void)A;
    return "uint8           ";
}

template <>
const char *Printer::getTypeString(const std::vector<uint16_t> &A) const
{
    (void)A;
    return "uint16          ";
}

template <>
const char *Printer::getTypeString(const std::vector<uint32_t> &A) const
{
    (void)A;
    return "uint32          ";
}

template <>
const char *Printer::getTypeString(const std::vector<uint64_t> &A) const
{
    (void)A;
    return "uint64          ";
}

template <>
const char *Printer::getTypeString(const std::vector<double> &A) const
{
    (void)A;
    return "double          ";
}

template <>
const char *Printer::getTypeString(const std::vector<std::string> &A) const
{
    (void)A;
    return "string          ";
}

template <>
void Printer::writeEntryToFile(std::ofstream &f, const std::string &entry) const
{
    // write string terminated by null character
    const char null_char[1] = {'\0'};
    f << entry;
    f.write(null_char, 1);
}

}
}
}
