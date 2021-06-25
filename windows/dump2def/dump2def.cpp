#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>

typedef std::set<std::string> SymbolMap;

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

// Written and placed in public domain by Jeffrey Walton

std::string exec(const char* cmd) {
    std::array<char, 512> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&_pclose)> pipe(_popen(cmd, "r"), _pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    //std::cout << result << std::endl;
    return result;
}

int main(int argc, char* argv[])
{
    // ******************** Handle Options ******************** //

    std::vector<std::string> opts;
    for (size_t i = 0; i < argc; ++i)
        opts.push_back(argv[i]);

    std::string moduleName = opts[1];
    std::string demumble = opts[2];

    // ******************** Read MAP file ******************** //

    SymbolMap symbols;

    try
    {
        std::ifstream infile(opts[3].c_str());
        std::string::size_type pos;
        std::string line;

        // Find start of the symbol table
        while (std::getline(infile, line))
        {
            pos = line.find("public symbols");
            if (pos == std::string::npos) { continue; }

            // Eat the whitespace after the table heading
            infile >> std::ws;
            break;
        }

        while (std::getline(infile, line))
        {
            // End of table
            if (line.empty()) { break; }

            std::istringstream iss(line);
            std::string address, symbol;

            iss >> address >> symbol;
            // Filter out anything we don't need
            if ((symbol.find(moduleName) == std::string::npos) &&
                symbol.find("@op_factory@sgpp@") == std::string::npos &&
                symbol.find("@json@") == std::string::npos)
            {
                continue;
            }

            if (symbol.rfind("?_", 0) == 0 || symbol.rfind("??$_", 0) == 0 || symbol.find("@?$_") != std::string::npos) {
              continue;
            }

            // Demumble symbol
            std::string readable = exec((demumble + " \"" + symbol + "\"").c_str());

            // Remove deleting destructors
            if ((readable.find("deleting d") != std::string::npos)
                || readable.find(symbol) != std::string::npos) // No symbol
			{
			    continue;
			}

            symbols.insert(symbol);
        }
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Unexpected exception:" << std::endl;
        std::cerr << ex.what() << std::endl;
        std::cerr << std::endl;
    }

    // ******************** Write DEF file ******************** //

    try
    {
        std::ofstream outfile(opts[4].c_str());

        std::string name = opts[4];
        std::string::size_type pos = name.find_last_of(".");

        if (pos != std::string::npos)
            name.erase(pos);

        outfile << "LIBRARY " << name << std::endl;
        outfile << "EXPORTS" << std::endl;
        outfile << std::endl;

        outfile << "\t;; " << symbols.size() << " symbols" << std::endl;

        auto it = symbols.begin();
        for (; it != symbols.end(); ++it)
            outfile << "\t" << *it << std::endl;
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Unexpected exception:" << std::endl;
        std::cerr << ex.what() << std::endl;
        std::cerr << std::endl;
    }

    return 0;
}