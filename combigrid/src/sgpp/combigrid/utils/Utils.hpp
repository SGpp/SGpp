/*
 * Utils.hpp
 *
 *  Created on: 05.11.2015
 *      Author: david
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <cstddef>
#include <sgpp/globaldef.hpp>
#include <vector>
#include <string>

namespace sgpp{
namespace combigrid {

template<typename T> T pow(T base, size_t exponent) {
    T result = 1;

    size_t mask = static_cast<size_t>(1) << (8 * sizeof(size_t) - 1);

    while(exponent != 0) {
        result *= result;

        if(exponent & mask) {
            result *= base;
        }

        exponent = exponent << 1;
    }

    return result;
}

long long binom(long long n, long long k);

std::vector<std::string> split(std::string str, std::string separator);
std::string join(std::vector<std::string> const &elements, std::string const &separator);
std::string escape(std::string str, char escapeCharacter, std::string avoidCharacters, std::string replaceCharacters);
std::string unescape(std::string str, char escapeCharacter, std::string avoidCharacters, std::string replaceCharacters);

}
} /* namespace utils */

#endif /* UTILS_HPP_ */
