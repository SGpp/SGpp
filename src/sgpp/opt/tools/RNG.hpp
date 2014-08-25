/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_TOOLS_RNG_HPP
#define SGPP_OPT_TOOLS_RNG_HPP

namespace sg
{
namespace opt
{
namespace tools
{

/**
 * Singleton class for generating pseudo-random numbers.
 * For comparability reasons, this class generates the same random numbers
 * as the C++11-only header <tt>&lt;random&gt;</tt> with the MT19937 random number generator.
 */
class RNG
{
public:
    /// type of the seed
    typedef size_t SeedType;
    
    /**
     * Constructor, initializes with zero seed.
     */
    RNG();
    
    /**
     * Generate a uniform pseudo-random number.
     * 
     * @param a lower bound
     * @param b upper bound
     * @return  uniform pseudo-random number in \f$[a, b]\f$
     */
    double getUniformRN(double a = 0.0, double b = 1.0);
    
    /**
     * Generate a uniform pseudo-random array index.
     * 
     * @param size  size of the array
     * @return      discrete uniform pseudo-random number in
     *              \f$\{0, \dotsc, \text{\texttt{size}} - 1\}\f$
     */
    size_t getUniformIndexRN(size_t size);
    
    /**
     * Generate a Gaussian pseudo-random number.
     * 
     * @param std_dev   standard deviation of the Gaussian distribution
     * @param mean      mean of the Gaussian distribution
     * @return          Gaussian pseudo-random number
     */
    double getGaussianRN(double std_dev = 1.0, double mean = 0.0);
    
    /**
     * @return      seed
     */
    SeedType getSeed();
    
    /**
     * Reseeds with time-dependent seed.
     */
    void setSeed();
    
    /**
     * Reseeds.
     * 
     * @param seed  seed to be used
     */
    void setSeed(SeedType seed);
    
protected:
    /// seed for MT19937
    SeedType seed;
    /// MT19937 state vector
    size_t mt_state[624];
    /// current MT19937 index
    size_t mt_index;
    
    /**
     * Initializes MT19937.
     */
    void initialize_mt();
    
    /**
     * Generates a new pseudo-random number via MT19937.
     * 
     * @return pseudo-random number
     */
    size_t run_mt();
};

/// singleton random number generator instance
extern RNG rng;

}
}
}

#endif
