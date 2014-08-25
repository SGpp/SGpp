/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_TOOLS_SCOPEDLOCK_HPP
#define SGPP_OPT_TOOLS_SCOPEDLOCK_HPP

#include "opt/tools/MutexType.hpp"

namespace sg
{
namespace opt
{
namespace tools
{

/**
 * Wrapper around MutexType which locks and unlocks upon construction/destruction.
 * 
 * Adopted from http://bisqwit.iki.fi/story/howto/openmp/#Locks.
 */
class ScopedLock
{
public:
    /**
     * Constructor, locks the MutexType object.
     * 
     * @param m     MutexType object to be wrapped
     */
    ScopedLock(MutexType &m) : mut(m), locked(true)
    {
        mut.lock();
    }
    
    /**
     * Destructor, unlocks the MutexType object.
     */
    ~ScopedLock()
    {
        unlock();
    }
    
    /**
     * Unlocks the MutexType object, if locked.
     */
    void unlock()
    {
        if (locked)
        {
            locked = false;
            mut.unlock();
        }
    }
    
    /**
     * Re-locks the MutexType object, if unlocked.
     */
    void lockAgain()
    {
        if (!locked)
        {
            mut.lock();
            locked = true;
        }
    }
    
protected:
    /// underlying MutexType object
    MutexType &mut;
    /// whether the MutexType object is locked or not
    bool locked;
    
private:
    /**
     * Custom copy constructor to prevent copying the lock.
     * 
     * @param other     object to be copied
     */
    ScopedLock(const ScopedLock &other);
    
    /**
     * Custom assignment operator to prevent copying the lock.
     * 
     * @param other     object to be assigned to
     */
    void operator=(const ScopedLock &other);
};

}
}
}

#endif
