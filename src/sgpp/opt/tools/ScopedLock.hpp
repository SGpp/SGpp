#ifndef SGPP_OPT_TOOLS_SCOPEDLOCK_HPP
#define SGPP_OPT_TOOLS_SCOPEDLOCK_HPP

#include "opt/tools/MutexType.hpp"

namespace sg
{
namespace opt
{
namespace tools
{

struct ScopedLock
{
    ScopedLock(MutexType &m) : mut(m), locked(true) { mut.lock(); }
    ~ScopedLock() { unlock(); }
    
    void unlock() { if (!locked) return; locked = false; mut.unlock(); }
    void lockAgain() { if (locked) return; mut.lock(); locked = true; }
private:
    MutexType& mut;
    bool locked;
    
    // prevent copying the scoped lock
    void operator=(const ScopedLock &other);
    ScopedLock(const ScopedLock &other);
};

}
}
}

#endif
