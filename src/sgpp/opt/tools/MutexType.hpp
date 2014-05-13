#ifndef SGPP_OPT_TOOLS_MUTEXTYPE_HPP
#define SGPP_OPT_TOOLS_MUTEXTYPE_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sg
{
namespace opt
{
namespace tools
{

#ifdef _OPENMP

class MutexType
{
public:
    MutexType() { omp_init_nest_lock(&_lock); }
    ~MutexType() { omp_destroy_nest_lock(&_lock); }
    void lock() { omp_set_nest_lock(&_lock); }
    void unlock() { omp_unset_nest_lock(&_lock); }
    
    MutexType(const MutexType &other) { (void)other; omp_init_nest_lock(&_lock); }
    MutexType &operator=(const MutexType &other) { (void)other; return *this; }
    
protected:
    omp_nest_lock_t _lock;
};

#else

class MutexType
{
public:
    void lock() {}
    void unlock() {}
};

#endif

}
}
}

#endif
