
#ifndef HASH_MAP_CONFIG
#define HASH_MAP_CONFIG

#include <ext/hash_map>

namespace std {

    using namespace __gnu_cxx;

}


namespace sg
{

template<class key>
struct hash { };

template<class key>
struct eqIndex { };

}
#endif
