#ifndef __LOGSUMSET_HPP
#define __LOGSUMSET_HPP

#include <cassert>
#include <set>

#include "logsum.h"
//#include "logger.hpp"

struct logsumset
{
    template < typename Float_Type >
    Float_Type operator () (std::multiset< Float_Type >& s) const
    {
        assert(not s.empty());
        while (s.size() > 1)
        {
            Float_Type a = *s.begin();
            s.erase(s.begin());
            Float_Type b = *s.begin();
            s.erase(s.begin());
            assert(b <= 0);
            assert(a <= b);
#ifdef LOG
            if (b - a > 16.0 and b > -80)
            {
                LOG("logsumset", warning)
                    << "precision loss: a=" << a << " b=" << b << std::endl;
            }
#endif
            s.insert(p7_FLogsum(a, b));
        }
        Float_Type res = *s.begin();
        s.erase(s.begin());
        return res;
    }
}; // struct logsumset

#endif
