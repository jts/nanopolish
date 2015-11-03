#ifndef __LOGSUMSET_HPP
#define __LOGSUMSET_HPP

#include <cassert>
#include <cmath>
#include <set>

#include "logsum.h"
//#include "logger.hpp"

template< typename Float_Type >
class logsumset
{
public:
    logsumset(bool use_set) : _val_set(), _val(-INFINITY), _use_set(use_set) {}

    template < typename Input_Iterator >
    logsumset(Input_Iterator it, Input_Iterator it_end, bool use_set)
        : logsumset(use_set)
    {
        while (it != it_end)
        {
            add(*it);
        }
    }

    template < typename Input_Range >
    logsumset(const Input_Range& rg, bool use_set)
        : logsumset(rg.begin(), rg.end(), use_set) {}

    const bool& use_set() const { return _use_set; }
    bool& use_set() { return _use_set; }

    void add(Float_Type v)
    {
        if (_use_set)
        {
            _val_set.insert(v);
        }
        else
        {
            _val = p7_FLogsum(_val, v);
        }
    }

    Float_Type val()
    {
        if (not _val_set.empty())
        {
            _val_set.insert(_val);
            while (_val_set.size() > 1)
            {
                Float_Type a = *_val_set.begin();
                assert(not std::isnan(a));
                _val_set.erase(_val_set.begin());
                Float_Type b = *_val_set.begin();
                assert(not std::isnan(b));
                _val_set.erase(_val_set.begin());
#ifdef LOG
                if (b - a > 15.7 and b > -80)
                {
                    LOG("logsumset", warning)
                        << "precision loss: a=" << a << " b=" << b << std::endl;
                }
#endif
                _val_set.insert(p7_FLogsum(a, b));
            }
            _val = *_val_set.begin();
            _val_set.erase(_val_set.begin());
        }
        return _val;
    }

private:
    std::multiset< Float_Type > _val_set;
    Float_Type _val;
    bool _use_set;
}; // class logsumset

#endif
