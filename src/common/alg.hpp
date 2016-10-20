//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Matei David (matei.david@oicr.on.ca)
// Part of https://github.com/mateidavid/hpptools
//---------------------------------------------------------
//

/// Extensions of various standard algorithms.
///
/// Contents:
///
/// for_each[|_it][|_advance]
///
///   Args: (iterator pair/range, functor)
///
///   Return: functor
///
///   Synopsis:
///     Apply functor to elements in range.
///
///   Versions:
///     _it : functor is given iterator to element
///     non-_it : functor is given element
///     _advance : save next iterator before applying function.
///                potential use: filter list in one pass
///
/// [min|max|minmax][|_value]_of
///
///   Args: (iterator pair/range, optional functor)
///
///   Return: min, max, pair of (min,max) in range
///
///   Synopsis:
///     Find minimum and/or maximum in range, using functor to extract keys.
///
///   Versions:
///     _value : return min/max values
///     non-_value : return iterators to min/max values
///
///   Note:
///     The range non-_value versions do not accept rvalue reference ranges,
///     because the iterator they return might be invalid.
///
/// mean_stdv_of
///
///   Args: (iterator pair/range, optional functor)
///
///   Synopsis:
///     Compute mean and sample stdv of a range, using functor to extract keys.
///
/// [equal|all|any]_of
///
///   Args:
///     (iterator pair/range, optional functor)
///
///   Synopsis:
///     Return true iff elements extracted with functor are all equal (equal_of),
///     or all are true (all_of), or at least one is true (any_of).
///
/// os_join
///
///   Args:
///     (iterator pair/range, separator, optional functor)
///
///   Synopsis:
///     Use operator << overloads to print the given range to stream, using the given
///     element separator, optionally passing elements through functor.

#ifndef __ALG_HPP
#define __ALG_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <type_traits>

namespace alg
{

namespace detail
{

struct Identity
{
    //template < typename U >
    //constexpr auto operator () (U&& v) const noexcept -> decltype(std::forward< U >(v))
    //{
    //    return std::forward< U >(v);
    //}
    template < typename U >
    U operator () (const U& v) const { return v; }
}; // struct Identity

/**
 * Temporary object that will cause the given range to be printed to an ostream.
 * Create such objects implicitly using the os_join() functions below.
 */
template < typename Input_Iterator, typename Unary_Function >
struct os_join_helper
{
    Input_Iterator it_begin;
    Input_Iterator it_end;
    const std::string& sep;
    Unary_Function& fn;

    os_join_helper(Input_Iterator _it_begin, Input_Iterator _it_end, const std::string& _sep, Unary_Function&& _fn)
        : it_begin(_it_begin), it_end(_it_end), sep(_sep), fn(_fn) {}

    friend std::ostream& operator << (std::ostream& os, const os_join_helper& o)
    {
        bool first = true;
        for (auto it = o.it_begin; it != o.it_end; ++it)
        {
            if (not first)
            {
                os << o.sep;
            }
            first = false;
            os << o.fn(*it);
        }
        return os;
    }

    operator std::string() const
    {
        std::ostringstream os;
        os << *this;
        return os.str();
    }
}; // struct os_join_helper

} //namespace detail

/**
 * for_each; iterator version.
 */
using std::for_each;

/**
 * for_each; range version.
 */
template < class Input_Range, class Unary_Function >
Unary_Function
for_each(Input_Range&& rg, Unary_Function&& fn)
{
    return for_each(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Apply function to elements in range, but advance iterator first; iterator version.
 * With this, it is safe to apply a function which might indirectly remove
 * the pointed element from a list.
 * @param first Start of range.
 * @param last End of range.
 * @param fn Unary function taking an element as argument.
 * @return fn
 */
template < class Input_Iterator, class Unary_Function >
Unary_Function
for_each_advance(Input_Iterator first, Input_Iterator last, Unary_Function&& fn)
{
    while (first != last)
    {
        fn(*(first++));
    }
    return fn;
}

/**
 * Apply function to elements in range, but advance iterator first; range version.
 * With this, it is safe to apply a function which might indirectly remove
 * the pointed element from a list.
 * @param rg Range.
 * @param fn Unary function taking an element as argument.
 * @return fn
 */
template < class Forward_Range, class Unary_Function >
Unary_Function
for_each_advance(Forward_Range&& rg, Unary_Function&& fn)
{
    return for_each_advance(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Apply function to iterators in range; iterator version.
 * @param first Start of range.
 * @param last End of range.
 * @param fn Unary function taking an iterator as argument.
 * @return fn
 */
template < class Input_Iterator, class Unary_Function >
Unary_Function
for_each_it(Input_Iterator first, Input_Iterator last, Unary_Function&& fn)
{
    while (first != last)
    {
        fn(first);
        ++first;
    }
    return fn;
}

/**
 * Apply function to iterators in range; range version.
 * @param first Start of range.
 * @param last End of range.
 * @param fn Unary function taking an iterator as argument.
 * @return fn
 */
template < class Forward_Range, class Unary_Function >
Unary_Function
for_each_it(Forward_Range&& rg, Unary_Function&& fn)
{
    return for_each_it(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Apply function to iterators in range, but advance iterator first; iterator version.
 * With this, it is safe to remove elements from a list on the go.
 * @param first Start of range.
 * @param last End of range.
 * @param fn Unary function taking an element as argument.
 * @return fn
 */
template < class Input_Iterator, class Unary_Function >
Unary_Function
for_each_it_advance(Input_Iterator first, Input_Iterator last, Unary_Function&& fn)
{
    while (first != last)
    {
        fn(first++);
    }
    return fn;
}

/**
 * Apply function to iterators in range, but advance iterator first; iterator version.
 * With this, it is safe to remove elements from a list on the go.
 * @param first Start of range.
 * @param last End of range.
 * @param fn Unary function taking an element as argument.
 * @return fn
 */
template < class Forward_Range, class Unary_Function >
Unary_Function
for_each_it_advance(Forward_Range&& rg, Unary_Function&& fn)
{
    return for_each_it_advance(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Return iterator to minimum element in range; iterator version.
 * @param first Range begin
 * @param last Range end
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Iterator, class Unary_Function = detail::Identity >
Forward_Iterator
min_of(Forward_Iterator first, Forward_Iterator last, Unary_Function&& fn = Unary_Function())
{
    if (first == last) return last;
    auto min_it = first++;
    for (auto it = first; it != last; ++it)
    {
        if (fn(*it) < fn(*min_it))
        {
            min_it = it;
        }
    }
    return min_it;
}

/**
 * Return iterator to minimum element in range; range version.
 * @param rg Range
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Range, class Unary_Function = detail::Identity >
auto
min_of(Forward_Range&& rg, Unary_Function&& fn = Unary_Function())
    -> decltype(rg.end())
{
    static_assert(not std::is_rvalue_reference< Forward_Range >::value,
                  "rvalue reference to range not allowed: returned iterator might not exist");
    return min_of(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Return minimum element in range; iterator version.
 * @param first Range begin
 * @param last Range end
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Iterator, class Unary_Function = detail::Identity >
typename std::result_of< Unary_Function(typename std::iterator_traits< Forward_Iterator >::value_type) >::type
min_value_of(Forward_Iterator first, Forward_Iterator last, Unary_Function&& fn = Unary_Function())
{
    auto it = min_of(first, last, std::forward< Unary_Function >(fn));
    return it != last
        ? fn(*it)
        : typename std::result_of< Unary_Function(typename std::iterator_traits< Forward_Iterator >::value_type) >::type();
}

/**
 * Return minimum element in range; range version.
 * @param rg Range
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Range, class Unary_Function = detail::Identity >
auto
min_value_of(Forward_Range&& rg, Unary_Function&& fn = Unary_Function())
    -> typename std::result_of< Unary_Function(typename std::iterator_traits< decltype(rg.begin()) >::value_type) >::type
{
    return min_value_of(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}


/**
 * Return iterator to maximum element in range; iterator version.
 * @param first Range begin
 * @param last Range end
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Iterator, class Unary_Function = detail::Identity >
Forward_Iterator
max_of(Forward_Iterator first, Forward_Iterator last, Unary_Function&& fn = Unary_Function())
{
    if (first == last) return last;
    auto max_it = first++;
    for (auto it = first; it != last; ++it)
    {
        if (fn(*max_it) < fn(*it))
        {
            max_it = it;
        }
    }
    return max_it;
}

/**
 * Return iterator to maximum element in range; range version.
 * @param rg Range
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Range, class Unary_Function = detail::Identity >
auto
max_of(Forward_Range&& rg, Unary_Function&& fn = Unary_Function())
    -> decltype(rg.end())
{
    static_assert(not std::is_rvalue_reference< Forward_Range >::value,
                  "rvalue reference to range not allowed: returned iterator might not exist");
    return max_of(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Return maximum element in range; iterator version.
 * @param first Range begin
 * @param last Range end
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Iterator, class Unary_Function = detail::Identity >
typename std::result_of< Unary_Function(typename std::iterator_traits< Forward_Iterator >::value_type) >::type
max_value_of(Forward_Iterator first, Forward_Iterator last, Unary_Function&& fn = Unary_Function())
{
    auto it = max_of(first, last, std::forward< Unary_Function >(fn));
    return it != last
        ? fn(*it)
        : typename std::result_of< Unary_Function(typename std::iterator_traits< Forward_Iterator >::value_type) >::type();
}

/**
 * Return maximum element in range; range version.
 * @param rg Range
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Range, class Unary_Function = detail::Identity >
auto
max_value_of(Forward_Range&& rg, Unary_Function&& fn = Unary_Function())
    -> typename std::result_of< Unary_Function(typename std::iterator_traits< decltype(rg.begin()) >::value_type) >::type
{
    return max_value_of(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Return iterator to minimum and maximum elements in range; iterator version.
 * @param first Range begin
 * @param last Range end
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Iterator, class Unary_Function = detail::Identity >
std::pair< Forward_Iterator, Forward_Iterator >
minmax_of(Forward_Iterator first, Forward_Iterator last, Unary_Function&& fn = Unary_Function())
{
    if (first == last) return std::make_pair(last, last);
    auto min_it = first++;
    auto max_it = min_it;
    for (auto it = first; it != last; ++it)
    {
        if (fn(*it) < fn(*min_it))
        {
            min_it = it;
        }
        if (fn(*max_it) < fn(*it))
        {
            max_it = it;
        }
    }
    return std::make_pair(min_it, max_it);
}

/**
 * Return iterator to minimum and maximum elements in range; range version.
 * @param rg Range
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Range, class Unary_Function = detail::Identity >
auto
minmax_of(Forward_Range&& rg, Unary_Function&& fn = Unary_Function())
    -> std::pair< decltype(rg.end()), decltype(rg.end()) >
{
    static_assert(not std::is_rvalue_reference< Forward_Range >::value,
                  "rvalue reference to range not allowed: returned iterator might not exist");
    return minmax_of(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Return minimum and maximum elements in range; iterator version.
 * @param first Range begin
 * @param last Range end
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Iterator, class Unary_Function = detail::Identity >
std::pair< typename std::result_of< Unary_Function(typename std::iterator_traits< Forward_Iterator >::value_type) >::type,
           typename std::result_of< Unary_Function(typename std::iterator_traits< Forward_Iterator >::value_type) >::type >
minmax_value_of(Forward_Iterator first, Forward_Iterator last, Unary_Function&& fn = Unary_Function())
{
    auto p = minmax_of(first, last, std::forward< Unary_Function >(fn));
    return p.first != last
        ? std::make_pair(fn(*p.first), fn(*p.second))
        : std::make_pair(typename std::result_of< Unary_Function(typename std::iterator_traits< Forward_Iterator >::value_type) >::type(),
                         typename std::result_of< Unary_Function(typename std::iterator_traits< Forward_Iterator >::value_type) >::type());
}

/**
 * Return minimum and maximum elements in range; range version.
 * @param rg Range
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Range, class Unary_Function = detail::Identity >
auto
minmax_value_of(Forward_Range&& rg, Unary_Function&& fn = Unary_Function())
    -> std::pair< typename std::result_of< Unary_Function(typename std::iterator_traits< decltype(rg.begin()) >::value_type) >::type,
                  typename std::result_of< Unary_Function(typename std::iterator_traits< decltype(rg.begin()) >::value_type) >::type >
{
    return minmax_value_of(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Compute mean and stdv of a sequence of samples; iterator version.
 * @tparam Float_Type Floating point type to use for calculations.
 * @param first Range begin
 * @param last Range end
 * @param fn Functor used to obtain key from value.
 */
template < typename Float_Type = double, class Input_Iterator, class Unary_Function = detail::Identity >
std::pair< Float_Type, Float_Type >
mean_stdv_of(Input_Iterator it_begin, Input_Iterator it_end, Unary_Function&& fn = Unary_Function())
{
    Float_Type s = 0.0;
    Float_Type s2 = 0.0;
    long unsigned n = 0;
    for (Input_Iterator it = it_begin; it != it_end; ++it)
    {
        s += fn(*it);
        s2 += fn(*it) * fn(*it);
        ++n;
    }
    Float_Type mean = n > 0? s / n : (Float_Type)0;
    Float_Type stdv = n > 1? std::sqrt((s2 - s * mean * 2.0 + mean * mean * (Float_Type)n)/(n - 1)) : (Float_Type)0;
    return std::make_pair(mean, stdv);
}

/**
 * Compute mean and stdv of a sequence of samples; range version.
 * @tparam Float_Type Floating point type to use for calculations.
 * @param rg Range
 * @param fn Functor used to obtain key from value.
 */
template < typename Float_Type = double, class Input_Range, class Unary_Function = detail::Identity >
std::pair< Float_Type, Float_Type >
mean_stdv_of(const Input_Range& rg, Unary_Function&& fn = Unary_Function())
{
    return mean_stdv_of< Float_Type >(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Check if all values in range are equal; iterator version.
 * @param first Range begin
 * @param last Range end
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Iterator, class Unary_Function = detail::Identity >
bool
equal_of(Forward_Iterator first, Forward_Iterator last, Unary_Function&& fn = Unary_Function())
{
    if (first == last) return true;
    for (auto it = std::next(first); it != last; ++it)
    {
        if (fn(*it) != fn(*first))
        {
            return false;
        }
    }
    return true;
}

/**
 * Check if all values in range are equal; range version.
 * @param rg Range
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Range, class Unary_Function = detail::Identity >
bool
equal_of(Forward_Range&& rg, Unary_Function&& fn = Unary_Function())
{
    return equal_of(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Check functor returns true for all values in range; iterator version.
 * @param first Range begin
 * @param last Range end
 * @param fn Functor used to obtain key from value.
 */
using std::all_of;

/**
 * Check functor returns true for all values in range; range version.
 * @param rg Range
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Range, class Unary_Function = detail::Identity >
bool
all_of(Forward_Range&& rg, Unary_Function&& fn = Unary_Function())
{
    return all_of(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Check functor returns true for at least one value in range; iterator version.
 * @param first Range begin
 * @param last Range end
 * @param fn Functor used to obtain key from value.
 */
using std::any_of;

/**
 * Check functor returns true for at least one value in range; range version.
 * @param rg Range
 * @param fn Functor used to obtain key from value.
 */
template < class Forward_Range, class Unary_Function = detail::Identity >
bool
any_of(Forward_Range&& rg, Unary_Function&& fn = Unary_Function())
{
    return any_of(rg.begin(), rg.end(), std::forward< Unary_Function >(fn));
}

/**
 * Accumulate; iterator version.
 */
using std::accumulate;

/**
 * Accumulate; range version.
 */
template < class Forward_Range, typename T >
T accumulate(Forward_Range&& rg, T init)
{
    return accumulate(rg.begin(), rg.end(), init);
}

template < class Forward_Range, typename T, class Binary_Operation >
T accumulate(Forward_Range&& rg, T init, Binary_Operation&& op)
{
    return accumulate(rg.begin(), rg.end(), init, std::forward< Binary_Operation >(op));
}

/**
 * Create os_join_helper given a pair of begin/end iterators.
 */
template < typename Input_Iterator, typename Unary_Function = detail::Identity >
detail::os_join_helper< Input_Iterator, Unary_Function >
os_join(Input_Iterator it_begin, Input_Iterator it_end, const std::string& sep,
        Unary_Function&& fn = Unary_Function())
{
    return detail::os_join_helper< Input_Iterator, Unary_Function >(
        it_begin, it_end, sep, std::forward< Unary_Function >(fn));
} // os_join()

/**
 * Create os_join_helper given a range.
 */
template < typename Input_Range, typename Unary_Function = detail::Identity >
auto
os_join(Input_Range&& rg, const std::string& sep,
        Unary_Function&& fn = Unary_Function())
    -> detail::os_join_helper< decltype(rg.end()), Unary_Function >
{
    return detail::os_join_helper< decltype(rg.end()), Unary_Function >(
        rg.begin(), rg.end(), sep, std::forward< Unary_Function >(fn));
} // os_join()

} // namespace alg

#endif

#ifdef SAMPLE_ALG

/*

Compile with:

g++ -std=c++11 -D SAMPLE_ALG -x c++ alg.hpp -o sample-alg

*/

#include <vector>
#include <array>

using namespace std;
using namespace alg;

int main()
{
    array< array< unsigned, 3 >, 3 > v2 = {{ {{ 42, 1, 15 }}, {{1, 2, 3}}, {{4, 5, 6}} }};
    array< unsigned, 3 >& v = v2[0];
    cout << "v: " << os_join(v.begin(), v.end(), ", ") << endl;
    cout << "v[0]: " << os_join(v.begin(), v.begin() + 1, ", ") << endl;
    cout << "v+1: " << os_join(v.begin(), v.end(), ", ", [] (unsigned i) { return i + 1; }) << endl;
    cout << "v_rg: " << os_join(v, ",") << endl;
    cout << "v+2_rg: " << os_join(v, ",", [] (unsigned i) { return i + 2; }) << endl;
    cout << "u: " << os_join(vector< int >{ 5, 17, 6 }, ";") << endl;
    string s = os_join(vector< char >{ 'a', 'b', 'c' }, "-");
    cout << "s: " << s << endl;
    cout << "min: " << min_value_of(v) << endl;
    cout << "max: " << max_value_of(v) << endl;
    cout << "minmax: " << minmax_value_of(v).first << "," << minmax_value_of(v).second << endl;
}

#endif
