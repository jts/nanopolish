/// @author    Matei David, Ontario Institute for Cancer Research
/// @version   2.0
/// @date      2013-2015
/// @copyright MIT Public License
///
/// A parallel for loop.
/// A parallel for loop that can re-sort its output, implemented in C++11.

#ifndef __PFOR_HPP
#define __PFOR_HPP

#include <cassert>
#include <chrono>
#include <functional>
#include <mutex>
#include <queue>
#include <vector>
#include <thread>
#include <type_traits>

namespace pfor
{

/// C++11 implementation of a parallel for loop.
///
/// Given a number of threads and a chunk size, the total work is split into
/// chunks, each chunk is processed independently by a thread. Optionally, the
/// outputs can be reordered in the same order they were read. More
/// specifically:
///   - The master thread starts `num_threads` worker threads, then waits for
///     all workers. Each worker proceeds as follows:
///   - Allocate a vector of `chunk_size` objects of type `Input` which are
///     default-initialized.
///   - Execute `do_before()`.
///   - Repeatedly, until done:
///     - Inside an input critical section: call `get_item()` repeatedly,
///       getting up to `chunk_size` items. The `get_item()` function should
///       save an item in the `Input` location given as parameter.
///     - If `Chunk_Output` is given:
///       - Allocate a new default-initialized Chunk_Output object on the heap.
///     - Repeatedly call `process_item()` in the order the items appear in the
///       thread local buffer. If `Chunk_Output` is given, the calls are passed
///       the same `Chunk_Output` object.
///     - If `Chunk_Output` is given: Inside an output critical section, the
///       `Chunk_Output` object is added to a heap. As long as the chunk that
///       must be output next is in found in the heap, `output_chunk` is called
///       on that object, after which it is destroyed.
///   - Execute `do_after()`.
///
/// To ENABLE output sorting:
/// - `Chunk_Output` should contain one (or more) ostringstream object(s);
/// - `process_item()` should take 2 arguments, the second being a
/// `Chunk_Output` object; the ostringstream inside that `Chunk_Output` should
/// be used for output inside `process_item()`;
/// - `output_chunk()` should flush the ostringstream to output.
///
/// To DISABLE output sorting:
/// - `Chunk_Output` should not be specified;
/// - `process_item()` should take 1 argument only;
/// - `output_chunk()` should not be specified.
///
/// @tparam Input Type of the input items.
/// @tparam Chunk_Output (optional) Type of an object created during the
/// processing of each chunk. Each such object is tagged with its chunk number,
/// and placed in a priority queue. The objects are destroyed inside an output
/// critical section in the order of their chunk numbers.
/// @param num_threads Number of worker threads to spawn.
/// @param chunk_size Number of items each thread should process in one chunk.
/// @param get_item Function executed inside an input critical section. If an
/// item is available, the function should save it in the `Input` location given
/// as parameter, and return true. If the items have been exhausted, the
/// function should return false.
/// @param process_item Function executed to process the given `Input` item. If
/// output sorting is enabled, the function takes as second argument a
/// `Chunk_Output` object that will be the same for all items in that chunk.
/// @param output_chunk (optional) Function called on a `Chunk_Output`
/// object inside an output critical section before it is destroyed.
/// @param do_before Function executed by each thread before any actual work.
/// @param do_after Function executed by each thread after all work.
/// @progress_report Function that will be called to report progress with 2
/// arguments: the number of items processed, and the number of seconds
/// elapsed so far.
/// @progress_count Number of items that should be processed between calls to
/// `progress_report`. (This will be rounded down to be divisible by
/// `chunk_size`).

template < typename Input >
void pfor(unsigned num_threads,
          size_t chunk_size,
          std::function< bool(Input&) > get_item,
          std::function< void(Input&) > process_item,
          std::function< void(size_t, size_t) > progress_report = nullptr,
          size_t progress_count = 0);

template < typename Input, typename Chunk_Output >
void pfor(unsigned num_threads,
          size_t chunk_size,
          std::function< bool(Input&) > get_item,
          std::function< void(Input&, Chunk_Output&) > process_item,
          std::function< void(Chunk_Output&) > output_chunk,
          std::function< void(size_t, size_t) > progress_report = nullptr,
          size_t progress_count = 0);

template < typename Input, typename Chunk_Output >
void pfor(unsigned num_threads,
          size_t chunk_size,
          std::function< bool(Input&) > get_item,
          std::function< void(Input&, Chunk_Output&) > process_item,
          std::function< void(Chunk_Output&) > output_chunk,
          std::function< void(void) > do_before,
          std::function< void(void) > do_after,
          std::function< void(size_t, size_t) > progress_report = nullptr,
          size_t progress_count = 0);

/// @cond
namespace detail
{

struct empty {};

template < typename Chunk_Output >
struct Chunk_Output_Wrapper
    : public Chunk_Output
{
    Chunk_Output_Wrapper(unsigned tid_, unsigned cid_)
        : Chunk_Output(), tid(tid_), cid(cid_) {}

    unsigned tid;
    unsigned cid;
}; // struct Chunk_Output_Wrapper

template < typename Chunk_Output >
struct Chunk_Output_Wrapper_Ptr_Comparator
{
    bool operator () (const Chunk_Output_Wrapper< Chunk_Output >* lhs_p,
                      const Chunk_Output_Wrapper< Chunk_Output >* rhs_p)
    {
        return lhs_p->cid > rhs_p->cid;
    }
};

template < typename Chunk_Output >
using Output_Heap = std::priority_queue< Chunk_Output_Wrapper< Chunk_Output >*,
                                    std::vector< Chunk_Output_Wrapper< Chunk_Output >* >,
                                    Chunk_Output_Wrapper_Ptr_Comparator< Chunk_Output > >;

template < typename Input, typename Chunk_Output >
struct Common_Storage
{
    std::function< bool(Input&) > get_item;
    std::function< void(Input&, Chunk_Output&) > process_item;
    std::function< void(Chunk_Output&) > output_chunk;
    std::function< void(void) > do_before;
    std::function< void(void) > do_after;
    std::function< void(size_t, size_t) > progress_report;
    size_t chunk_size;
    size_t progress_count;
    size_t item_count;
    size_t cid_in;
    size_t cid_out;
    std::chrono::system_clock::time_point start_time;
    std::mutex input_mutex;
    std::mutex output_mutex;
    Output_Heap< Chunk_Output > h;
}; // struct Common_Storage

template < typename Input, typename Chunk_Output, bool sort_output >
void do_work(unsigned tid, std::reference_wrapper< Common_Storage< Input, Chunk_Output > > cs_wrap)
{
    Common_Storage< Input, Chunk_Output >& cs = cs_wrap;
    Chunk_Output_Wrapper< Chunk_Output > _empty_wrapper(tid, 0);
    std::vector< Input > buff(cs.chunk_size);
    size_t load = 0;
    size_t cid;
    bool done = false;
    if (cs.do_before) cs.do_before();
    while (not done)
    {
        load = 0;
        // input critical section
        {
            std::lock_guard< std::mutex > input_lock(cs.input_mutex);
            cid = cs.cid_in++;
            while (load < cs.chunk_size and cs.get_item(buff[load]))
            {
                ++load;
            }
            done = (load < cs.chunk_size);
        }
        // parallel work
        if (load == 0)
        {
            break;
        }
        Chunk_Output_Wrapper< Chunk_Output >* cow_p = &_empty_wrapper;
        if (sort_output)
        {
            cow_p = new Chunk_Output_Wrapper< Chunk_Output >(tid, cid);
        }
        for (size_t i = 0; i < load; ++i)
        {
            cs.process_item(buff[i], *cow_p);
        }
        // output critical section
        {
            std::lock_guard< std::mutex > output_lock(cs.output_mutex);
            cs.item_count += load;
            if (cs.progress_report and cs.item_count % cs.progress_count == 0)
            {
                auto crt_time = std::chrono::system_clock::now();
                auto elapsed = std::chrono::duration_cast< std::chrono::seconds >(crt_time - cs.start_time);
                cs.progress_report(cs.item_count, elapsed.count());
            }
            if (sort_output)
            {
                cs.h.push(cow_p);
                while (cs.h.size() > 0)
                {
                    cow_p = cs.h.top();
                    assert(cow_p->cid >= cs.cid_out);
                    if (cow_p->cid > cs.cid_out)
                    {
                        break;
                    }
                    cs.h.pop();
                    if (cs.output_chunk)
                    {
                        cs.output_chunk(*cow_p);
                    }
                    delete cow_p;
                    cs.cid_out++;
                }
            }
        }
    }
    if (cs.do_after) cs.do_after();
} // do_work()

} // namespace detail
/// @endcond

template < typename Input, typename Chunk_Output >
void pfor(unsigned num_threads,
          size_t chunk_size,
          std::function< bool(Input&) > get_item,
          std::function< void(Input&, Chunk_Output&) > process_item,
          std::function< void(Chunk_Output&) > output_chunk,
          std::function< void(void) > do_before,
          std::function< void(void) > do_after,
          std::function< void(size_t, size_t) > progress_report,
          size_t progress_count)
{
    std::vector< std::thread > thread_v;
    detail::Common_Storage< Input, Chunk_Output > cs;
    cs.get_item = get_item;
    cs.process_item = process_item;
    cs.output_chunk = output_chunk;
    cs.do_before = do_before;
    cs.do_after = do_after;
    cs.progress_report = progress_report;
    cs.chunk_size = chunk_size;
    if (progress_count == 0) progress_count = 10 * num_threads * chunk_size;
    cs.progress_count = std::max(progress_count / chunk_size, (size_t)1) * chunk_size;
    cs.item_count = 0;
    cs.cid_in = 0;
    cs.cid_out = 0;
    cs.start_time = std::chrono::system_clock::now();
    static const bool sort_output = not std::is_same< Chunk_Output, detail::empty >::value;
    for (unsigned i = 0; i < num_threads; ++i)
    {
        thread_v.emplace_back(detail::do_work< Input, Chunk_Output, sort_output >, i, std::ref(cs));
    }
    for (auto& t : thread_v)
    {
        t.join();
    }
    if (cs.progress_report and cs.item_count % cs.progress_count != 0)
    {
        auto crt_time = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast< std::chrono::seconds >(crt_time - cs.start_time);
        cs.progress_report(cs.item_count, elapsed.count());
    }
} // pfor()

template < typename Input >
void pfor(unsigned num_threads,
          size_t chunk_size,
          std::function< bool(Input&) > get_item,
          std::function< void(Input&) > process_item,
          std::function< void(size_t, size_t) > progress_report,
          size_t progress_count)
{
    pfor< Input, detail::empty >(
        num_threads,
        chunk_size,
        get_item,
        [&] (Input& in, detail::empty&) { process_item(in); },
        nullptr,
        nullptr,
        nullptr,
        progress_report,
        progress_count);
} // pfor()

template < typename Input, typename Chunk_Output >
void pfor(unsigned num_threads,
          size_t chunk_size,
          std::function< bool(Input&) > get_item,
          std::function< void(Input&, Chunk_Output&) > process_item,
          std::function< void(Chunk_Output&) > output_chunk,
          std::function< void(size_t, size_t) > progress_report,
          size_t progress_count)
{
    pfor< Input, Chunk_Output >(
        num_threads,
        chunk_size,
        get_item,
        process_item,
        output_chunk,
        nullptr,
        nullptr,
        progress_report,
        progress_count);
} // pfor()

} // namespace pfor

#endif

#ifdef PFOR_SAMPLE
/*

Compile:

g++ -std=c++11 -D PFOR_SAMPLE -x c++ pfor.hpp -o pfor

*/

#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

int main()
{
    //
    // pfor loop WITHOUT output sorting
    //
    vector< unsigned > v(1000);
    unsigned crt_idx = 0;
    pfor::pfor< unsigned >(
        // num_threads
        4,
        // chunk_size
        10,
        // get_item
        [&] (unsigned& i) {
            if (crt_idx >= v.size()) return false;
            i = crt_idx++;
            return true;
        },
        // process_item
        [&] (unsigned& i) {
            v[i] = i*i;
        },
        //
        // The following 2 parameters are optional
        //
        // progress_report
        [&] (size_t items, size_t seconds) {
            cout << "processed " << items << " items in " << seconds << " seconds" << endl;
        },
        // progress_count
        110);
    //
    // pfor loop WITH output sorting
    //
    crt_idx = 0;
    pfor::pfor< unsigned, std::ostringstream >(
        // num_threads
        4,
        // chunk_size
        1,
        // get_item
        [&] (unsigned& i) {
            if (crt_idx >= 10) return false;
            i = crt_idx++;
            return true;
        },
        // process_item
        [&] (unsigned& i, std::ostringstream& os) {
            thread::id this_id = this_thread::get_id();
            os << setw(3) << right << this_id << ": " << i << "*" << i << "=" << i*i << endl;
        },
        // output_chunk
        [&] (std::ostringstream& os) {
            cout << os.str();
        });
}

#endif
