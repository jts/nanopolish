/// @author    Matei David, Ontario Institute for Cancer Research
/// @version   2.0
/// @date      2013-2015
/// @copyright MIT Public License
///
/// A parallel for loop.
/// A parallel for loop that can re-sort its output, implemented
/// in C++11.

#ifndef __PFOR_HPP
#define __PFOR_HPP

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <functional>
#include <thread>
#include <mutex>
#include <chrono>

/// C++11 implementation of a parallel for loop that can re-sort its output.
///
/// Given a number of threads and a chunk size, the total work is split into
/// chunks, each chunk is processed independently by a thread, and the outputs
/// are reordered and processed in the same order they were read. More
/// specifically:
///   - The master thread starts `num_threads` worker threads, then waits for
///     all workers. Each worker proceeds as follows:
///   - Allocate a vector of `chunk_size` objects of type `Input` which are
///     default-initialized.
///   - Execute `do_before()`.
///   - Repeatedly, until done:
///     - Inside an input critical section: call `get_item()` repeatedly, getting
///       up to `chunk_size` items. The `get_item()` function is expected to save
///       an item in the Input given to it.
///     - Allocate a new default-initialized Chunk_Output object (on the heap).
///     - Repeatedly call `process_item()` in the order the items appear in the
///       thread local buffer. The calls are given the same Chunk_Output object.
///     - Inside an output critical section: the Chunk_Output object is added to
///       a heap. As long as the chunk that must be output is in found in the
///       heap, the corresponding Chunk_Output object is destroyed.
///   - Execute `do_after()`.
///
/// Note: To re-sort the output, Chunk_Output should contain an
/// ostringstream object that is used for output in `process_item()`, and that
/// is flushed to std::cout upon destruction (which happens inside the output
/// critical section).
///
/// @tparam Input Type of the input items.
/// @tparam Chunk_Output Type of an object created during the processing of each
/// chunk. Each such object is tagged with its chunk number, and placed in a
/// priority queue. The objects are destroyed inside an output critical section
/// in the order of their chunk numbers.
/// @param do_before Function executed by each thread before any actual work.
/// @param get_item Function executed inside an input critical section. If an item is
/// available, the function is expected to save it in the `Input` location given
/// as parameter and return true. If the items have been exhausted, the function
/// is execpted to return false.
/// @param process_item Function executed to process the given `Input` item. The
/// `Output_Chunk` object will be the same for all items in that chunk.
/// @param do_after Function executed by each thread after all work.
/// @param num_threads Number of worker threads to spawn.
/// @param chunk_size Number of items each thread should process in one chunk.
template < typename Input, typename Chunk_Output >
void pfor(unsigned num_threads,
          size_t chunk_size,
          std::function< bool(Input&) > get_item,
          std::function< void(Input&, Chunk_Output&) > process_item,
          std::function< void(size_t, size_t) > progress_report = nullptr,
          size_t progress_count = 0);

template < typename Input, typename Chunk_Output >
void pfor(unsigned num_threads,
          size_t chunk_size,
          std::function< bool(Input&) > get_item,
          std::function< void(Input&, Chunk_Output&) > process_item,
          std::function< void(void) > do_before,
          std::function< void(void) > do_after,
          std::function< void(size_t, size_t) > progress_report = nullptr,
          size_t progress_count = 0);

/// @cond
namespace detail
{

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

template < typename Input, typename Chunk_Output >
void do_work(unsigned tid, std::reference_wrapper< Common_Storage< Input, Chunk_Output > > cs_wrap)
{
    Common_Storage< Input, Chunk_Output >& cs = cs_wrap;
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
        Chunk_Output_Wrapper< Chunk_Output >* cow_p = new Chunk_Output_Wrapper< Chunk_Output >(tid, cid);
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
                delete cow_p;
                cs.cid_out++;
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
          std::function< void(void) > do_before,
          std::function< void(void) > do_after,
          std::function< void(size_t, size_t) > progress_report,
          size_t progress_count)
{
    std::vector< std::thread > thread_v;
    detail::Common_Storage< Input, Chunk_Output > cs;
    cs.get_item = get_item;
    cs.process_item = process_item;
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
    for (unsigned i = 0; i < num_threads; ++i)
    {
        thread_v.emplace_back(detail::do_work< Input, Chunk_Output >, i, std::ref(cs));
    }
    for (auto& t : thread_v)
    {
        t.join();
    }
    if (cs.progress_report)
    {
        auto crt_time = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast< std::chrono::seconds >(crt_time - cs.start_time);
        cs.progress_report(cs.item_count, elapsed.count());
    }
} // pfor()

template < typename Input, typename Chunk_Output >
void pfor(unsigned num_threads,
          size_t chunk_size,
          std::function< bool(Input&) > get_item,
          std::function< void(Input&, Chunk_Output&) > process_item,
          std::function< void(size_t, size_t) > progress_report,
          size_t progress_count)
{
    pfor< Input, Chunk_Output >(
        num_threads, chunk_size, get_item, process_item, nullptr, nullptr, progress_report, progress_count);
} // pfor()

#endif
