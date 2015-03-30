///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Profiler.h -- Lightweight macro-based function profiler.
//
#ifndef PROFILER_H
#define PROFILER_H

#include <time.h>
#include <iostream>

//#define USE_PROFILER 1

#if USE_PROFILER

// Change this to determine how often the profile should print
#define PROFILE_TICKS_BEFORE_PRINT 1000

// This class writes the lifespan of the object
// to the output variable, in nanoseconds
class TimeTracker
{
    public:
        TimeTracker(size_t & output) : m_output(output)
        {
             timespec start;
             clock_gettime(CLOCK_REALTIME, &start);
             m_start_ns = start.tv_sec * 1000000000 + start.tv_nsec;
        }

        ~TimeTracker()
        {
             timespec end;
             clock_gettime(CLOCK_REALTIME, &end);
             size_t end_ns = end.tv_sec * 1000000000 + end.tv_nsec;

             // Update the result using an atomic compare and swap
             size_t diff = end_ns - m_start_ns;
             while(!__sync_bool_compare_and_swap(&m_output, m_output, m_output + diff)) {}
        }

    private:
        size_t m_start_ns;
        size_t& m_output;
};

// Place this macros at the start of the function you wish the profile
// The static variable updates are done via atomic compare and swaps so
// the profiling should be threadsafe
#define PROFILE_FUNC(x) static std::string __profile_name = x; \
                        static size_t __profile_iterations = 0; \
                        static size_t __profile_total_nanoseconds = 0; \
                        double micro_seconds = (double)__profile_total_nanoseconds / 1000.0f; \
                        double avg_per_iteration = micro_seconds / __profile_iterations; \
                        while(!__sync_bool_compare_and_swap(&__profile_iterations, __profile_iterations, __profile_iterations + 1)) { } \
                        if(__profile_iterations % PROFILE_TICKS_BEFORE_PRINT == 0) \
                            fprintf(stderr, "[Profile] count: %zu time: %.0lf ms avg: %.0lf us func: %s\n", __profile_iterations, micro_seconds / 1000, avg_per_iteration, __profile_name.c_str()); \
                        TimeTracker __profile_timer(__profile_total_nanoseconds);
                         
#else

// Eliminate the macro
#define PROFILE_FUNC(x)

#endif // #ifdef HAVE_CLOCK_GETTIME

#endif // #ifndef PROFILER_H
