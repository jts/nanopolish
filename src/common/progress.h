//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// progress -- Small class to write the progress
// of a long running task to stderr
//
#ifndef PROGRESS_H
#define PROGRESS_H

#include <time.h>

class Progress
{
    public:
        
        Progress(const std::string message) : m_os(std::cerr), m_message(message)
        {
#if HAVE_CLOCK_GETTIME            
            timespec start;
            clock_gettime(CLOCK_REALTIME, &start);
            m_start_time = start.tv_sec;
#else
            m_start_time = 0;
#endif
        }

        // derived from: http://stackoverflow.com/a/14539953/378881
        void print(float progress) const
        {

            // print 
            unsigned max_leader = 40;
            unsigned bar_width = 50;

            std::string leader;
            if(m_message.size() > max_leader) {
                leader = m_message.substr(0, max_leader - 3) + "..."; // truncate
            } else {
                leader = m_message + std::string(max_leader - m_message.size(), ' '); // pad
            }
            
            m_os << leader << " [";
            
            unsigned pos = bar_width * progress;
            for (unsigned i = 0; i < bar_width; ++i) {
                if (i < pos) m_os << "=";
                else if (i == pos) m_os << ">";
                else m_os << " ";
            }
            m_os << "] " << unsigned(progress * 100.0) << "% in " << get_elapsed_seconds() << "s\r";
            m_os.flush();
        }
        
        // 
        void end() const
        {
            print(1.0);
            std::cerr << std::endl;
        }

        double get_elapsed_seconds() const
        {
            // get current time
#if HAVE_CLOCK_GETTIME            
            timespec now;
            clock_gettime(CLOCK_REALTIME, &now);
            return now.tv_sec - m_start_time;
#else
            return 0;
#endif
        }

    private:

        std::ostream& m_os;
        std::string m_message;
        size_t m_start_time;
};

#endif
