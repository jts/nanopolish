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
        
        Progress(const std::string message) : m_message(message), m_os(std::cerr)
        {
            timespec start;
            clock_gettime(CLOCK_REALTIME, &start);
            m_start_time = start.tv_sec;
        }

        // derived from: http://stackoverflow.com/a/14539953/378881
        void print(float progress) const
        {
            // get current time
            timespec now;
            clock_gettime(CLOCK_REALTIME, &now);

            // print 
            int max_leader = 40;
            int bar_width = 50;

            std::string leader;
            if(m_message.size() > max_leader) {
                leader = m_message.substr(0, max_leader - 3) + "..."; // truncate
            } else {
                leader = m_message + std::string(max_leader - m_message.size(), ' '); // pad
            }
            
            m_os << leader << " [";
            
            int pos = bar_width * progress;
            for (int i = 0; i < bar_width; ++i) {
                if (i < pos) m_os << "=";
                else if (i == pos) m_os << ">";
                else m_os << " ";
            }
            m_os << "] " << int(progress * 100.0) << "% in " << now.tv_sec - m_start_time << "s\r";
            m_os.flush();
        }
        
        // 
        void end() const
        {
            print(1.0);
            std::cerr << std::endl;
        }

    private:

        std::ostream& m_os;
        std::string m_message;
        size_t m_start_time;
};

#endif
