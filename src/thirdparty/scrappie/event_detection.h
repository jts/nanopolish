#ifndef EVENT_DETECTION_H
#    define EVENT_DETECTION_H

#    include "scrappie_structures.h"

typedef struct {
    size_t window_length1;
    size_t window_length2;
    float threshold1;
    float threshold2;
    float peak_height;
} detector_param;


static detector_param const event_detection_defaults = {
    .window_length1 = 3,
    .window_length2 = 6,
    .threshold1 = 1.4f,
    .threshold2 = 9.0f,
    .peak_height = 0.2f
};

static detector_param const event_detection_rna = {
    .window_length1 = 7,
    .window_length2 = 14,
    .threshold1 = 2.5f,
    .threshold2 = 9.0f,
    .peak_height = 1.0f
};


event_table detect_events(raw_table const rt, detector_param const edparam);

#endif                          /* EVENT_DETECTION_H */
