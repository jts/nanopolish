#pragma once
#ifndef SCRAPPIE_STDLIB
#    define SCRAPPIE_STDLIB
#    include <assert.h>
#    include <err.h>
#    include <stdlib.h>
#    include <string.h>

#    if defined(CHAOSMONKEY) && ! defined(BANANA)
#        define FUZZTHRESH (long int)(CHAOSMONKEY * RAND_MAX)
#        define malloc(A) (warnx("Opportunity for chaos at %s:%d", __FILE__, __LINE__),  rand() > FUZZTHRESH) \
                             ? malloc(A) \
                             : (warnx("Injected NULL at %s:%d", __FILE__, __LINE__), NULL)

#        define calloc(A, B) (warnx("Opportunity for chaos at %s:%d", __FILE__, __LINE__), rand() > FUZZTHRESH) \
                             ? calloc(A, B) \
                             : (warnx("Injected NULL at %s:%d", __FILE__, __LINE__), NULL)

#        define aligned_alloc(A, B) (warnx("Opportunity for chaos at %s:%d", __FILE__, __LINE__), rand() > FUZZTHRESH) \
                             ? aligned_alloc(A, B) \
                             : (warnx("Injected NULL at %s:%d", __FILE__, __LINE__), NULL)
#    else
#        define malloc(A) malloc(A)
#        define calloc(A, B) calloc(A, B)
#        define aligned_alloc(A, B) aligned_alloc(A, B)
#    endif /* CHAOSMONKEY */

#    ifdef ABORT_ON_NULL
#       define RETURN_NULL_IF(A, B) \
		if (A) {  	\
			warnx("Failure at %s : %d", __FILE__, __LINE__);	\
            abort(); \
		}
#    else
#        define RETURN_NULL_IF(A, B) if (A) { return B; }
#    endif


//  strlen on NULL is undefined.  Define it.
#    define strlen(A) (NULL != A) ? strlen(A) : 0
#endif /* SCRAPPIE_STDLIB */
