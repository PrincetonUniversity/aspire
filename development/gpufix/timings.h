/* Timing routines
*
* Yoel Shkolnisky, July 2016
*/

#ifndef __CU_TIMINGS__

#define __CU_TIMINGS__

#include <time.h>

clock_t start;
clock_t time_diff;
int msec;

#define TIC (start=clock())
#define TOC time_diff = clock() - start; msec = time_diff * 1000 / CLOCKS_PER_SEC
#define TOCM(X) TOC; printf("%s %d milliseconds\n", X,msec);

#endif