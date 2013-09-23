// Debug functions
//
// Usage:
// 1. To enable debug, uncomment the sybmol DEBUG.
// 2. Call SETUPLOG
// 3. Call OPENLOG with the filename of the log file as argument. 
//    For example OPENLOG("log.txt"). 
//    Call OPENLOG(NULL) to print messages to STDIN.
// 4. Call CLOSELOG when finished using logging.
// 5. Use the functions DUMPVAR, DUMPMSG, DUMPARR, DUMPIARR to print log messages.
//		DUMPVAR		Print the value of a given variable.
//		DUMPMSG		Print text message.
//      DUMPARR		Print double array.
//		DUMPIARR	Print integer array.
//
// Yoel Shkolnisky, February 2010.

#ifndef __DEBUG_H__
#define __DEBUG_H__

//#undef DEBUG
#define DEBUG
//#define MATLAB

#include <stdio.h>
#include <assert.h>

/*******************************
 * Macros for debugging messages
 ********************************/

extern FILE *logfid;

#ifdef DEBUG
	void printarr(FILE* logfid, const double *A, int m, int n);
	void printiarr(FILE* logfid, const int *A, int m, int n);
    #define SETUPLOG\
        FILE *logfid;\
        char logfname[80];
#else
    #define SETUPLOG
#endif


#ifdef DEBUG	
    #define OPENLOG(fname)\
		if (fname!=NULL){\
			sprintf(logfname,"%s.log",fname);\
			logfid=fopen(logfname,"w");\
			if (logfid==NULL){\
			    printf("Cannot open log file ""%s"".\n",logfname);\
				return;\
			}\
		}\
		else{printf("Using stdin...\n"); logfid=stdin;}
#else
    #define OPENLOG(fname) printf("YYY\n");
#endif


#ifdef DEBUG
    #define CLOSELOG\
        fclose(logfid);
#else
    #define CLOSELOG
#endif


#ifdef DEBUG
    #define DUMPVAR(prefix, x, fmt)\
        fprintf(logfid,"%s (L%u) %s%s=" fmt "\n",__FILE__,__LINE__, prefix,#x, x);\
        fflush (logfid);
    #define DUMPMSG(msg)\
        fprintf(logfid,"%s (L%u) %s",__FILE__,__LINE__,msg);\
        fflush(logfid);
    #define DUMPARR(a,m,n)\
        fprintf(logfid,"%s (L%u) %s (%d X %d)=\n",__FILE__,__LINE__,#a,m,n);\
        printarr(logfid,a,m,n);\
        fflush(logfid);
    #define DUMPIARR(a,m,n)\
        fprintf(logfid,"%s (L%u) %s (%d X %d)=\n",__FILE__,__LINE__,#a,m,n);\
        printiarr(logfid,a,m,n);\
        fflush(logfid);    
#else
    #define DUMPVAR(prefix,x,fmt)
    #define DUMPMSG(msg)
    #define DUMPARR(a,m,n)
    #define DUMPIARR(a,m,n)
#endif

#endif
