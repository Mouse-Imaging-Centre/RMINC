#ifndef __MINC2_IO_H__
#define __MINC2_IO_H__

#include <minc2.h>
#include <stdio.h>

// maximum number of frames permitted in a volume
#define MAX_FRAMES 100
// maximum buffer size for attributes
#define MAX_STRING_LEN   256
// R return list size (number of list elements)
#define R_RTN_LIST_LEN 15
// debug switch
#define R_DEBUG_mincIO 0

//#define FALSE 0
//#define TRUE 1
//typedef int boolean


//typedef enum { FALSE, TRUE } boolean

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>



#endif
