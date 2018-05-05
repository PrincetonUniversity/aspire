// QuickSort.h
//
// Simple template-based quicksort implementation which uses randomized pivots
// and switches to an insertion sort for small partitions.  These routines
// order an auxillary array along with the array being sorted.  This could
// have been accomplished using the STL library, but it would have entailed
// some overhead.
//
// Obviously, this is used to sort lists of (index, value) pairs for
// sparse vectors and matrices.
//
// Dependencies:        none
//
// (c) 2005 James Bremer (james.bremer@yale.edu)

#ifndef QUICKSORT__H
#define QUICKSORT__H

#include <stdlib.h>  // rand
#include <math.h>

// switch to insertion sort when partitions become smaller than this
#define QSORT_MIN_PART 16

inline int our_round(double x)
{
   int flr = (int)floor(x);

   if (x - (double)flr > .5)
      return flr+1;
   else
      return flr;
}

// get a pseudorandom number in the range left, ..., right
inline int Rand(int left, int right)
{
   return left+ our_round((double)rand()/(double)RAND_MAX * ((double)right - (double)left));
}


// template-based insert sort
template<class T, class I>
inline void InsertSort(T* data, I *aux_data, int left, int right )
{
   for(int j=left+1; j <= right; j++) {

      // find the correct position for data[j] in the sorted list
      // data[left ... j-1]
      for(int i=left; i < j; i++) {

         // check to see if i is the right position
         if(data[j] < data[i]) {
            T temp_t = data[j];
            I temp_i = aux_data[j];
            for(int k=j; k > i; k--) {
               data[k]     = data[k-1];
               aux_data[k] = aux_data[k-1];
            }
            data[i]     = temp_t;
            aux_data[i] = temp_i;
            break;
         }
      }
   }
}
template<class T, class I>
void QuickSort(T *values, I *aux_data, int left, int right)
{
   if(right-left+1 > QSORT_MIN_PART) {

      int random = Rand(left, right);

      T pivot = values[random];  // random pivot
      int i=left-1;
      int j=right+1;

      while(1) {

         do { i++; } while(values[i] < pivot);
         do { j--; } while(values[j] > pivot);

         if(i<j) {

            // swap now
            T  swap_d   = values[i];
            I  swap_i   = aux_data[i];
            values[i]   = values[j];
            aux_data[i] = aux_data[j];
            values[j]   = swap_d;
            aux_data[j] = swap_i;
         } else {
            QuickSort(values, aux_data, left, j);
            QuickSort(values, aux_data, j+1, right);
            break;
         }
      }

   } else
      InsertSort(values, aux_data, left, right);
}


// template-based insert sort
template<class T>
inline void InsertSort(T* data, int left, int right )
{
   for(int j=left+1; j <= right; j++) {

      // find the correct position for data[j] in the sorted list
      // data[left ... j-1]
      for(int i=left; i < j; i++) {

         // check to see if i is the right position
         if(data[j] < data[i]) {
            T temp_t = data[j];
            for(int k=j; k > i; k--) {
               data[k]     = data[k-1];
            }
            data[i]     = temp_t;
            break;
         }
      }
   }
}
template<class T>
void QuickSort(T *values, int left, int right)
{
   if(right-left+1 > QSORT_MIN_PART) {

      int random = Rand(left, right);

      T pivot = values[random];  // random pivot
      int i=left-1;
      int j=right+1;

      while(1) {

         do { i++; } while(values[i] < pivot);
         do { j--; } while(values[j] > pivot);

         if(i<j) {

            // swap now
            T  swap_d   = values[i];
            values[i]   = values[j];
            values[j]   = swap_d;
         } else {
            QuickSort(values, left, j);
            QuickSort(values, j+1, right);
            break;
         }
      }

   } else
      InsertSort(values, left, right);
}

#endif
