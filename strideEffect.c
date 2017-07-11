// Based on numatest.cpp by James Brock
// http://stackoverflow.com/questions/7259363/measuring-numa-non-uniform-memory-access-no-observable-asymmetry-why
//
// Changes by Andreas Kloeckner, 10/2012:
// - Rewritten in C + OpenMP
// - Added contention tests



#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double mysecond()
{
        struct timeval tp;
        int i;

        i = gettimeofday(&tp,NULL);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


double measure_access(void *x, int array_size, int stride)
{
  double t1, t2;
  
  t1 = mysecond();
  for(int j = 0; j < array_size ; j += stride)
    {
      *(((char*)x) + (j*16 % array_size)) += 1;
    }
  
  t2 = mysecond();

  return t2 - t1;
}


int main(int argc, const char **argv)
{
 // char *x;
  const int cache_line_size = 64;
  const int array_size = 100*1024*1024;
  double t_base;
//  x = (char *) malloc(array_size);
  
  for( int i = 1; i < 2  ; i *= 2 )
  {
    char *x;

    x = (char *) malloc(array_size);
    double t = measure_access(x, array_size, i);
    if ( i == 1 ){
    t_base = t;
    }
    printf("Array Size: %d Stride %d Time: %f SpeedUP: %f \n",array_size, i, t, t_base/t);
    free(x);

  } 
  //free(x);
  return 0;

}
