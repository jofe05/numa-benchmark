// Based on numatest.cpp by James Brock
// http://stackoverflow.com/questions/7259363/measuring-numa-non-uniform-memory-access-no-observable-asymmetry-why
//
// Changes by Andreas Kloeckner, 10/2012:
// - Rewritten in C + OpenMP
// - Added contention tests

#define _GNU_SOURCE


#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

double mysecond()
{
        struct timeval tp;
        int i;

        i = gettimeofday(&tp,NULL);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

double col_wise(int **m, int dim)
{
  double t1, t2;
  int i, j;
  t1 = mysecond();
  for (i = 0;i< dim; i++){
        for (j=0;j <dim;j++){
            m[i][j]= i*j;
        }
    }
    t2 = mysecond();
    return t2-t1;
}

double row_wise(int **m, int dim)
{
  double t1, t2;
  int i, j;
  t1 = mysecond();
  for (i = 0;i< dim; i++){
        for (j=0;j <dim;j++){
            m[j][i]= i*j;
        }
    }
    t2 = mysecond();
    return t2-t1;
}








int main(int argc, const char **argv)
{
    int t;
    size_t i,j;
    size_t dim = 1024 ;
    typedef int **Matrix;
    Matrix m = malloc(sizeof(int *) * dim);

	  for ( i = 0; i < dim; i++ )
    		m[i] = malloc(sizeof(int) * dim);
    t = col_wise(m, dim);
    printf("Execution Colomun wise %d: \n",t);
      

  free(m);


  return 0;

}
