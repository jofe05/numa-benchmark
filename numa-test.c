// Based on numatest.cpp by James Brock
// http://stackoverflow.com/questions/7259363/measuring-numa-non-uniform-memory-access-no-observable-asymmetry-why
//
// Changes by Andreas Kloeckner, 10/2012:
// - Rewritten in C + OpenMP
// - Added contention tests

#define _GNU_SOURCE
#include <numa.h>
#include <numaif.h>
#include <sched.h>
#include <stdio.h>
#include <pthread.h>
#include <omp.h>
#include <assert.h>
#include "timing.h"

void pin_to_core(size_t core)
{
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(core, &cpuset);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
}

void print_bitmask(const struct bitmask *bm)
{
  for(size_t i=0; i<bm->size; ++i)
    printf("%d", numa_bitmask_isbitset(bm, i));
}


double measure_access(void *x, size_t array_size, size_t ntrips)
{
  timestamp_type t1;
  get_timestamp(&t1);

  for (size_t i = 0; i<ntrips; ++i)
    for(size_t j = 0; j<array_size; ++j)
    {
      *(((char*)x) + ((j * 1009) % array_size)) += 1;
    }

  timestamp_type t2;
  get_timestamp(&t2);

  return timestamp_diff_in_seconds(t1, t2);
}

double measure_access_sequential(void *x, size_t array_size, size_t ntrips)
{ 
  timestamp_type t1;
  get_timestamp(&t1);
  
  for (size_t i = 0; i<ntrips; ++i) 
    for(size_t j = 0; j<array_size; ++j)
    { 
      *(((char*)x) + ((j*64) % array_size))  += 1;
    }
  
  timestamp_type t2;
  get_timestamp(&t2);
  
  return timestamp_diff_in_seconds(t1, t2);
}


int main(int argc, const char **argv)
{
  numa_available();
  int num_cpus = numa_num_task_cpus();
  int num_numa = numa_max_node() + 1;
  printf("num cpus: %d\n", num_cpus);
  printf("numa available: %d\n", num_numa);
  //numa_set_localalloc();

  /*struct bitmask *bm = numa_bitmask_alloc(num_cpus);
  for (int i=0; i<=numa_max_node(); ++i)
  {
    numa_node_to_cpus(i, bm);
    printf("numa node %d ", i);
    print_bitmask(bm);
    printf(" - %g GiB\n", numa_node_size(i, 0) / (1024.*1024*1024.));
  }
  numa_bitmask_free(bm);*/

  puts("");

  char *x;
  const size_t cache_line_size = 64;
  const size_t array_size = 100*1000*1000;
  size_t ntrips = 2;
  int numa_node = -1;

      //int numa_node = -1;
      //get_mempolicy(&numa_node, NULL, 0, (void*)x, 0);
      //printf("NUMA node: %d\n", numa_node);

#pragma omp parallel 
  {
    omp_get_num_threads() == num_cpus;
    int tid = omp_get_thread_num();
    pin_to_core(tid);
    if(tid == 0){
      //x = (char *) numa_alloc_local(array_size);
      x = (char *) numa_alloc_onnode(array_size,num_numa-1);
      x[0] = 0;
      get_mempolicy(&numa_node, NULL, 0, (void*)x, MPOL_F_NODE | MPOL_F_ADDR);
      printf("NUMA node where data is allocated: %d\n", numa_node);
      /*
     
   All numa memory allocation policy only takes effect when a
   page is actually faulted into the address space of a process
   by accessing it. The numa_alloc_* functions take care of this
   automatically. 

      int numa_node = -1;
      get_mempolicy(&numa_node, NULL, 0, (void*)x, MPOL_F_NODE | MPOL_F_ADDR);
      printf("NUMA node: %d\n", numa_node);
      */
}



    // {{{ single access
#pragma omp barrier
    for (size_t i = 0; i<num_cpus; ++i)
    {
      if (tid == i)
      {
      double t = measure_access(x, array_size, ntrips);
      printf("Tid %d from CPU %d -> bank %d : BW %g MB/s\n",tid, sched_getcpu(),numa_node,array_size*ntrips*cache_line_size / t / 1e6);
      }
#pragma omp barrier
    }
/*
#pragma omp barrier
    for (size_t i = 0; i<num_cpus; ++i)
    {
      if (tid == i)
      {
      double t = measure_access_sequential(x, array_size, ntrips);
      printf("Sequential: Tid %d from CPU %d -> bank %d : BW %g MB/s\n",tid, sched_getcpu(),numa_node,array_size*ntrips*cache_line_size / t / 1e6);
      }
#pragma omp barrier
    }

*/

  numa_free(x, array_size);
}

  return 0;

}
