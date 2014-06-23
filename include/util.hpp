#ifndef UTIL_HPP__
#define UTIL_HPP__

#include <execinfo.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <algorithm>

//Calling this macro prints a stack backtrace of N lines
#define BACKTRACE(N)						\
  {								\
    void* tracePtrs[N];						\
    int count = backtrace(tracePtrs, N);			\
    char** funcNames = backtrace_symbols(tracePtrs, count);	\
    printf("Backtrace:\n");					\
    if(funcNames != NULL)					\
      for(int i = 0; i < count; i++)				\
	if(funcNames[i][0] != '\0')				\
	  printf("%s\n", funcNames[i]);				\
    free(funcNames);						\
  }

//Asserts that 's' is true, otherwise prints msg..., some runtime information (source of error and backtrace) and exits
#define ASSERT(s, msg...)						\
  {									\
    if(!(s)) {								\
      printf("Failed to assert '%s' in %s on line %d\n\n", #s, __FILE__, __LINE__); \
      printf(msg);							\
      printf("\n");							\
      BACKTRACE(20);							\
      fflush(stdout);							\
      raise(SIGINT);							\
    }									\
  }

/*//Isn't this a part of std::?
template <class T>
void swap(T& a, T& b) {
  T tmp = a;
  a = b;
  b = tmp;
}*/

//This enforces periodic boundary conditions on n
//  n == -1 is the same as n == N - 1
inline double adjust(double n, double N)
{
  if(n < 0)
    {
      n = n + (floor(-n / N) + 1) * N;
    }

  if(n >= N)
    {
      n = fmod(n, N);
    }

  return n;
}

/*inline Index adjust(Index q, int N)
{
  return Index(adjust(q.getn(), N), adjust(q.getm(), N), q.getl());
  }*/

/*template<class T>
inline typename T::iterator select(double rn, T &array)
{
  double sum = 0.0;
  for(typename T::iterator it = array.begin(); it != array.end(); it++)
    sum += (*it).first;
  
  double v = rn * sum;
  double tsum = 0.0;

  for(typename T::iterator it = array.begin(); it != array.end(); it++)
    {
      tsum += (*it).first;
      if(tsum >= v)
	{
	  if(DEBUG)
	    std::cout << "findR: Selecting event, rn = " << rn << ", v = " << v << ", sum = " << sum << std::endl;

	  return it;
	}
    }
  
  return array.end();
  }*/

 /*inline int select(double rn, double *array, int size)
{
  double sum = 0.0;
  for(int i = 0; i < size; i++)
    sum += array[i];
  
  double v = rn * sum;
  sum = 0.0;
  int val = 0;
  for(int i = 0; i < size; i++)
    {
      sum += array[i];
      if(sum >= v)
	{
	  val = i;
	  break;
	}
    }

  return val;
  }*/

template<class T, class U>
U copy_it(T itBegin, T itEnd, U out)
{
  for(; itBegin != itEnd; itBegin++, out++)
    {
      *out = itBegin;
    }

  return out;
}

template <class T>
void align_allocate(T **ptr, void **memory, size_t count)
{
  uint64_t mptr = (uint64_t)*memory;

  *ptr = (T *)((mptr + 255ull) & ~255ull);

  *memory = *ptr + count;
}
 
//This computes a periodic distance [0, N - 1]
template <class T>
inline T cyclicDistance(T n1, T n2, T N)
{
  /*if(std::abs(n2 - n1) > N / 2)
    {
      return std::abs(N - std::(
    }
  else
    {
      return std::abs(n2 - n1);
      }*/

  if(n2 > n1)
    {
      if(n2 - n1 < (N - n2 + n1))
	{
	  return n1 - n2;
	}
      else
	{
	  return (N - n2 + n1);
	}
    }
  else
    {
      if(n1 - n2 < (N - n1 + n2))
	{
	  return n1 - n2;
	}
      else
	{
	  return -(N - n1 + n2);
	}
    }
  //return std::min(std::min(std::abs(n2 - n1), std::abs((N - n1) + (n2 - 0))), std::abs((N - n2) + (n1 - 0)));
}

/*template <>
inline Index cyclicDistance<Index>(Index n1, Index n2, int N)
{
  return Index(cyclicDistance<int>(n1.getn(), n2.getn(), N),
	       cyclicDistance<int>(n1.getm(), n2.getm(), N), std::abs(n1.getl() - n2.getl()));
               }*/

#endif
