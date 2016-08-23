#if !defined(__print_array_h)
#define __print_array_h

#include<iostream>

// For printing C-type vector to std::ostream
template <class T> struct W {
  int n;
  const T *a;
  const char *sep;
  W(const int &n, const T *a, const char *sep) : n(n), a(a), sep(sep) {}
};


template <class T> std::ostream & operator << (std::ostream& os, const W<T> & w){
  for (int i = 0; i < w.n-1; ++i) os << w.a[i] << w.sep;
  return (os << w.a[w.n-1]);
}

#endif
