#pragma once

/*
  Copy array -- simple but vectorizable
*/
template<class T> void cp(T* dest, const T* src, const int &dim) {
  for (int i = 0; i < dim; ++i) dest[i] = src[i];
}

/*
  Zero array -- simple but vectorizable
*/
template<class T> void zero(T* a, const int &dim) {
  for (int i = 0; i < dim; ++i) a[i] = T(0);
}



/*
  Library for allocating matrices of different shapes

  Author: Martin Horvat, March 2021
*/

/*
  Allocate a C-style matrix m[nrow][ncol]

  Input:
    nrow - number of rows,
    ncol - number of columns

  Return:
    m
*/
template <class T> T** matrix(int nrow, int ncol) {

  T **m = new T* [nrow];

  m[0] = (T *) new char [nrow*ncol*sizeof(T)];

  for( int i = 1; i < nrow; i++) m[i] = m[i-1]+ ncol;

  return m;
}


/*
  Free a C-style matrix.

  Input:
    m - C-style matrike

  Output:
    m = 0, m is put to zero
*/
template <class T> void free_matrix(T **&m) {

  delete [] (char*) m[0];

  delete [] m;

  m = 0;
}

// ============================================================================

/*
  Allocate C-style tensor of 3. order m[n1][n2][n3]

  Input:
    n1, n2, n3 - dimensions

  Return
    m
*/
template <class T> T*** tensor3 (const int &n1, const int &n2, const int &n3){

  int i, j;

  T ***m = new T** [n1], *p = new T[n1*n2*n3];

  for (i = 0; i < n1; ++i) {

    m[i] = new T* [n2];

    for (j = 0; j < n2; ++j){
      m[i][j] = p;
      p += n3;
    }
  }

  return m;
}

/*
  Free C-style tensor of 3.oder

  Input:
    m - C-style matrike
    n1, n2, n3 - dimensions, n2, n3 are optional

  Output:
    m = 0, m is put to zero
*/
template <class T> void free_tensor3(T ***&m, const int &n1, const int &n2=0, const int &n3=0){

  delete [] m[0][0];

  for (int i = 0; i < n1; ++i) delete [] m[i];

  delete [] m;

  m = 0;
}

// ============================================================================

/*
  Reserve space of triangle C-style matrix A for degrees n having elements

    A[0][0]
    A[1][0] A[1][1]
    ...
    A[n][0] A[n][1] .... A[n][n]

  Input:

    n - degree ( >= 0)

  Return:
    A - C-style triangle matrix
*/
template <class T>
T **tmatrix (int n){

  int size = (n + 1)*(n + 2)/2;

  T **m = new T* [n + 1];

  *m = new T [size];

  for (int i = 1; i <= n; ++i) m[i] = m[i-1] + i;

  return m;
}

/*
  Free the space of the triangle matrix

  Input:
    A - C-style triangle matrix

  Output:
    A  is put to zero
*/
template <class T>
void free_tmatrix(T ** & m){

  delete [] *m;
  delete [] m;

  m = 0;
}

// ============================================================================

/*
  Reserve space of triangle C-style matrix A for degrees n in [nmin,nmax]
  having elements

    A[nmin][0] ... A[nmin][nmin]
    ...
    A[n][0] A[n][1] .... A[n][n]
    ...
    A[nmax][0] A[nmax][1] .... A[nmax][max]

  Input:

    nmin, nmax - degree ( >= 0)

  Return:
    A - C-style triangle matrix
*/
template <class T>
T **tmatrix2 (int nmin, int nmax){

  int size = ((nmax + 1)*(nmax + 2) - nmin*(nmin + 1)) >> 1;

  T **m = new T* [nmax - nmin + 1];

  *m = new T [size];

  m -= nmin;

  for (int i = nmin + 1; i <= nmax; ++i) m[i] = m[i-1] + i;

  return m;
}

/*
  Free the space of the triangle matrix

  Input:
    A - C-style triangle matrix

  Output:
    A  is put to zero
*/
template <class T>
void free_tmatrix2(T ** & m, int nmin, int nmax){

  m += nmin;

  delete [] *m;
  delete [] m;

  m = 0;
}
