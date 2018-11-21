#ifndef __JPLEWA_AGHMATRIX_H
#define __JPLEWA_AGHMATRIX_H

#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <numeric>
#include <algorithm>

#define K_MAX 1000
#define EPSILON 1e-10

template <typename T>
class AGHMatrix
{
  private:
    unsigned rows;
    unsigned cols;

  public:
    std::vector<std::vector<T>> matrix;

    AGHMatrix(const std::vector<std::vector<T>> &matrix);
    AGHMatrix(unsigned _rows, unsigned _cols, const T &_initial);
    AGHMatrix(const AGHMatrix<T> &rhs);
    virtual ~AGHMatrix() = default;

    // Operator overloading, for "standard" mathematical matrix operations
    AGHMatrix<T> &operator=(const AGHMatrix<T> &rhs);

    // Matrix mathematical operations
    AGHMatrix<T> operator+(const AGHMatrix<T> &rhs);
    AGHMatrix<T> operator*(const AGHMatrix<T> &rhs);

    // Access the individual elements
    T &operator()(const unsigned &row, const unsigned &col);
    const T &operator()(const unsigned &row, const unsigned &col) const;

    // Printing matrix
    std::ostream &operator<<(const AGHMatrix<T> &matrix);

    // Access the row and column sizes
    unsigned get_rows() const;
    unsigned get_cols() const;


    AGHMatrix<T> jacobiMethod(int equation_id, std::ostream &convergence_file);
    AGHMatrix<T> gaussSeidelMethod(int equation_id, std::ostream &convergence_file);
    AGHMatrix<T> sorMethod(T omega, int equation_id, std::ostream &convergence_file);
};

#include "JPlewa_AGHMatrix.cpp"

#endif