#ifndef __JPLEWA_AGHMATRIX_H
#define __JPLEWA_AGHMATRIX_H

#include <vector>
#include <iostream>
#include <math.h>
#include <numeric>
#include <algorithm>

template <typename T>
class AGHMatrix
{
  private:
    unsigned rows;
    unsigned cols;
    T laplaceExpansion(int start_col, std::vector<int> rows);
    void swap_rows(int h, int i_max);

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

    bool isSymmetric();
    T getDeterminant();
    AGHMatrix<T> transpose();

    std::pair<AGHMatrix<T>, AGHMatrix<T>> luDecomposition();
    AGHMatrix<T> choleskyDecomposition();

    void gaussianElimination();
    AGHMatrix<T> backwardSubstitution();
    AGHMatrix<T> forwardSubstitution();

    AGHMatrix<T> jacobiMethod();
};

#include "JPlewa_AGHMatrix.cpp"

#endif