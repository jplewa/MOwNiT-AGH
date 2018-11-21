#include "JPlewa_AGHMatrix.h"

namespace exceptions
{
	class MismatchedMatricesException : public std::logic_error
	{
	  public:
		  MismatchedMatricesException() : std::logic_error{ "Matrices have mismatched dimensions" } {}
	};

  class MatrixNotSquareException : public std::logic_error
	{
	  public:
		  MatrixNotSquareException() : std::logic_error{ "Matrix must be square" } {}
	};

  class UnderdeterminedSystemException : public std::logic_error
	{
	  public:
		  UnderdeterminedSystemException() : std::logic_error{ "The system of equations is underdetermined" } {}
	};
  class OverdeterminedSystemException : public std::logic_error
	{
	  public:
		  OverdeterminedSystemException() : std::logic_error{ "The system of equations is overdetermined" } {}
	};
}

// Parameter Constructor                                                                                                                                                      
template<typename T>
AGHMatrix<T>::AGHMatrix(const std::vector<std::vector<T>>& mat) 
{
  matrix.resize(mat.size());
  for (unsigned i = 0; i < mat.size(); i++) 
  {
    matrix[i].resize(mat[i].size());
    for(unsigned j = 0; j < mat[i].size(); j++)
    {
      matrix[i][j] = mat[i][j];
    }
  }
  rows = matrix.size();
  cols = matrix[0].size();
}

// Parameter Constructor                                                                                                                                                      
template<typename T>
AGHMatrix<T>::AGHMatrix(unsigned _rows, unsigned _cols, const T& _initial) 
{
  matrix.resize(_rows);
  for (unsigned i=0; i<matrix.size(); i++) 
  {
    matrix[i].resize(_cols, _initial);
  }
  rows = _rows;
  cols = _cols;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
AGHMatrix<T>::AGHMatrix(const AGHMatrix<T>& rhs) 
{
  matrix = rhs.matrix;
  rows = rhs.get_rows();
  cols = rhs.get_cols();
}

// Get the number of rows of the matrix                                                                                                                                       
template<typename T>
unsigned AGHMatrix<T>::get_rows() const 
{
  return this->rows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<typename T>
unsigned AGHMatrix<T>::get_cols() const 
{
  return this->cols;
}

// Assignment Operator                                                                                                                                                        
template<typename T>
AGHMatrix<T>& AGHMatrix<T>::operator=(const AGHMatrix<T>& rhs) 
{
  if (&rhs == this)
    return *this;

  unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols();

  matrix.resize(new_rows);
  for (unsigned i = 0; i < matrix.size(); i++) 
  {
    matrix[i].resize(new_cols);
  }

  for (unsigned i = 0; i < new_rows; i++) 
  {
    for (unsigned j = 0; j < new_cols; j++) 
    {
      matrix[i][j] = rhs(i, j);
    }
  }
  rows = new_rows;
  cols = new_cols;

  return *this;
}

// Access the individual elements                                                                                                                                             
template<typename T>
T& AGHMatrix<T>::operator()(const unsigned& row, const unsigned& col) 
{
  return this->matrix[row][col];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& AGHMatrix<T>::operator()(const unsigned& row, const unsigned& col) const 
{
  return this->matrix[row][col];
}

// Addition of two matrices                                                                                                                                                   
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator+(const AGHMatrix<T>& rhs) 
{
  if ((get_cols() != rhs.cols) || (get_rows() != rhs.rows))
  {
  			throw exceptions::MismatchedMatricesException();
  }
  AGHMatrix<T> result(*this);
  for (int i = 0; i < get_rows(); i++)
  {
    for (int j = 0; j < get_cols(); j++){
      result.matrix[i][j] = matrix[i][j] + rhs.matrix[i][j];
    }
  }
  return(result);
}

// Left multiplication of this matrix and another                                                                                                                              
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator*(const AGHMatrix<T>& rhs) 
{
  // the number of columns in the first matrix needs to be equal to the number of rows in the second one
    if ((get_cols() != rhs.rows))
    {
    			throw exceptions::MismatchedMatricesException();
    }
    // the resulting matrix has the same number of rows as the first one and the same number of columns as the second one
    AGHMatrix<T> result(get_rows(), rhs.cols, 0);

    for (int i = 0; i < get_rows(); i++)
    {
      for (int j = 0; j < rhs.cols; j++)
      {
        for (int k = 0; k < get_cols(); k++)
        {
          result.matrix[i][j] += matrix[i][k] * rhs.matrix[k][j];
        }
      }
    }
    return(result);
}

// Printing matrix                                                                                                                        
template<typename T>
std::ostream& operator<<(std::ostream& stream, const AGHMatrix<T>& matrix) 
{
  for (int i = 0; i < matrix.get_rows(); i++) 
  { 
    for (int j = 0; j < matrix.get_cols(); j++) 
    {
        stream << matrix(i,j) << ", ";
    }
    stream << std::endl;
  }
  stream << std::flush;
}


template<typename T>
bool AGHMatrix<T>::isSymmetric()
{
  if (get_cols() != get_rows())
  {
    return false;  
  }
  for (int i = 0; i < get_rows(); i++)
  {
    for (int j = 0; j < i; j++)
    {
      if (matrix[i][j] != matrix[j][i])
      {
        return false;
      }
    }
  }
  return true;
}

// Laplace expansion of the matrix along the first column
template<typename T>
T AGHMatrix<T>::laplaceExpansion(int start_col, std::vector<int> rows)
{ 
  if (start_col + 1 == get_cols() && rows.size() == 1)
  {
    return matrix[rows[0]][start_col];
  }
  T det = 0;
  for (int i = 0; i < rows.size(); i++)
  {
    std::vector<int> new_rows(rows);
    new_rows.erase(std::remove(new_rows.begin(), new_rows.end(), rows[i]), new_rows.end());
    det += matrix[rows[i]][start_col] * pow(-1, i + 2) * laplaceExpansion(start_col + 1, new_rows);
  }
  return det;
}

template<typename T>
T AGHMatrix<T>::getDeterminant()
{
  // the matrix needs to be square
  if (get_cols() != get_rows())
  {
    throw exceptions::MatrixNotSquareException();
  }
  // to avoid having to copy matrix values we shall always use Laplace expansion along the first column
  // this way we need to remember the column we're starting from and the rows that compose the matrix in question
  std::vector<int> rows(get_rows());
  std::iota(rows.begin() + 1, rows.end(), 1);
  return laplaceExpansion(0, rows);
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::transpose()
{
  AGHMatrix<T> transposed(get_cols(), get_rows(), 0);
  for (int i = 0; i < get_rows(); i++)
  {
    for (int j = 0; j < get_cols(); j++)
    {
      (transposed.matrix)[j][i] = matrix[i][j];
    }
  } 
  return transposed;
}

template<typename T>
std::pair<AGHMatrix<T>, AGHMatrix<T>> AGHMatrix<T>::luDecomposition()
{
  if (get_cols() != get_rows()){
    throw exceptions::MatrixNotSquareException();
  }
  T tmp_U = 0;
  T tmp_L = 0;
  AGHMatrix<T> U(get_rows(), get_cols(), 0);
  AGHMatrix<T> L(get_rows(), get_cols(), 0);
  for (int i = 0; i < get_rows(); i++)
  {
    L.matrix[i][i] = 1;
    tmp_U = 0;
    for (int k = 0; k < i; k++)
    {
      tmp_U += L.matrix[i][k] * U.matrix[k][i];
    }
    U.matrix[i][i] = matrix[i][i] - tmp_U;
    for (int j = i + 1; j < get_rows(); j++)
    {
      tmp_U = 0;
      tmp_L = 0;
      for (int k = 0; k < i; k++)
      {
        tmp_U += L.matrix[i][k] * U.matrix[k][j];
        tmp_L += L.matrix[j][k] * U.matrix[k][i];
      }
      U.matrix[i][j] = matrix[i][j] - tmp_U;
      L.matrix[j][i] = (matrix[j][i] - tmp_L) / U.matrix[i][i];
    }
  }
  return std::make_pair(L, U);
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::choleskyDecomposition()
{
  if (get_cols() != get_rows()){
    throw exceptions::MatrixNotSquareException();
  }
  T tmp_L;
  AGHMatrix<T> L(get_rows(), get_cols(), 0);
  for (int i = 0; i < get_rows(); i++){
    tmp_L = 0;
    for (int k = 0; k < i; k++)
    {
      tmp_L += pow(L.matrix[i][k], 2);
    }
    L.matrix[i][i] = sqrt(matrix[i][i] - tmp_L);
    for (int j = i + 1; j < get_rows(); j++)
    {
      tmp_L = 0;
      for (int k = 0; k < i; k++)
      {
        tmp_L = L.matrix[j][k] * L.matrix[i][k];
      }
      L.matrix[j][i] = (matrix[j][i] - tmp_L) / L.matrix[i][i];
    }
  }
  return L;
}

// swap matrix rows
template<typename T>
void AGHMatrix<T>::swap_rows(int h, int i_max)
{
  T tmp;
  for (int i = 0; i <  get_cols(); i++){
    tmp = matrix[h][i];
    matrix[h][i] = matrix[i_max][i];
    matrix[i_max][i] = tmp;
  }
}

// transform matrix into row echelon form 
template<typename T>
void AGHMatrix<T>::gaussianElimination()
{
  // underdetermined systems have no solutions
  if (get_cols() > (get_rows()) + 1)
  {
      throw exceptions::UnderdeterminedSystemException();
  }
  // overdetermined systems might have a solution
  else if (get_cols() != (get_rows()) + 1)
  {
      throw exceptions::OverdeterminedSystemException();
  }
  int h = 0;  // pivot row
  int k = 0;  // pivot column
  while(h < get_rows() && k < get_cols())
  {
    // looking for the k-th partial pivot
    int i_max = h;
    for (int i = h + 1; i < get_rows(); i++)
    {
      if (abs(matrix[i][k]) > abs(matrix[i_max][k])){
        i_max = i;
      }
    }
    // if no pivot is found in the k-th column, go to the next column
    if (matrix[i_max][k] == 0)
    {
      k += 1;
    }
    else
    {
      swap_rows(h, i_max);
      // all rows below the pivot
      for (int i = h + 1; i < get_rows(); i++)
      {
        T coeff = matrix[i][k] / matrix[h][k];
        // fill the lower part of pivot column with zeros
        matrix[i][k] = 0;
        // continue going down the current row
        for (int j = k + 1; j < get_cols(); j++)
        {
          matrix[i][j] -= matrix[h][j] * coeff;
        }
      }
      h += 1;
      k += 1;
    }
  }
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::backwardSubstitution()
{
  AGHMatrix<T> results(get_rows(), 1, 0);
  for (int i = get_rows() - 1; i >= 0; i--)
  {
    results.matrix[i][0] = matrix[i][get_cols() - 1];
    for (int j = get_cols() - 2; j > i; j--)
    {
      results.matrix[i][0] -= matrix[i][j] * results.matrix[j][0];
    }
    // nan -> infinite solutions
    // inf -> no solutions
    results.matrix[i][0] /= matrix[i][i];
  }
  return results;
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::forwardSubstitution()
{
  AGHMatrix<T> results(get_rows(), 1, 0);
  for (int i = 0; i < get_rows(); i++)
  {
    results.matrix[i][0] = matrix[i][get_cols() - 1];
    for (int j = 0; j < i; j++)
    {
      results.matrix[i][0] -= matrix[i][j] * results.matrix[j][0];
    }
    // nan -> infinite solutions
    // inf -> no solutions
    results.matrix[i][0] /= matrix[i][i];
  }
  return results;
}

#define K_MAX 100
#define EPSILON 1e-7

template<typename T>
AGHMatrix<T> AGHMatrix<T>::jacobiMethod()
{
  // result vector - the initial estimation is 0
  AGHMatrix<T> results(get_rows(), 1, 0);
  // the next estimation of the result
  std::vector<double> y(get_rows(), 0);
  // K_MAX denotes the maximum number of iterations
  for (int k = 0; k < K_MAX; k++)
  {
    // copy current result estimations into a temporary vector
    for (int i = 0; i < get_rows(); i++)
    {
      y[i] = results.matrix[i][0];
    }
    // calculate a new estimation of the variables based on the previous one
    for (int i = 0; i < get_rows(); i++)
    {
      T sum = matrix[i][get_cols() - 1];  
      T diag = matrix[i][i];
      for (int j = 0; j < get_rows(); j++)
      {
        if (j != i)
        {
          sum -= matrix[i][j] * y[j];
        }
      }
      results.matrix[i][0] = sum / diag;
    }
    std::cout << "ITERATION #" << k << std::endl << results;
    double norm = 0;
    for (int i = 0; i < get_rows(); i++)
    {
      if (fabs(y[i] - results.matrix[i][0]) > norm)
      {
        norm = fabs(y[i] - results.matrix[i][0]);
      }
    }
    std::cout << "ERR: " << norm << std::endl << std::endl;
    // return if convergence was reached
    if (norm < EPSILON)
    {
        return results;
    }
  }
  std::cout << "ITERATION LIMIT REACHED" << std::endl;
  return results;
}

