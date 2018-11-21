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
AGHMatrix<T> AGHMatrix<T>::jacobiMethod(int equation_id, std::ostream &convergence_file)
{
  // result vector - the initial estimation is 0
  AGHMatrix<T> x(get_rows(), 1, 0);
  // the next estimation of the result
  std::vector<double> y(get_rows(), 0);
  // K_MAX denotes the maximum number of iterations
  for (int k = 0; k < K_MAX; k++)
  {
    // copy current result estimations into a temporary vector
    for (int i = 0; i < get_rows(); i++)
    {
      y[i] = x.matrix[i][0];
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
      x.matrix[i][0] = sum / diag;
    }
    T norm = 0;
    T max_y = 0;
    for (int i = 0; i < get_rows(); i++)
    {
      T diff = fabs(y[i] - x.matrix[i][0]);
      if (diff > norm)
      {
        norm = diff;
      }
      if (y[i] > max_y)
      {
        max_y = y[i];
      }
    }
    if (equation_id != -1)
    {
      convergence_file << equation_id << ", " << k << ", " << norm << std::endl;
    } 
    // return if convergence was reached
    if (norm < EPSILON * max_y)
    {
        std::cout << "CONVERGENCE REACHED IN ITERATION: " << k << std::endl;
        return x;
    }
  }
  std::cout << "ITERATION LIMIT REACHED WITH NO CONVERGENCE" << std::endl;
  return x;
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::sorMethod(T omega, int equation_id, std::ostream &convergence_file)
{
  // result vector - the initial estimation is 0
  AGHMatrix<T> x(get_rows(), 1, 0);
  // the next estimation of the result
  std::vector<double> y(get_rows(), 0);
  // K_MAX denotes the maximum number of iterations
  for (int k = 0; k < K_MAX; k++)
  {
    // copy current result estimations into a temporary vector
    for (int i = 0; i < get_rows(); i++)
    {
      y[i] = x.matrix[i][0];
    }
    // calculate a new estimation of the variables based on the previous one
    for (int i = 0; i < get_rows(); i++)
    {
      T sum = matrix[i][get_cols() - 1];  
      T diag = matrix[i][i];
      for (int j = 0; j < i; j++)
      {
        sum -= matrix[i][j] * x.matrix[j][0];
      }
      for (int j = i + 1; j < get_rows(); j++)
      {
        sum -= matrix[i][j] * x.matrix[j][0]; 
      }
      x.matrix[i][0] = sum / diag;
      x.matrix[i][0] = omega * x.matrix[i][0] + (1 - omega) * y[i];
    }
    T norm = 0;
    T max_y = 0;
    for (int i = 0; i < get_rows(); i++)
    {
      T diff = fabs(y[i] - x.matrix[i][0]);
      if (diff > norm)
      {
        norm = diff;
      }
      if (y[i] > max_y)
      {
        max_y = y[i];
      }
    }
    if (equation_id != -1)
    {
      convergence_file << equation_id << ", " << k << ", " << norm << std::endl;
    } 
    // return if convergence was reached
    if (norm < EPSILON * max_y)
    {
        std::cout << "CONVERGENCE REACHED IN ITERATION: " << k << std::endl;
        return x;
    }
  }
  std::cout << "ITERATION LIMIT REACHED WITH NO CONVERGENCE" << std::endl;
  return x;
}

template<typename T>
AGHMatrix<T> AGHMatrix<T>::gaussSeidelMethod(int equation_id, std::ostream &convergence_file){
  return sorMethod(1.0, equation_id, convergence_file);
}