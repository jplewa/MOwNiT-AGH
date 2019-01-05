#include "JPlewa_AGHMatrix.h"

std::vector<std::vector<std::vector<double>>> MATRICES
{ 
    { { 2, 1, 11 }, 
      { 5, 7, 13 } },

    { { 4,  1, 3, 17 },
      { 1,  5, 1, 14 },
      { 2, -1, 8, 12 } }, 

    { { 12,  3, -5,  1 },
      {  1,  5,  3, 28 },
      {  3,  7, 13, 76 } }, 

    { { 3.0, -0.1, -0.2,  7.85 },
      { 0.1,    7, -0.3, -19.3 },
      { 0.3, -0.2,   10,  71.4 } }, 

    { { 10, -1,  2,  0,   6 }, 
      { -1, 11, -1,  3,  25 },
      {  2, -1, 10, -1, -11 },
      {  0,  3, -1,  8,  15 } },

    { { 10,  -1,  2,  0,   6,    1 }, 
      { -1,  11, -1,  3,   5,   -5 },
      {  2,  -1, 15, -1, -11,   10 },
      {  0,   3, -1,  8, 5.5,  3.5 },
      { -4, 7.5, -9,  3,  17, -0.8 } },

    { { 10,  -1,  2,  0,   6,    1, 1 }, 
      { -1,  11, -1,  3,   5,   -5, 2 },
      {  2,  -1, 15, -1, -11,   10, 3 },
      {  0,   3, -1,  8, 5.5,  3.5, 4 },
      { -4, 7.5, -9,  3,  17, -0.8, 5 },
      {  0,   1,  2,  3,   4,   19, 6 } },

    { { 10,  -1,  2,  0,   6,    1,  1, 6 }, 
      { -1,  11, -1,  3,   5,   -5,  2, 5 },
      {  2,  -1, 15, -1, -11,   10,  3, 4 },
      {  0,   3, -1,  8, 5.5,  3.5,  4, 3 },
      { -4, 7.5, -9,  3,  17, -0.8,  5, 2 },
      {  0,   1,  2,  3,   4,   19,  6, 1 },
      {  6,   5,  4,  3,   2,    1, 10, 0 } },

};

void ex1()
{
    for (int i = 0; i < MATRICES.size(); i++)
    {
      AGHMatrix<double> J1(MATRICES[i]);
      std::cout << J1.jacobiMethod(-1, std::cout) << std::endl;
    }
}

void ex2()
{
    for (int i = 0; i < MATRICES.size(); i++)
    {
      AGHMatrix<double> J1(MATRICES[i]);
      std::cout << J1.gaussSeidelMethod(-1, std::cout) << std::endl;
    }
}

void ex3()
{
    for (int i = 0; i < MATRICES.size(); i++)
    {
      AGHMatrix<double> J1(MATRICES[i]);
      std::cout << J1.sorMethod(1.1, -1, std::cout) << std::endl;
    }
}

void ex5()
{
    std::ofstream jacobi;
    jacobi.open("csv/ex5_jacobi.csv");

    std::ofstream gauss_seidel;
    gauss_seidel.open("csv/ex5_gauss_seidel.csv");

    std::ofstream sor;
    sor.open("csv/ex5_sor.csv");


    for (int i = 0; i < MATRICES.size(); i++)
    {
      AGHMatrix<double> J1(MATRICES[i]);
      std::cout << J1.jacobiMethod(i, jacobi) << std::endl;
      std::cout << J1.gaussSeidelMethod(i, gauss_seidel) << std::endl;
      std::cout << J1.sorMethod(1.1, i, sor) << std::endl;
    }

    jacobi.close();
    gauss_seidel.close();
    sor.close();
}

int main() 
{
    ex1();
    ex2();
    ex3();
    ex5();
    return 0;
}