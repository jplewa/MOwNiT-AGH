#include "JPlewa_AGHMatrix.h"
#include <iostream>

void ex1_addition()
{
    AGHMatrix<double> m1(5, 5, 1.2);
    AGHMatrix<double> m2(5, 5, 2.8);
    AGHMatrix<double> m3 = m1 + m2;
    //std::cout << m1 << "+" << std::endl << m2 <<  "=" << std::endl << m3 << std::endl;
    std::cout << m3 << std::endl;
}

void ex1_multiplication_1()
{
    std::vector<std::vector<double>> m1_v { { 1, 2 }, 
                                            { 3, 4 } };
    std::vector<std::vector<double>> m2_v { { 2, 0 }, 
                                            { 1, 2 } };
    AGHMatrix<double> m1(m1_v);
    AGHMatrix<double> m2(m2_v);
    AGHMatrix<double> m3 = m1 * m2;
    //std::cout << m1 << "*" << std::endl << m2 <<  "=" << std::endl << m3 << std::endl;
    std::cout << m3 << std::endl;
}

void ex1_multiplication_2()
{
    std::vector<std::vector<double>> m1_v { { 3, 4, 2 } };
    std::vector<std::vector<double>> m2_v { { 13, 9, 7, 15 },
                                            { 8, 7, 4, 6 },
                                            { 6, 4, 0, 3 } };
    AGHMatrix<double> m1(m1_v);
    AGHMatrix<double> m2(m2_v);
    AGHMatrix<double> m3 = m1 * m2;
    //std::cout << m1 << "*" << std::endl << m2 <<  "=" << std::endl << m3 << std::endl;
    std::cout << m3 << std::endl;
}

void ex1()
{
    std::cout << "EXERCISE1" << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::cout << "MATRIX ADDITION" << std::endl;
    ex1_addition();
    std::cout << "MATRIX MULTIPLICATION" << std::endl;
    ex1_multiplication_1();
    std::cout << "MATRIX MULTIPLICATION" << std::endl;
    ex1_multiplication_2();
    std::cout << "---------------------------" << std::endl;
}

void print_and_check_symmetry(std::vector<std::vector<double>> m_v) 
{
    AGHMatrix<double> m(m_v);
    std::cout << "A = " << std::endl << m;
    if (m.isSymmetric())
    {
        std::cout << "SYMMETRIC" <<  std::endl << std::endl;
    }
    else
    {
        std::cout << "NOT SYMMETRIC" <<  std::endl << std::endl;
    }
}

void ex2a()
{
    std::cout << "EXERCISE2a" << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::vector<std::vector<double>> m1_v { { 13, 9, 7, 15 },
                                           { 8, 7, 4, 6 },
                                           { 6, 4, 0, 3 } };
    print_and_check_symmetry(m1_v);

    std::vector<std::vector<double>> m2_v { { 2, 0 }, 
                                            { 1, 2 } };
    print_and_check_symmetry(m2_v);

    std::vector<std::vector<double>> m3_v { { 2, 1 }, 
                                            { 1, 7 } };
    print_and_check_symmetry(m3_v);

    std::vector<std::vector<double>> m4_v { { 2,  1,  3, 13, 8}, 
                                            { 1,  2,  7, 0, -3},
                                            { 3,  7, -7, 10, 12},
                                            { 13, 0, 10, 3, 19}, 
                                            { 8, -3, 12, 19, 100} };
    print_and_check_symmetry(m4_v);

    std::vector<std::vector<double>> m5_v { { 2,  1,  8, 13, 8}, 
                                            { 1,  2,  7, 0, -3},
                                            { 3,  7, -7, 10, 12},
                                            { 13, 0, 10, 3, 19}, 
                                            { 8, -3, 12, 19, 100} };
    print_and_check_symmetry(m5_v);
    
    std::vector<std::vector<double>> m6_v { { 13 } };
    print_and_check_symmetry(m6_v);

    std::cout << "---------------------------" << std::endl;
}

void print_and_get_determinant(std::vector<std::vector<double>> m_v)
{
    AGHMatrix<double> m(m_v);
    std::cout << "A = " << std::endl << m;
    std::cout << "det A = " << m.getDeterminant() << std::endl << std::endl;
}

void ex2b()
{
    std::cout << "EXERCISE2b" << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::vector<std::vector<double>> m1_v { { 2 } };
    print_and_get_determinant(m1_v);

    std::vector<std::vector<double>> m2_v { { 2, 0 }, 
                                            { 1, 2 } };
    print_and_get_determinant(m2_v);

    std::vector<std::vector<double>> m3_v { { -5, 0, -1 }, 
                                            {  1, 2, -1 },
                                            { -3, 4,  1} };
    print_and_get_determinant(m3_v);

    std::vector<std::vector<double>> m4_v { {  0, 1,  2, 7 }, 
                                            {  1, 2,  3, 4 },
                                            {  5, 6,  7, 8 },
                                            { -1, 1, -1, 1 } };
    print_and_get_determinant(m4_v);

    std::vector<std::vector<double>> m5_v { { 0,  1, 0, -2, 1 }, 
                                            { 1,  0, 3,  1, 1 },
                                            { 1, -1, 1,  1, 1 },
                                            { 2,  2, 1,  0, 1 },
                                            { 3,  1, 1,  1, 2 } };
    print_and_get_determinant(m5_v);
    std::cout << "---------------------------" << std::endl;
}

void print_m_and_mT(std::vector<std::vector<double>> m_v)
{
    AGHMatrix<double> m(m_v);
    std::cout << "A = " << std::endl << m;
    std::cout << "AT = " << std::endl << m.transpose() << std::endl;
}

void ex2c()
{
    std::cout << "EXERCISE2c" << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::vector<std::vector<double>> m1_v { { 2 } };
    print_m_and_mT(m1_v);

    std::vector<std::vector<double>> m2_v { { 2, 7 }, 
                                            { 1, 2 } };
    print_m_and_mT(m2_v);

    std::vector<std::vector<double>> m3_v { {  0, 1,  2, 7 }, 
                                            {  1, 2,  3, 4 },
                                            {  5, 6,  7, 8 } };
    print_m_and_mT(m3_v);

    std::vector<std::vector<double>> m4_v { { 2,  1,  3, 13, 8}, 
                                            { 1,  2,  7, 0, -3},
                                            { 3,  7, -7, 10, 12},
                                            { 13, 0, 10, 3, 19}, 
                                            { 8, -3, 12, 19, 100} };
    print_m_and_mT(m4_v);
    std::cout << "---------------------------" << std::endl;
}

void ex3()
{
    std::cout << "EXERCISE3" << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::vector<std::vector<double>> init_LU {{ 5.0, 3.0, 2.0 }, 
                                              { 1.0, 2.0, 0.0 }, 
                                              { 3.0, 0.0, 4.0 }};
    AGHMatrix<double> m(init_LU);
    std::cout << "A = " << std::endl << m << std::endl;
    std::pair<AGHMatrix<double>, AGHMatrix<double>> LU = m.luDecomposition();
    std::cout << "L = " << std::endl << LU.first << std::endl;
    std::cout << "U = " << std::endl << LU.second << std::endl;
    std::cout << "---------------------------" << std::endl;
}


void ex4()
{
    std::cout << "EXERCISE4" << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::vector<std::vector<double>> init_Cholsky {{   4,  12, -16 }, 
                                                   {  12,  37, -43 }, 
                                                   { -16, -43,  98 }};
    AGHMatrix<double> m(init_Cholsky);
    std::cout << "A = " << std::endl << m << std::endl;
    AGHMatrix<double> cholsky = m.choleskyDecomposition();
    std::cout << "L = " << std::endl << cholsky << std::endl;
    std::cout << "LT = " << std::endl << cholsky.transpose() << std::endl;
    std::cout << "---------------------------" << std::endl;
}

void ex5()
{
    std::cout << "EXERCISE5" << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::vector<std::vector<double>> init_Gauss_1 { { 0.0001, -5.0300, 5.8090, 7.8320, 9.5740 },
                                                    { 2.2660, 1.9950,  1.2120, 8.0080, 7.2190 },
                                                    { 8.8500, 5.6810,  4.5520, 1.3020, 5.7300 },
                                                    { 6.7750, -2.253,  2.9080, 3.9700, 6.2910 } };
    AGHMatrix<double> m_1(init_Gauss_1);
    std::cout << "MATRIX" << std::endl;
    std::cout << m_1;
    m_1.gaussianElimination();
    std::cout << "MATRIX AFTER ELIMINATION " << std::endl <<  m_1;

    std::cout << "RESULT" << std::endl << m_1.backwardSubstitution()  << std::endl;

    std::vector<std::vector<double>> init_Gauss_2 { {  2,  1, -1,   8 }, 
                                                    { -3, -1,  2, -11 }, 
                                                    { -2,  1,  2,  -3 } };
    AGHMatrix<double> m_2(init_Gauss_2);
    std::cout << "MATRIX" << std::endl;
    std::cout << m_2;
    m_2.gaussianElimination();
    std::cout << "MATRIX AFTER ELIMINATION" << std::endl <<  m_2;
    std::cout << "RESULT" << std::endl << m_2.backwardSubstitution() << std::endl;

    std::vector<std::vector<double>> init_Gauss_3 { { 1, 1, 1, 4 }, 
                                                    { 2, 2, 2, 8 },
                                                    { 0, 0, 5, 5 } };
    AGHMatrix<double> m_3(init_Gauss_3);
    std::cout << "MATRIX" << std::endl;
    std::cout << m_3;
    m_3.gaussianElimination();
    std::cout << "MATRIX AFTER ELIMINATION" << std::endl <<  m_3;
    std::cout << "RESULT" << std::endl << m_3.backwardSubstitution() << std::endl;
    
    std::vector<std::vector<double>> init_Gauss_4 { { 1, 1,  4 }, 
                                                    { 2, 2, 10 } };
    AGHMatrix<double> m_4(init_Gauss_4);
    std::cout << "MATRIX" << std::endl;
    std::cout << m_4;
    m_4.gaussianElimination();
    std::cout << "MATRIX AFTER ELIMINATION" << std::endl <<  m_4;
    std::cout << "RESULT" << std::endl << m_4.backwardSubstitution() << std::endl;
    std::cout << "---------------------------" << std::endl;
}

void ex6()
{
    std::cout << "---------------------------" << std::endl;
    std::cout << "EXERCISE6" << std::endl;

    std::vector<std::vector<double>> init_Jacobi_1 { { 2, 1, 11 }, 
                                                    { 5, 7, 13 } };
    AGHMatrix<double> m1(init_Jacobi_1);
    std::cout << "MATRIX" << std::endl;
    std::cout << m1 << std::endl;
    AGHMatrix<double> result1(m1.jacobiMethod());
    std::cout << "RESULT" << std::endl << result1 << std::endl << std::endl;


    std::vector<std::vector<double>> init_Jacobi_2 { { 10, -1,  2,  0,   6 }, 
                                                     { -1, 11, -1,  3,  25 },
                                                     {  2, -1, 10, -1, -11 },
                                                     {  0,  3, -1,  8,  15 } };
    AGHMatrix<double> m2(init_Jacobi_2);
    std::cout << "MATRIX" << std::endl;
    std::cout << m2 << std::endl;
    AGHMatrix<double> result2(m2.jacobiMethod());
    std::cout << "RESULT" << std::endl << result2 << std::endl;

    std::vector<std::vector<double>> init_Jacobi_3 { {  2, -1, 0, 1},
                                                     { -1, 3, -1, 8},
                                                     { 0, -1, 2, -5}} ;

    AGHMatrix<double> m3(init_Jacobi_3);
    std::cout << "MATRIX" << std::endl;
    std::cout << m3 << std::endl;
    AGHMatrix<double> result3(m3.jacobiMethod());
    std::cout << "RESULT" << std::endl << result3 << std::endl;
    std::cout << "---------------------------" << std::endl;
}

int main() 
{
    //ex1();
    //ex2a();
    //ex2b();
    //ex2c();
    //ex3();
    //ex4();
    //ex5();
    ex6();
    return 0;
}