#include "JPlewa_AGHMatrix.h"
#include "gtest/gtest.h"

TEST(LinearEquationSolver, overdeterminedSystem) { 
    std::vector<std::vector<double>> matrix_vals { { 0.0001, -5.0300, 5.8090, 7.8320 },
                                                   { 2.2660, 1.9950,  1.2120, 8.0080 },
                                                   { 8.8500, 5.6810,  4.5520, 1.3020 },
                                                   { 6.7750, -2.253,  2.9080, 3.9700 } };
    AGHMatrix<double> matrix(matrix_vals);
    try
    {
        matrix.gaussianElimination();
        FAIL();
    }
    catch(const exceptions::OverdeterminedSystemException& err)
    {
        ASSERT_STREQ("The system of equations is overdetermined" , err.what());
    }
}

TEST(LinearEquationSolver, underdeterminedSystem) { 
    std::vector<std::vector<double>> matrix_vals { { 0.0001, -5.0300, 5.8090, 7.8320 },
                                                   { 2.2660, 1.9950,  1.2120, 8.0080 } };
    AGHMatrix<double> matrix(matrix_vals);
    try
    {
        matrix.gaussianElimination();
        FAIL();
    }
    catch(const exceptions::UnderdeterminedSystemException& err)
    {
        ASSERT_STREQ("The system of equations is underdetermined" , err.what());
    }
}

TEST(LinearEquationSolver, uniqueSolution1) { 
    std::vector<double> expected_result { 0.21602477, -0.00791511,  0.63524333, 0.74617428 };

    std::vector<std::vector<double>> matrix_vals {{ 0.0001, -5.0300, 5.8090, 7.8320, 9.5740 },
                                                  { 2.2660, 1.9950,  1.2120, 8.0080, 7.2190 },
                                                  { 8.8500, 5.6810,  4.5520, 1.3020, 5.7300 },
                                                  { 6.7750, -2.253,  2.9080, 3.9700, 6.2910 }};
    AGHMatrix<double> matrix(matrix_vals);
    matrix.gaussianElimination();
    std::vector<double> result(matrix.backwardSubstitution().transpose().matrix[0]);
    
    EXPECT_EQ(result.size(), expected_result.size());
    for (int i = 0; i < result.size(); i++){
        EXPECT_NEAR(expected_result[i], result[i], 1e-8);
    }
}

TEST(LinearEquationSolver, uniqueSolution2) { 
    std::vector<double> expected_result { 2, 3, -1 };

    std::vector<std::vector<double>> matrix_vals  {{  2,  1, -1,   8 }, 
                                                   { -3, -1,  2, -11 }, 
                                                   { -2,  1,  2,  -3 }};  
    AGHMatrix<double> matrix(matrix_vals);
    matrix.gaussianElimination();
    std::vector<double> result(matrix.backwardSubstitution().transpose().matrix[0]);
    
    EXPECT_EQ(result.size(), expected_result.size());
    for (int i = 0; i < result.size(); i++){
        EXPECT_NEAR(expected_result[i], result[i], 1e-8);
    }
}

TEST(LinearEquationSolver, uniqueSolutionPivotingRequired) { 
    std::vector<double> expected_result { -10, 4, 11 };

    std::vector<std::vector<double>> matrix_vals  {{  1,  -1,  2,   8 }, 
                                                   {  0,  0,  -1, -11 }, 
                                                   { 0,  2,  -1,  -3 }};  
    AGHMatrix<double> matrix(matrix_vals);
    matrix.gaussianElimination();
    std::vector<double> result(matrix.backwardSubstitution().transpose().matrix[0]);
    
    EXPECT_EQ(result.size(), expected_result.size());
    for (int i = 0; i < result.size(); i++){
        EXPECT_NEAR(expected_result[i], result[i], 1e-8);
    }
}

TEST(LinearEquationSolver, infinitelyManySolutions) { 
    std::vector<std::vector<double>> matrix_vals  {{  2,  1, -1,  8 }, 
                                                   { 2,  1, -1,  8 }, 
                                                   { -2,  1,  2, -3 }};
    AGHMatrix<double> matrix(matrix_vals);
    matrix.gaussianElimination();
    std::vector<double> result(matrix.backwardSubstitution().transpose().matrix[0]);
    
    EXPECT_TRUE(isnan(result[result.size() - 1]));
}

TEST(LinearEquationSolver, noSolutions) { 
    std::vector<std::vector<double>> matrix_vals  {{  2,  1, -1,  8 }, 
                                                   {  2,  1, -1,  7 }, 
                                                   { -2,  1,  2, -3 }};
    
    AGHMatrix<double> matrix(matrix_vals);
    matrix.gaussianElimination();
    std::vector<double> result(matrix.backwardSubstitution().transpose().matrix[0]);
    
    EXPECT_TRUE(isinf(result[result.size() - 1]));
}

  
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}







