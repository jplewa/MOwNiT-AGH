import numpy as np

def printSolution(A, b):
    x = np.linalg.solve(A, b)
    np.allclose(np.dot(A, x), b)
    print(x.transpose())
    print('------------------------------')

def uniqueSolution1():
    print('uniqueSolution1')
    A = np.matrix([[0.0001, -5.0300, 5.8090, 7.8320],
                    [2.2660, 1.9950,  1.2120, 8.0080],
                    [8.8500, 5.6810,  4.5520, 1.3020],
                    [6.7750, -2.253,  2.9080, 3.9700]])

    b = np.matrix([9.5740, 7.2190, 5.7300, 6.2910]).transpose()
    printSolution(A, b)

def uniqueSolution2():
    print('uniqueSolution2')   
    A = np.matrix([[  2,  1, -1], 
                   [ -3, -1,  2],  
                   [ -2,  1,  2]])
    b = np.matrix([8, -11, -3]).transpose()
    printSolution(A, b)
    

def uniqueSolutionPivotingRequired():
    print('uniqueSolutionPivotingRequired')   
    A = np.matrix([[ 1,  -1,  2 ],
                   [ 0,  0,  -1 ], 
                   [ 0,  2,  -1 ]])
    b = np.matrix([8, -11, -3]).transpose()
    printSolution(A, b)


uniqueSolution1()
uniqueSolution2()
uniqueSolutionPivotingRequired()
