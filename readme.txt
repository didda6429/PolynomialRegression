##-----------------------------------
## README - Polynomial Regression
##-----------------------------------

PART 1 - Introduction
-----------------------------------
This module implements a polynomial regression algorithm which works for polynomials
of arbitrary order. It is suited for use in basic machine learning and statistics applications.

Note that I have created a method called testMain which acts as a de facto main method,
whilst still maintaining that it is a module rather than a standalone program.

This was part of my coursework for the year 1 imperative programming module at the 
University of Bristol.

PART 2 - Usage
-----------------------------------
The module only exposes the calculateCoefficients function as the other functions
are implementation specific and are only there because they are needed by the calculateCoefficients
function.

The calculateCoefficients functions expects 5 inputs:
  --> int n:                               the number of x and y data points being input into the algorithm.
                                           This must be at least equal to the desired polynomial order.

  --> double xInput[n]:                    an array containing all the input x values

  --> double yInput[n]:                    an array containing all the input y values

  --> int polynomialOrder:                 the order of the desired output polynomial

  --> double outputArray[polynomialOrder]: a dummy array which is passed by reference so that we can
                                           return the calculated coefficients.

PART 3 - Implementation
-----------------------------------
This module implements polynomial regression by the method of least squares, which seeks
to minimise the variance of the estimates of the coefficients. This can be expressed
in terms of matrices by b = (((X^t)X)^(-1))(X^t)y where:
          [1 x1 (x1^2) ... (x1^m)
           1 x2 (x2^2) ... (x2^m)
 --> x =   . .     .         .      --> m is the degree of the polynomial,
           : :     :         :          n is the number of input values given
           1 xn (xn^2) ... (xn^m)]

 --> y = the corresponding input vector of y values

 --> b = the estimates for the coefficients for the corresponding power in the polynomial

To implement this, I implemented several methods which perform several matrix operations
on 2d-arrays:
  --> multiply matrices
  --> transpose a matrix
  --> calculate a determinant using Laplace's formula
  --> create a co-factor matrix
  --> find the inverse of a matrix using the co-factor matrix and determinant
However, as these methods are implementation specific, they are not exposed to the user.

PART 4 - Limitations
-----------------------------------
One of the major limitations of this module is it's memory usage. Since Laplace's formula
(which is used to calculate the determinants and inverses) is recursive, the memory
required will scale by O(n^n) (I think) because it will have to allocate space
for new matrices multiple times in every recursive iteration.

Several of the methods in this module also have runtime complexities of O(n^2) (multiplyMatrices,
matrixTranspose, and inverseMatrix), whilst the determinant function has a complexity
of O(n^3) and the co-factor matrix function has a complexity of O(n^4).

This means that, whilst theoretically the module is accurate and can be used for arbitrarily large
data set, in practice it is only useful for relatively small datasets.
