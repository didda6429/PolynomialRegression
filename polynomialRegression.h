/* This module implements polynomial regression for polynomials of arbitrary order.

   Note that this module only exposes a single function with the aim of improving
   it's modularity.

   Internally, this is done with a matrix-based backend, which was achieved by
   implementing functions (namely determinant, inverse, multiplication, transpose,
   and co-factor matrices) which perform matrix arithmetic on 2d-arrays (which
   are used to represent arrays).
*/


void calculateCoefficients(int n, double xInput[n], double yInput[n], int polynomialOrder, double outputArray[polynomialOrder]);
