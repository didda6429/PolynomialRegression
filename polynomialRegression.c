#include "polynomialRegression.h"
#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "math.h"

void multiplyMatrices(int n, int y, int z, double firstInput[n][y], double secondInput[y][z], double output[n][z]){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < z; j++){
      double dotProduct = 0;
      for(int k = 0; k < y; k++){
        dotProduct += firstInput[i][k] * secondInput[k][j];
      }
      output[i][j] = dotProduct;
    }
  }
}

void matrixTranspose(int n, int k, double input[n][k], double output[k][n]){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < k; j++){
      output[j][i] = input[i][j];
    }
  }
}

double determinant(int n, double input[n][n]){
  if(n == 2){
    return ((input[0][0] * input[1][1]) - (input[0][1] * input[1][0]));
  } else {
    double returnVal = 0;
    for(int i = 0; i < n; i++){
      //extract the minors and recursively calculate their determinants
      double iMinor[n-1][n-1];
      int p = 0;
      int q = 0;
      for(int k = 0; k < n; k++){
        for(int j = 0; j < n; j++){
          if(k!= 0 && j != i){
            iMinor[p][q] = input[k][j];
            if(q < (j-2)){
              q++;
            } else {
              p = 0;
              q++;
            }
          }
        }
      }
      returnVal += pow(-1, i) * input[0][i] * determinant(n-1, iMinor);
    }
    return returnVal;
  }
}

void coFactorMatrix(int n, double input[n][n], double output[n][n]){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      //extract the minors and calculate their determinants
      double ijMinor[n-1][n-1];
      int p = 0;
      int q = 0;
      for(int k = 0; k < n; k++){
        for(int l = 0; l < n; l++){
          if(k!= i && l != j){
            ijMinor[p][q] = input[k][l];
            if(q < (l-2)){
              q++;
            } else {
              p = 0;
              q++;
            }
          }
        }
      }
      output[i][j] = pow(-1, i+j) * determinant(n-1, ijMinor);
    }
  }
}

void inverseMatrix(int n, double input[n][n], double output[n][n]){
  double cofactors[3][3];
  coFactorMatrix(n, input, cofactors);
  matrixTranspose(n, n, cofactors, output);
  double inputDeterminant = determinant(n, input);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      output[i][j] = output[i][j] / inputDeterminant;
    }
  }
}

void calculateCoefficients(int n, double xInput[n], double yInput[n], int polynomialOrder, double outputArray[polynomialOrder]){
  assert(n >= polynomialOrder); //to ensure we have enough data to find all coeffiecients

  //create the matrix containing the powers of xInput as a 2-d array
  double xMatrixT[n][polynomialOrder];
  for(int i = 0; i < polynomialOrder; i++){
    for(int j = 0; j < n; j++){
      xMatrixT[i][j] = pow(xInput[j], i);
    }
  }

  //convert the yInput into a 2d matrix we can work with
  double yMatrix[n][1];
  for(int i = 0; i < n; i++){
    yMatrix[i][0] = yInput[i];
  }

  //calculate the coeffiecients by least squares
  double xMatrix[polynomialOrder][n];
  matrixTranspose(n, polynomialOrder, xMatrixT, xMatrix);

  double xMatrixProduct[polynomialOrder][polynomialOrder];
  multiplyMatrices(n, polynomialOrder, n, xMatrixT, xMatrix, xMatrixProduct);

  double xMatrixProductInverse[polynomialOrder][polynomialOrder];
  inverseMatrix(polynomialOrder, xMatrixProduct, xMatrixProductInverse);

  double finalMatrix[polynomialOrder][n];
  multiplyMatrices(polynomialOrder, n, n, xMatrixProductInverse, xMatrixT, finalMatrix);

  double finalCoefficients[polynomialOrder][1];
  multiplyMatrices(polynomialOrder, n, 1, finalMatrix, yMatrix, finalCoefficients);

  for(int i = 0; i < polynomialOrder; i++){
    outputArray[i] = finalCoefficients[i][0];
  }
}

void testMultiplyMatrices(){
  double firstInput[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
  double secondInput[3][3] = {{10,11,12},{13,14,15},{16,17,18}};
  double output[3][3];
  multiplyMatrices(3,3,3,firstInput,secondInput,output);
  double testData[3][3] = {{84,90,96},{201,216,231},{318,342,366}};
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      assert(output[i][j] == testData[i][j]);
    }
  }
  printf("testMultiplyMatrices passed.\n");
}

void testMatrixTranspose(){
  double input[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
  double output[3][3];
  matrixTranspose(3, 3, input, output);
  double testData[3][3] = {{1,4,7},{2,5,8},{3,6,9}};
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      assert(output[i][j] == testData[i][j]);
    }
  }
  printf("testMatrixTranspose passed.\n");
}

void testDeterminant(){
  double input[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  assert(determinant(3, input) == 1);
  double firstInput[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
  assert(determinant(3, firstInput) == 0);
  double secondInput[4][4] = {{1,7,5,8},{5,1,4,6},{4,9,6,8},{4,9,6,2}};
  assert((int)determinant(4, secondInput) == -462);
  printf("testDeterminant passed.\n");
}

void testCoFactorMatrix(){
  double input[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
  double output[3][3];
  coFactorMatrix(3, input, output);
  double testData[3][3] = {{-3,6,-3},{6,-12,6},{-3,6,-3}};
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      assert(output[i][j] == testData[i][j]);
    }
  }
  printf("testCoFactorMatrix passed.\n");
}

void testInverseMatrix(){
  //We'll test this case by multiply the matrix with it's inverse (because we've already
  //tested multiplication works) and checking that the result is the identity matrix
  double input[3][3] = {{1,7,5},{5,1,4},{4,9,6}};
  double inverseInput[3][3];
  double output[3][3];
  inverseMatrix(3, input, inverseInput);
  multiplyMatrices(3,3,3,input, inverseInput, output);
  double testData[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      assert((output[i][j] - testData[i][j]) < 0.1);
    }
  }
  printf("testInverseMatrix passed.\n");
}

void testCalculateCoefficients(){
  //y=2x ==> coefficients are 0,1,0
  double xInput1[3] = {1,2,3};
  double yInput1[3] = {2,4,6};
  double output1[3];
  calculateCoefficients(3, xInput1, yInput1, 3, output1);
  assert(output1[0] == 0);
  assert(output1[1] == 2);
  assert(output1[2] == 0);

  //y=(x^2)+1 ==> coefficients are 1,0,1
  double xInput2[3] = {1,2,3};
  double yInput2[3] = {2,5,10};
  double output2[3];
  calculateCoefficients(3, xInput2, yInput2, 3, output2);
  assert(output2[0] == 1);
  assert(output2[1] == 0);
  assert(output2[2] == 1);
  printf("testCalculateCoefficients passed.\n");
}

void test(){
  testMultiplyMatrices();
  testMatrixTranspose();
  testDeterminant();
  testCoFactorMatrix();
  testInverseMatrix();
  testCalculateCoefficients();
  printf("all tests passed.\n");
}

int testMain(int n, char *args[n]){
  test();
  return 0;
}
