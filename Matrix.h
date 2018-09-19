/******************************************************
 * Dioghenes
 * Polytechnic of Turin
 * 2017
 * Matrix v0.2
 ******************************************************/

#ifndef MATRIX_H
#define MATRIX_H

#define  _MATRIX_ERR      -1
#define  _MATRIX_VOID     0
#define  _MATRIX_OK       +1
#define  _PREC_LOW        100
#define  _PREC_MEDIUM     500
#define  _PREC_HIGH       1000

const int _PRECISION  =  _PREC_LOW;

using namespace std;

#include <iostream>
#include <cstdlib>
#include <cmath>

class Matrix
{
    public:
        // Constructors
        Matrix();                       // Create an empty object
        Matrix(int rows, int columns);  // Create a zero matrix
        Matrix(const Matrix& toCopy);   // Create a copy matrix

        // Cleaners
        void clear();                   // Set all elements to zero
        void reset();                   // Create a VOID_MATRIX
        void err();                     // Create a ERR_MATRIX

        //Printer
        void print() const;             // Print the matrix

        //Getters
        int getStatus() const;          // Get the error flag value
        int* size() const;              // Get a vector of two integers in the format rowsxcolumns
        bool isSquare() const;			// Is the matrix a square one
        bool isDominant() const;
		bool isStrictlyDominant() const;

        //Assignment
        void assign(int rows, int columns, double* vector);     // Assign a matrix to the Matrix object

        //Operations
        void operator=(const Matrix& assignMatrix);             // Assign a matrix to a new one. All previous values are discarded.
        double& operator()(int row, int column) const;          // Get the address of an element
        Matrix operator+(const Matrix& addMatrix) const;        // Addiction between matrices
        Matrix operator-(const Matrix& subtractMatrix) const;   // Subtraction between matrices
        Matrix operator*(const Matrix& multMatrix) const;       // Multiplication between matrices


    private:
        int _rows;
        int _columns;
        double** _matrix;
        int _statusFlag;
};



Matrix eye(int size);
Matrix zeros(int rows, int columns);
Matrix ones(int rows, int columns);
Matrix transpose(const Matrix& matrix);
Matrix Vander(int size, double* vector, int vectorlength);
Matrix Jacobi(const Matrix& matrixA, double* vectorb);
Matrix GaussSeidel(const Matrix& matrixA, double* vectorb);
Matrix operator*(double scalar, const Matrix& multMatrix);
Matrix operator*(const Matrix& multMatrix, double scalar);
Matrix operator^(const Matrix& matrixA, int exponent);
Matrix operator/(const Matrix& matrixA, double scalar);
Matrix operator/(double scalar, const Matrix& matrixA);
bool isDominant(const Matrix& matrix);
bool isStrictlyDominant(const Matrix& matrix);
double det(const Matrix& matrix);

#endif



