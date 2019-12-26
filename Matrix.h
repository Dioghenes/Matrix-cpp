/******************************************************
 * Dioghenes
 * Polytechnic of Turin
 * 2019
 * Matrix v0.3
 ******************************************************/

#ifndef MATRIX_H
#define MATRIX_H

#define  _MATRIX_ERR      		-1
#define  _MATRIX_VOID    		0
#define  _MATRIX_OK   			+1
#define  _MATRIX_PREC_LOW       100
#define  _MATRIX_PREC_MEDIUM    500
#define  _MATRIX_PREC_HIGH      1000

#define	 _MATRIX_PRECISION	  	_MATRIX_PREC_LOW

using namespace std;

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>

static double _errmatrix = 0.0;

class Matrix
{
    public:
        // Constructors
        Matrix();                       // Create an empty object
        Matrix(int rows, int columns);  // Create a zero matrix
        Matrix(const Matrix& toCopy);   // Create a copy matrix
		Matrix(string name);						// Create a named empty object
		Matrix(int rows, int columns, string name);	// Create a named zero matrix
		Matrix(const Matrix& toCopy, string name);	// Create a named copy matrix

        // Cleaners
        void clear();                   // Set all elements to zero
        void reset();                   // Set matrix to VOID
		void _err();					// Set matrix to ERR

		//Printer
        void print() const;             // Print the matrix

        //Getters
        int getStatus() const;          // Get the error flag value
        string getName() const;			// Get the matrix name
		int* size() const;              // Get a vector of two integers in the format rowsxcolumns
        bool isSquare() const;			// Is the matrix a square one
        bool isDominant() const;		// Is the matrix dominant by row
		bool isStrictlyDominant() const;// Is the matrix strictly dominant by row

        //Assignment
        void setName(string name);								// Assign a name to the matrix
		void assign(int rows, int columns, double* vector);     // Assign a matrix to the Matrix object

        //Operations
        void operator=(const Matrix& assignMatrix);             // Assign a matrix to a new one. All previous values are discarded.
        double& operator()(int row, int column) const;          // Get the address of an element
        Matrix operator+(const Matrix& addMatrix) const;        // Addiction between matrices
        Matrix operator-(const Matrix& subtractMatrix) const;   // Subtraction between matrices
        Matrix operator*(const Matrix& multMatrix) const;       // Multiplication between matrices


    private:
		string _name;
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
Matrix operator^(const Matrix& expMatrix, int exponent);
Matrix operator/(const Matrix& divMatrix, double scalar);
Matrix operator/(double scalar, const Matrix& divMatrix);
bool isDominant(const Matrix& matrix);
bool isStrictlyDominant(const Matrix& matrix);
double det(const Matrix& matrix);

#endif



