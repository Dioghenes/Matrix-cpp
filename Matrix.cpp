/******************************************************
 * Dioghenes
 * Polytechnic of Turin
 * 2019
 * Matrix v0.3
 ******************************************************/

#include "Matrix.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;


/****************************************************
 * Inner functions
 ****************************************************/
double _det(const Matrix& matrix, int dim);


/****************************************************
 * Constructors
 ****************************************************/

// Allocate space for an empty matrix
Matrix::Matrix():
    _rows(0),
    _columns(0),
    _matrix(NULL),
    _statusFlag(_MATRIX_VOID)
{}

// Allocate space for a generic matrix initialized to zero
Matrix::Matrix(int rows, int columns):
    _rows(rows),
    _columns(columns),
	_matrix(NULL),
    _statusFlag(_MATRIX_OK)
{
    _matrix = (double**)malloc(_rows*sizeof(double*));
    if(_matrix == NULL){
        this->_err();
        return;
    }
    for(int i=0;i<_rows;i++){
        _matrix[i] = (double*)malloc(_columns*sizeof(double));
        if(_matrix[i] == NULL){
            this->_err();
            return;
        }
    }
    this->clear();
	return;
}

// Copy constructor for a generic matrix
Matrix::Matrix(const Matrix& toCopy):
    _rows(toCopy._rows),
    _columns(toCopy._columns),
    _statusFlag(toCopy._statusFlag)
{
    _matrix = (double**)malloc(_rows*sizeof(double*));
    if(_matrix == NULL){
		this->_err();
        return;
    }
    for(int i=0;i<_rows;i++){
        _matrix[i] = (double*)malloc(_columns*sizeof(double));
        if(_matrix[i] == NULL){
            this->_err();
            return;
        }
    }
    for(int i=0;i<_rows;i++){
        for(int j=0;j<_columns;j++){
            _matrix[i][j] = toCopy._matrix[i][j];
        }
    }
	return;
}


/****************************************************
 * Named Constructors
 ****************************************************/

// Allocate space for an empty matrix
Matrix::Matrix(string name):
	_name(name),
    _rows(0),
    _columns(0),
    _matrix(NULL),
    _statusFlag(_MATRIX_VOID)
{}

// Allocate space for a generic matrix initialized to zero
Matrix::Matrix(int rows, int columns, string name):
    _name(name),
	_rows(rows),
    _columns(columns),
    _matrix(NULL),
    _statusFlag(_MATRIX_OK)
{
    _matrix = (double**)malloc(_rows*sizeof(double*));
    if(_matrix == NULL){
        this->_err();
        return;
    }
    for(int i=0;i<_rows;i++){
        _matrix[i] = (double*)malloc(_columns*sizeof(double));
        if(_matrix[i] == NULL){
            this->_err();
            return;
        }
    }
    this->clear();
    return;
}

// Copy constructor for a generic matrix
Matrix::Matrix(const Matrix& toCopy, string name):
	_name(name),
    _rows(toCopy._rows),
    _columns(toCopy._columns),
    _statusFlag(toCopy._statusFlag)
{
    _matrix = (double**)malloc(_rows*sizeof(double*));
    if(_matrix == NULL){
        this->_err();
        return;
    }
    for(int i=0;i<_rows;i++){
        _matrix[i] = (double*)malloc(_columns*sizeof(double));
        if(_matrix[i] == NULL){
            this->_err();
            return;
        }
    }
    for(int i=0;i<_rows;i++){
        for(int j=0;j<_columns;j++){
            _matrix[i][j] = toCopy._matrix[i][j];
        }
    }
    return;
}


/****************************************************
 * Cleaners
 ****************************************************/

// Set all the elements of the matrix to 0
void Matrix::clear(){
	for(int i=0;i<_rows;i++){
        for(int j=0;j<_columns;j++){
            _matrix[i][j] = 0;
        }
    }
	return;
}

// Reset the matrix, deallocating memory and resetting flags
void Matrix::reset(){
    for(int i=0;i<_rows;i++){
		if(_matrix[i]!=NULL) free(_matrix[i]);
    }
    if(_matrix!=NULL) free(_matrix);
    _statusFlag = _MATRIX_VOID;
    _rows = 0;
    _columns = 0;
	_matrix = NULL;
	_name = "";
	return;
}

// Create an _MATRIX_ERR matrix with flag set to _MATRIX_ERR
void Matrix::_err(){
    for(int i=0;i<_rows;i++){
        if(_matrix[i]!=NULL) free(_matrix[i]);
    }
    if(_matrix!=NULL) free(_matrix);
    _statusFlag = _MATRIX_ERR;
    _rows = 0;
    _columns = 0;
	_matrix = NULL;
	_name = "";
	return;
}


/****************************************************
 * Printer
 ****************************************************/
 
// Print the matrix on the screen
void Matrix::print() const{
    if(_statusFlag == _MATRIX_ERR){
        cout << "Matrix " << _name << " is invalid." << endl;
        return;
    }
    if(_statusFlag == _MATRIX_VOID){
        cout << "Matrix " << _name << "is void." << endl;
        return;
    }
    cout << endl << "Matrix " << _name << endl;
    for(int i=0;i<_rows;i++){
        for(int j=0;j<_columns;j++){
            cout << _matrix[i][j] << "\t\t";
        }
        cout << endl;
    }
    cout << endl;
    return;
}


/****************************************************
 * Getters
 ****************************************************/
 
// Get the matrix flag
int Matrix::getStatus() const{
    return _statusFlag;
}

// Get the matrix name
string Matrix::getName() const{
	return _name;
}

// Get the size: [0] is the rows size, [1] is the columns one
int* Matrix::size() const{
	int* size;
    size = (int*)malloc(2*sizeof(int));
    size[0] = _rows;
    size[1] = _columns;
    return size;
}

// Is the matrix a square one?
bool Matrix::isSquare() const{
	return _rows==_columns;
}

// Is the matrix dominant-by-row?
bool Matrix::isDominant() const{
	if(!this->isSquare()) return 0;
    bool isDom = 1;
    double tmp = 0;
    for(int i=0;i<_rows;i++){
        for(int j=0;j<_columns;j++){
            if(i==j) continue;
            tmp = tmp + abs(_matrix[i][j]);
        }
        if(abs(_matrix[i][i]) < tmp){
            isDom = 0;
            break;
        }
        tmp = 0;
    }
    return isDom;
}

// Is the matrix strictly dominant-by-row?
bool Matrix::isStrictlyDominant() const{
	if(!this->isSquare()) return 0;
    bool isDom = 1;
    double tmp = 0;
    for(int i=0;i<_rows;i++){
        for(int j=0;j<_columns;j++){
            if(i==j) continue;
            tmp = tmp + abs(_matrix[i][j]);
        }
        if(abs(_matrix[i][i]) <= tmp){
            isDom = 0;
            break;
        }
        tmp = 0;
    }
    return isDom;
}


/****************************************************
 * Assignment
 *****************************************************/

// Assign a name to the matrix
void Matrix::setName(string name){
	_name = name;
	return;
} 

// Reset and set to new values
void Matrix::assign(int rows, int columns, double* vector){
    this->reset();
    _rows = rows;
    _columns = columns;

    _matrix = (double**)malloc(_rows*sizeof(double*));
    if(_matrix == NULL){
        this->_err();
        return;
    }
    for(int i=0;i<_rows;i++){
        _matrix[i] = (double*)malloc(_columns*sizeof(double));
        if(_matrix[i] == NULL){
			this->_err();
			return;
		}
    }

    int k = 0;
    for(int i=0;i<_rows;i++){
        for(int j=0;j<_columns;j++){
            _matrix[i][j] = vector[k];
            k+=1;
        }
    }
    _statusFlag = _MATRIX_OK;
    return;
}


/****************************************************
 * Operations
 ****************************************************/

// Direct read of an element of the matrix in the form myMatrix(r,c)
double& Matrix::operator()(int row, int column) const{
	if(_statusFlag == _MATRIX_OK && (row<_rows && column<_columns)) return _matrix[row][column];
	else return _errmatrix;
}

// Copy a matrix into another one
void Matrix::operator=(const Matrix& assignMatrix){
	this->reset();

	// ERR dest matrix in case of ERR source matrix
	if(assignMatrix.getStatus() != _MATRIX_OK){
        this->_err();
        return;
    }

    _rows = assignMatrix._rows;
    _columns = assignMatrix._columns;
    _statusFlag = assignMatrix._statusFlag;

    _matrix = (double**)malloc(_rows*sizeof(double*));
    if(_matrix == NULL){
		this->_err();
		return;
	}
    for(int i=0;i<_rows;i++){
        _matrix[i] = (double*)malloc(_columns*sizeof(double));
        if(_matrix[i] == NULL){
			this->_err();
			return;
		}
    }

    for(int i=0;i<_rows;i++){
        for(int j=0;j<_columns;j++){
            _matrix[i][j] = assignMatrix._matrix[i][j];
		}
    }
}

// Add two matrix and create a new one with the result
Matrix Matrix::operator+(const Matrix& addMatrix) const{
    Matrix retMatrix(addMatrix._rows,addMatrix._columns);

	// ERR dest matrix in case of ERR source matrix
	if(addMatrix.getStatus() != _MATRIX_OK){
        retMatrix._err();
        return retMatrix;
    }

	// Check the correctness of dimensions
    if(_rows != addMatrix._rows || _columns != addMatrix._columns){
        retMatrix._err();
        return retMatrix;
    }

    for(int i=0;i<_rows;i++){
		for(int j=0;j<_columns;j++){
			retMatrix(i,j) = _matrix[i][j]+addMatrix._matrix[i][j];
		}
    }

    return retMatrix;
}

// Subtract two matrix and create a new one with the result
Matrix Matrix::operator-(const Matrix& subtractMatrix) const{
    Matrix retMatrix(subtractMatrix._rows,subtractMatrix._columns);

	// ERR dest matrix in case of ERR source matrix
	if(subtractMatrix.getStatus() != _MATRIX_OK){
        retMatrix._err();
        return retMatrix;
    }

	// Check the correctness of dimensions
    if(_rows != subtractMatrix._rows || _columns != subtractMatrix._columns){
        retMatrix._err();
        return retMatrix;
    }

    for(int i=0;i<_rows;i++){
		for(int j=0;j<_columns;j++){
			retMatrix(i,j)= _matrix[i][j]-subtractMatrix._matrix[i][j];
		}
    }

    return retMatrix;
}

// Multiply two matrix and create a new one with the result
Matrix Matrix::operator*(const Matrix& multMatrix) const{
    Matrix retMatrix(_rows,multMatrix._columns);

	// ERR dest matrix in case of ERR source matrix
	if(multMatrix.getStatus() != _MATRIX_OK){
        retMatrix._err();
        return retMatrix;
    }

	 // Check the correctness of dimensions
    if(_columns != multMatrix._rows){
        retMatrix._err();
        return retMatrix;
    }

    double tmp = 0;
    for(int i=0;i<_rows;i++){
        for(int k=0;k<multMatrix._columns;k++){
            for(int j=0;j<_columns;j++){
                tmp = tmp + _matrix[i][j]*multMatrix._matrix[j][k];
            }
            retMatrix(i,k) = tmp;
            tmp = 0;
        }
    }
    return retMatrix;
}


/****************************************************
 * Functions
 ****************************************************/

// Create an eye matrix of size x size
Matrix eye(int size){
    Matrix retMatrix(size,size);
    for(int i=0;i<size;i++){
        retMatrix(i,i) = 1;
    }
    return retMatrix;
}

// Create a zero matrix of rows x columns
Matrix zeros(int rows, int columns){
    Matrix retMatrix(rows,columns);
    for(int i=0;i<rows;i++){
        for(int j=0;j<columns;j++){
            retMatrix(i,j) = 0;
        }
    }
    return retMatrix;
}

// Create a matrix of ones of rows x columns
Matrix ones(int rows, int columns){
    Matrix retMatrix(rows,columns);
    for(int i=0;i<rows;i++){
        for(int j=0;j<columns;j++){
            retMatrix(i,j) = 1;
        }
    }
    return retMatrix;
}

// Returns the transposed matrix of a given matrix
Matrix transpose(const Matrix& matrix){
    int *dimension;
    dimension = (int*)malloc(2*sizeof(int));
    dimension = matrix.size();
    Matrix retMatrix(dimension[1],dimension[0]);

	// ERR dest matrix in case of ERR source matrix
	if(matrix.getStatus() != _MATRIX_OK){
        retMatrix._err();
        return retMatrix;
    }

    for(int i=0;i<dimension[1];i++){
        for(int j=0;j<dimension[0];j++){
            retMatrix(i,j) = matrix(j,i);
        }
    }
    return retMatrix;
}

// Create a Vandermonde matrix from a vector
Matrix Vander(int size, double* vector, int vectorlength){
    Matrix retMatrix(vectorlength,size);
    for(int i=0;i<vectorlength;i++){
        for(int j=0;j<size;j++){
            retMatrix(i,j) = pow(vector[i],double(j));
        }
    }
    return retMatrix;
}

// Solve a linear system using Jacobi method
Matrix Jacobi(const Matrix& matrixA, double* vectorb){
    int* dimensions = (int*)malloc(2*sizeof(int));
    dimensions = matrixA.size();
    Matrix solutionVector(dimensions[0],1);

	// ERR dest matrix in case of ERR source matrix
	if(matrixA.getStatus() != _MATRIX_OK){
        solutionVector._err();
        return solutionVector;
    }

    if(isStrictlyDominant(matrixA) == 0){
        solutionVector._err();
        return solutionVector;
    }
    Matrix copyVector = solutionVector;

    int kMax = _MATRIX_PRECISION;
    double tmp;
    for(int k=0;k<kMax;k++){
        for(int i=0;i<dimensions[0];i++){
            tmp = 0;
            for(int j=0;j<dimensions[1];j++){
                if(i==j) continue;
                tmp = tmp + copyVector(j,0)*matrixA(i,j);
            }
            solutionVector(i,0) = (vectorb[i] - tmp)/matrixA(i,i);
        }
        copyVector = solutionVector;
    }
    return solutionVector;
}

// Solve a linear system using Gauss-Seidel method
Matrix GaussSeidel(const Matrix& matrixA, double* vectorb){
    int* dimensions = (int*)malloc(2*sizeof(int));
    dimensions = matrixA.size();
    Matrix solutionVector(dimensions[0],1);

	// ERR dest matrix in case of ERR source matrix
	if(matrixA.getStatus() != _MATRIX_OK){
        solutionVector._err();
        return solutionVector;
    }

    if(isStrictlyDominant(matrixA) == 0){
        solutionVector._err();
        return solutionVector;
    }

    int kMax = _MATRIX_PRECISION;
    double tmp;
    for(int k=0;k<kMax;k++){
        for(int i=0;i<dimensions[0];i++){
            tmp = 0;
            for(int j=0;j<dimensions[1];j++){
                if(i==j) continue;
                tmp = tmp + solutionVector(j,0)*matrixA(i,j);
            }
            solutionVector(i,0) = (vectorb[i] - tmp)/matrixA(i,i);
        }
    }
    return solutionVector;
}

//Multiplication scalar * matrix
Matrix operator*(double scalar, const Matrix& multMatrix){
    Matrix retMatrix(multMatrix);

	// ERR dest matrix in case of ERR source matrix
    if(multMatrix.getStatus() != _MATRIX_OK){
        retMatrix._err();
        return retMatrix;
    }

	int* dimension = (int*)malloc(2*sizeof(int));
    dimension = retMatrix.size();
    for(int i=0;i<dimension[0];i++){
        for(int j=0;j<dimension[1];j++){
            retMatrix(i,j) = retMatrix(i,j)*scalar;
        }
    }
    return retMatrix;
}

// Multiplication matrix * scalar
Matrix operator*(const Matrix& multMatrix, double scalar){
    Matrix retMatrix(multMatrix);

	// ERR dest matrix in case of ERR source matrix
	if(multMatrix.getStatus() != _MATRIX_OK){
        retMatrix._err();
        return retMatrix;
    }

	int* dimension = (int*)malloc(2*sizeof(int));
    dimension = retMatrix.size();
    for(int i=0;i<dimension[0];i++){
        for(int j=0;j<dimension[1];j++){
            retMatrix(i,j) = retMatrix(i,j)*scalar;
        }
    }
    return retMatrix;
}

// Power of a matrix ^ scalar
Matrix operator^(const Matrix& expMatrix, int exponent){
    Matrix retMatrix(expMatrix);

	// ERR dest matrix in case of ERR source matrix
    if(expMatrix.getStatus() != _MATRIX_OK){
        retMatrix._err();
        return retMatrix;
    }

	if(!expMatrix.isSquare()){
		retMatrix._err();
		return retMatrix;
	}
    for(int i=0;i<exponent-1;i++){
        retMatrix = retMatrix*expMatrix;
    }
    return retMatrix;
}

// Divide a each matrix element by a scalar
Matrix operator/(const Matrix& divMatrix, double scalar){
    Matrix retMatrix(divMatrix);

	// ERR dest matrix in case of ERR source matrix
	if(divMatrix.getStatus() != _MATRIX_OK){
        retMatrix._err();
        return retMatrix;
    }

    int* dimensions = (int*)malloc(2*sizeof(int));
    dimensions = divMatrix.size();
    for(int i=0;i<dimensions[0];i++){
        for(int j=0;j<dimensions[1];j++){
            retMatrix(i,j) = divMatrix(i,j)/scalar;
        }
    }
    return retMatrix;
}

// Divide a scalar by each matrix element 
Matrix operator/(double scalar, const Matrix& divMatrix){
    Matrix retMatrix(divMatrix);

	// ERR dest matrix in case of ERR source matrix
    if(divMatrix.getStatus() != _MATRIX_OK){
		retMatrix._err();
		return retMatrix;
	}

	int* dimensions = (int*)malloc(2*sizeof(int));
    dimensions = divMatrix.size();
    for(int i=0;i<dimensions[0];i++){
        for(int j=0;j<dimensions[1];j++){
            retMatrix(i,j) = scalar/divMatrix(i,j);
        }
    }
    return retMatrix;
}

// Is the matrix dominant-by-row?
bool isDominant(const Matrix& matrix){
	// ERR source matrix return value
	if(matrix.getStatus() != _MATRIX_OK) return 0;

	if(!matrix.isSquare()) return 0;
    bool isDom = 1;
    double tmp = 0;
    int* dimensions = (int*)malloc(2*sizeof(int));
    dimensions = matrix.size();
    for(int i=0;i<dimensions[0];i++){
        for(int j=0;j<dimensions[1];j++){
            if(i==j) continue;
            tmp = tmp + abs((float)matrix(i,j));
        }
        if(abs(matrix(i,i)) < tmp){
            isDom = 0;
            break;
        }
        tmp = 0;
    }
    return isDom;
}

// Is the matrix strictly dominant-by-row?
bool isStrictlyDominant(const Matrix& matrix){
	// ERR source matrix return value
	if(matrix.getStatus() != _MATRIX_OK) return 0;

	if(!matrix.isSquare()) return 0;
    bool isDom = 1;
    double tmp = 0;
    int* dimensions = (int*)malloc(2*sizeof(int));
    dimensions = matrix.size();
    for(int i=0;i<dimensions[0];i++){
        for(int j=0;j<dimensions[1];j++){
            if(i==j) continue;
            tmp = tmp + abs((double)matrix(i,j));
        }
        if(abs(matrix(i,i)) <= tmp){
            isDom = 0;
            break;
        }
        tmp = 0;
    }
    return isDom;
}

// Calculate the determinant of the matrix
double det(const Matrix& matrix){
	// ERR source matrix return value
    if(matrix.getStatus() != _MATRIX_OK) return 0;

	if(!matrix.isSquare()) return 0;
	int* dimensions = (int*)malloc(2*sizeof(int));
    dimensions = matrix.size();
    return _det(matrix,dimensions[0]);
}
double _det(const Matrix& matrix, int dim){
	double determ = 0;
	Matrix temp(dim,dim);
    if(dim==1) {
        return matrix(0,0);
    }
    else if(dim==2) {
        determ=(matrix(0,0)*matrix(1,1)-matrix(0,1)*matrix(1,0));
        return determ;
    }
    else{
        for(int p=0;p<dim;p++){
			int h = 0;
			int k = 0;
			for(int i=1;i<dim;i++){
				for(int j=0;j<dim;j++){
					if(j==p){
						continue;
					}
					temp(h,k) = matrix(i,j);
					k++;
					if(k==dim-1){
						h++;
						k = 0;
					}
				}
			}
			determ=determ+matrix(0,p)*pow(-1,(double)p)*_det(temp,dim-1);
        }
		return determ;
	}
}



