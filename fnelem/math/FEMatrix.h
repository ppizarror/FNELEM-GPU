/**
FNELEM-GPU MATRIX UTILITARY FUNCTIONS
Performs matrix inversion using Gauss Jordan algorithm.
Based on: https://github.com/ZhengzhongSun/Matrix-Inversion-with-CUDA

@package fnelem.math
@author ppizarror
@date 20/11/2018
@license
    MIT License
    Copyright (c) 2018 Pablo Pizarro R.

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

// Library imports
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string>

// Constant definition
#define FEMATRIX_MIN_INVERSION_VALUE 0.0005

/**
 * Matrix class for working with CUDA. Stores matrix in an array [1..n*m].
 */
class FEMatrix {
private:

    // Number of rows
    int n = 0;

    // Number of columns
    int m = 0;

    // Matrix data
    double *mat;

    // Origin from one
    int origin = 0;

public:

    // Constructor
    FEMatrix(int n, int m);

    // Create matrix from array
    FEMatrix(int n, int m, double *matrix);

    // Destructor
    ~FEMatrix();

    // Set origin
    void set_origin(int o);

    // Fill matrix with zeros
    void fill_zeros();

    // Fill matrix with ones
    void fill_ones();

    // Display matrix in console
    void disp() const;

    // Update matrix A[i][j] = val
    void set(int i, int j, double val);

    // Returns value A[i][j]
    double get(int i, int j) const;

    // Save matrix to file
    void save_to_file(std::string filename) const;

    // Get matrix array
    double *get_array() const;

    // Return matrix dimension
    int *get_dimension() const;

    // Get square dimension
    int get_square_dimension() const;

    // Check if matrix is square nxn
    bool is_square() const;

    // Assign
    FEMatrix &operator=(const FEMatrix &matrix);

    // Adds a matrix with self
    FEMatrix &operator+=(const FEMatrix &matrix);

    // Add and return a new matrix
    FEMatrix operator+(const FEMatrix &matrix) const;

    // Substract a matrix with self
    FEMatrix &operator-=(const FEMatrix &matrix);

    // Substract and return a new matrix
    FEMatrix operator-(const FEMatrix &matrix) const;

    // Unary substract
    FEMatrix operator-() const;

    // Matrix transpose
    void transpose();

    // Matrix multiplication with self
    FEMatrix &operator*=(const FEMatrix &matrix);

    // Clone object
    FEMatrix clone() const;

};

/**
 * Creates matrix.
 *
 * @param n Number of rows
 * @param m Number of columns
 */
FEMatrix::FEMatrix(int n, int m) {
    this->n = n;
    this->m = m;
    this->mat = new double[n * m];
    this->fill_zeros();
}

/**
 * Generates matrix from array.
 *
 * @param matrix Array
 * @param n Number of rows
 * @param m Number of columns
 */
FEMatrix::FEMatrix(int n, int m, double *matrix) {
    this->n = n;
    this->m = m;
    this->mat = new double[n * m];
    for (int i = 0; i < n; i++) { // Rows
        for (int j = 0; j < m; j++) { // Columns
            this->mat[i * m + j] = matrix[i * m + j];
        }
    }
}

/**
 * Destroy matrix.
 */
FEMatrix::~FEMatrix() {
    delete[] this->mat;
}

/**
 * Fill matrix with zeros.
 */
void FEMatrix::fill_zeros() {
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            this->mat[i * this->m + j] = 0;
        }
    }
}

/**
 * Fill matrix with ones.
 */
void FEMatrix::fill_ones() {
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            this->mat[i * this->m + j] = 1;
        }
    }
}

/**
 * Display matrix in console.
 */
void FEMatrix::disp() const {

    // Get max number length
    

    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            std::cout << this->mat[i * this->m + j] << "\t";
        }
        std::cout << "" << std::endl;
    }
    std::cout << "" << std::endl;
}

/**
 * Updates matrix value
 * @param i Row position
 * @param j Column position
 * @param val Value
 */
void FEMatrix::set(int i, int j, double val) {
    if (i >= this->n + this->origin || j >= this->m + this->origin) {
        throw std::logic_error("[FEMATRIX] Column or row position overflow matrix");
    }
    this->mat[(i - this->origin) * this->m + (j - this->origin)] = val;
}

/**
 * Returns matrix value.
 *
 * @param i Row position
 * @param j Column position
 * @return Value at matrix[i][j]
 */
double FEMatrix::get(int i, int j) const {
    if (i >= this->n + this->origin || j >= this->m + this->origin) {
        throw std::logic_error("[FEMATRIX] Get value from matrix olumn or row position overflow");
    }
    return this->mat[(i - this->origin) * this->m + (j - this->origin)];
}


/**
 * Save matrix to file.
 *
 * @param filename Filename.
 */
void FEMatrix::save_to_file(std::string filename) const {
    std::ofstream plik;
    plik.open(filename);
    for (int j = 0; j < this->n; j++) {
        for (int i = 0; i < this->m; i++) {
            plik << this->mat[j * this->m + i] << "\t";
        }
        if (j < this->n - 1) {
            plik << std::endl;
        }
    }
    plik.close();
}

/**
 * Return matrix array.
 *
 * @return
 */
double *FEMatrix::get_array() const {
    double *mat_array = new double[this->n * this->m];
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            mat_array[i * this->m + j] = this->mat[i * this->m + j];
        }
    }
    return mat_array;
}

/**
 * Return matrix dimension.
 *
 * @return
 */
int *FEMatrix::get_dimension() const {
    int *dim = new int[2];
    dim[0] = this->n;
    dim[1] = this->m;
    return dim;
}

/**
 * Return true/false if matrix is square.
 * @return
 */
bool FEMatrix::is_square() const {
    return this->n == this->m;
}

/**
 * Returns square dimension.
 *
 * @return
 */
int FEMatrix::get_square_dimension() const {
    if (!is_square()) {
        return 0;
    } else {
        return this->n;
    }
}

/**
 * Asignation operation.
 *
 * @param matrix
 * @return
 */
FEMatrix &FEMatrix::operator=(const FEMatrix &matrix) {
    this->n = matrix.n;
    this->m = matrix.m;
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            this->mat[i * this->m + j] = matrix.get(i, j);
        }
    }
    return *this;
}

/**
 * Adds a matrix with self.
 *
 * @param matrix Matrix to add
 * @return
 */
FEMatrix &FEMatrix::operator+=(const FEMatrix &matrix) {

    // Get matrix dimension
    if (matrix.n != this->n || matrix.m != this->m) {
        throw std::logic_error("[FEMATRIX] Matrix dimension must be the same");
    }

    // Adds
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            this->mat[i * this->m + j] += matrix.get(i, j);
        }
    }

    // Returns self pointer
    return *this;

}

/**
 * Adds with a matrix and return new.
 *
 * @param matrix Matrix to add
 * @return
 */
FEMatrix FEMatrix::operator+(const FEMatrix &matrix) const {

    // Get matrix dimension
    if (matrix.n != this->n || matrix.m != this->m) {
        throw std::logic_error("[FEMATRIX] Matrix dimension must be the same");
    }

    FEMatrix newMatrix = FEMatrix(this->n, this->m);
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            newMatrix.set(i, j, this->get(i, j) + matrix.get(i, j));
        }
    }

    // Return matrix pointer
    return newMatrix;

}

/**
 * Substract a matrix with self.
 *
 * @param matrix Matrix to add
 * @return
 */
FEMatrix &FEMatrix::operator-=(const FEMatrix &matrix) {

    // Get matrix dimension
    if (matrix.n != this->n || matrix.m != this->m) {
        throw std::logic_error("[FEMATRIX] Matrix dimension must be the same");
    }

    // Adds
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            this->mat[i * this->m + j] -= matrix.get(i, j);
        }
    }

    // Returns self pointer
    return *this;

}

/**
 * Substract with a matrix and return new.
 *
 * @param matrix Matrix to add
 * @return
 */
FEMatrix FEMatrix::operator-(const FEMatrix &matrix) const {

    // Get matrix dimension
    if (matrix.n != this->n || matrix.m != this->m) {
        throw std::logic_error("[FEMATRIX] Matrix dimension must be the same");
    }

    FEMatrix newMatrix = FEMatrix(this->n, this->m);
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            newMatrix.set(i, j, this->get(i, j) - matrix.get(i, j));
        }
    }

    // Return matrix pointer
    return newMatrix;

}

/**
 * Unary substract.
 *
 * @return
 */
FEMatrix FEMatrix::operator-() const {
    FEMatrix newMatrix = FEMatrix(this->n, this->m);
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            newMatrix.set(i, j, -this->get(i, j));
        }
    }
    return newMatrix;
}

/**
 * Matrix transpose.
 */
void FEMatrix::transpose() {

    // Save temporal dimension
    int tempdim = this->n;

    // Create new transposed matrix
    double newMat[this->n * this->m];
    for (int i = 0; i < this->m; i++) { // Rows
        for (int j = 0; j < this->n; j++) { // Columns
            newMat[i * this->n + j] = this->get(j, i);
        }
    }

    // Update self matrix
    for (int i = 0; i < this->m; i++) { // Rows
        for (int j = 0; j < this->n; j++) { // Columns
            this->mat[i * this->n + j] = newMat[i * this->n + j];
        }
    }

    // Update matrix
    this->n = this->m;
    this->m = tempdim;

}

/**
 * Matrix multiplication with self.
 *
 * @param matrix Matrix to multiply
 * @return
 */
FEMatrix &FEMatrix::operator*=(const FEMatrix &matrix) {

    // Check dimension
    if (this->m != matrix.n) {
        throw std::logic_error("[FEMATRIX] Can't multiply matrix, dimension doest not agree");
    }

    // Create new auxiliar matrix AXB = (this) AXN * (matrix) NXB
    int a = this->n;
    int b = matrix.m;
    double auxMatrix[a * b];
    double sum = 0; // Stores partial sum

    // Multiply
    for (int i = 0; i < a; i++) { // Rows of new matrix
        for (int j = 0; j < b; j++) { // Columns of new matrix
            sum = 0;
            for (int k = 0; k < this->m; k++) {
                sum += this->get(i, k) * matrix.get(k, j);
            }
            auxMatrix[i * b + j] = sum;
        }
    }

    // Update matrix
    this->n = a;
    this->m = b;
    for (int i = 0; i < a; i++) { // Rows of new matrix
        for (int j = 0; j < b; j++) { // Columns of new matrix
            this->mat[i * b + j] = auxMatrix[i * b + j];
        }
    }

    // Return self
    return *this;

}

/**
 * Clones matrix.
 *
 * @return
 */
FEMatrix FEMatrix::clone() const {
    FEMatrix newMatrix = FEMatrix(this->n, this->m);
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            newMatrix.set(i, j, this->get(i, j));
        }
    }
    return newMatrix;
}

/**
 * Set matrix origin.
 *
 * @param o New origin
 */
void FEMatrix::set_origin(int o) {
    this->origin = o;
}