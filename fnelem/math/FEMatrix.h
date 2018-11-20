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

public:

    // Constructor
    FEMatrix(int n, int m);

    // Create matrix from array
    FEMatrix(double *matrix, int n, int m);

    // Destructor
    ~FEMatrix();

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
FEMatrix::FEMatrix(double *matrix, int n, int m) {
    this->n = n;
    this->m = m;
    this->mat = new double[n * m];
    for (int i = 0; i < n; i++) { // Rows
        for (int j = 0; j < m; j++) { // Columns
            this->mat[i * n + j] = matrix[i * n + j];
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
            this->mat[i * this->n + j] = 0;
        }
    }
}

/**
 * Fill matrix with ones.
 */
void FEMatrix::fill_ones() {
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            this->mat[i * this->n + j] = 1;
        }
    }
}

/**
 * Display matrix in console.
 */
void FEMatrix::disp() const {
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            std::cout << this->mat[i * this->n + j] << "\t";
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
    if (i >= this->n || j >= this->m) {
        throw std::logic_error("Column or row position overflow matrix");
    }
    this->mat[i * this->n + j] = val;
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
            plik << this->mat[j * this->n + i] << "\t";
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
            mat_array[i * this->n + j] = this->mat[i * this->n + j];
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
 * Returns matrix value.
 *
 * @param i Row position
 * @param j Column position
 * @return Value at matrix[i][j]
 */
double FEMatrix::get(int i, int j) const {
    if (i >= this->n || j >= this->m) {
        throw std::logic_error("Column or row position overflow matrix");
    }
    return this->mat[i * this->n + j];
}
