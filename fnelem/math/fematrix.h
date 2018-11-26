/**
FNELEM-GPU MATRIX DEFINITION
FEMatrix implement a full matrix manipulation environment.

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

// Init header file
#ifndef FNELEM_GPU_MATH_FEMATRIX_H
#define FNELEM_GPU_MATH_FEMATRIX_H

// Constant definition
#define __FEMATRIX_ZERO_TOL 1e-12

// Library imports
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string>

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

    // Saves origin
    int origin_temp = 0;

    // Update matrix A[i][j] = val, no origin
    void _set(int i, int j, double val);

    // Returns value A[i][j], no origin
    double _get(int i, int j) const;

    // Uses pad or not to display matrix on console
    bool apply_pad = false;

    // Calcualtes determinant recursive
    double _det_recursive(double *matrix, int d) const;

    // Display a matrix in console
    void disp_matrix(double *matrix, int dim_n, int dim_m) const;

    // Check if matrix has been deleted
    bool deleted = false;

    // Matrix name
    std::string mat_name = "";

public:

    // Default constructor
    FEMatrix();

    // Constructor
    FEMatrix(int n, int m);

    // Create matrix from array
    FEMatrix(int n, int m, double *matrix);

    // Destructor
    ~FEMatrix();

    // Set origin
    void set_origin(int o);

    // Disables origin
    void disable_origin();

    // Enable origin
    void enable_origin();

    // Fill matrix with value
    void fill(double value);

    // Fill matrix with zeros
    void fill_zeros();

    // Fill matrix with ones
    void fill_ones();

    // Display matrix in console
    void disp() const;

    // Update matrix A[i][j] = val
    void set(int i, int j, double val);

    // Set value for vector A[i] = val
    void set(int i, double val);

    // Returns value A[i][j]
    double get(int i, int j) const;

    // Returns value for vector A[i]
    double get(int i) const;

    // Get row
    FEMatrix *get_row(int i, int from, int to) const;

    // Get full row
    FEMatrix *get_row(int i) const;

    // Get column
    FEMatrix *get_column(int j, int from, int to) const;

    // Get full column
    FEMatrix *get_column(int j) const;

    // Save matrix to file
    void save_to_file(std::string filename) const;

    // Get matrix array
    double *get_array() const;

    // Return matrix dimension [N, M]
    int *size() const;

    // Return max dimension, usefull for vector
    int length() const;

    // Get square dimension
    int get_square_dimension() const;

    // Assign
    FEMatrix &operator=(const FEMatrix *matrix);

    // Adds a matrix with self
    FEMatrix &operator+=(const FEMatrix &matrix);

    // Adds a matrix with self
    FEMatrix &operator+=(const FEMatrix *matrix);

    // Add and return a new matrix
    FEMatrix *operator+(const FEMatrix &matrix) const;

    // Substract a matrix with self
    FEMatrix &operator-=(const FEMatrix &matrix);

    // Substract a matrix with self
    FEMatrix &operator-=(const FEMatrix *matrix);

    // Substract and return a new matrix
    FEMatrix *operator-(const FEMatrix &matrix) const;

    // Unary substract
    FEMatrix *operator-() const;

    // Matrix transpose
    void transpose_self();

    // Matrix transpose and return new
    FEMatrix *transpose() const;

    // Matrix multiplication with self
    FEMatrix &operator*=(const FEMatrix &matrix);

    // Matrix multiplication and return new matrix
    FEMatrix *operator*=(const FEMatrix &matrix) const;

    // Multiply self matrix by a constant
    FEMatrix &operator*=(double a);

    // Multiply matrix by a constant and return new matrix
    FEMatrix *operator*=(double a) const;

    // Create new matrix
    FEMatrix *clone() const;

    // Return max value of matrix
    double max() const;

    // Return min value of matrix
    double min() const;

    // Check if matrix is identity
    bool is_identity() const;

    // Check if matrix is symmetric
    bool is_symmetric() const;

    // Check if matrix is square nxn
    bool is_square() const;

    // Check if matrix is vector
    bool is_vector() const;

    // Make matrix symmetric
    void make_symmetric(bool upper);

    // Uses upper by default
    void make_symmetric();

    // Sum all matrix
    double sum() const;

    // Equal operator
    bool operator==(const FEMatrix &matrix) const;

    // Not equal operator
    bool operator!=(const FEMatrix &matrix) const;

    // Calculate determinant
    double det() const;

    // Get norm of vector
    double norm() const;

    // Check if matrix is diagonal
    bool is_diag() const;

    // Check if matrix is certain value
    bool is_double(double a) const;

    // Check if matrix is equal
    bool is_equal() const;

    // Check if matrix is only zeros
    bool is_zeros() const;

    // Check if matrix is only ones
    bool is_ones() const;

    // Transform matrix to string
    std::string to_string(bool matlab_like, std::string sep, bool to_int) const;

    // Transform matrix to string, to integer disabled
    std::string to_string(bool matlab_like, bool to_int) const;

    // Transform matrix to string, to integer disabled
    std::string to_string(bool matlab_like) const;

    // Transform matrix to string line separated by tab
    std::string to_string_line(bool to_int) const;

    // Transform matrix to string line separated by tab, integer disabled
    std::string to_string_line() const;

    // Set matrix name
    void set_name(std::string name);

    // Get matrix name
    std::string get_name() const;

    // Check if matrix is same as other
    bool equals(FEMatrix *mat) const;

};

#endif // FNELEM_GPU_MATH_FEMATRIX_H