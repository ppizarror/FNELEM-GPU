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

// Include header
#include "fematrix.h"

/**
 * Default constructor, used by other class that initialize
 * pointer but no implicit constructor is called.
 */
FEMatrix::FEMatrix() {
    this->n = 0;
    this->m = 0;
    this->mat = new double[0];
}

/**
 * Creates matrix.
 *
 * @param n Number of rows
 * @param m Number of columns
 */
FEMatrix::FEMatrix(int n, int m) {
    if (n < 1 || m < 1) {
        throw std::logic_error("[FEMATRIX] Invalid matrix dimension");
    }
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
    if (n < 1 || m < 1) {
        throw std::logic_error("[FEMATRIX] Invalid matrix dimension");
    }
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
 * Fill matrix with a certain value.
 *
 * @param value Value to fill
 */
void FEMatrix::fill(double value) {
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            this->mat[i * this->m + j] = value;
        }
    }
}

/**
 * Fill matrix with zeros.
 */
void FEMatrix::fill_zeros() {
    this->fill(0.0);
}

/**
 * Fill matrix with ones.
 */
void FEMatrix::fill_ones() {
    this->fill(1.0);
}

void FEMatrix::disp_matrix(double *matrix, int dim_n, int dim_m) const {
    double disp_val; // Stores display value
    if (this->apply_pad) {
        int maxn = 0, snuml = 0;
        std::string snum;
        for (int i = 0; i < dim_n; i++) { // Rows
            for (int j = 0; j < dim_m; j++) { // Columns
                snuml = static_cast<int>(std::to_string(matrix[i * dim_m + j]).length());
                if (snuml > maxn) {
                    maxn = snuml;
                };
            }
        }
        for (int i = 0; i < dim_n; i++) { // Rows
            for (int j = 0; j < dim_m; j++) { // Columns
                disp_val = matrix[i * dim_m + j];
                if (abs(disp_val) < __FEMATRIX_ZERO_TOL) {
                    disp_val = 0;
                }
                std::cout << std::noshowpoint << std::setprecision(4) << std::setw(maxn) << disp_val;
                if (j < dim_m - 1) std::cout << " ";
            }
            std::cout << "" << std::endl;
        }
    } else {
        for (int i = 0; i < dim_n; i++) { // Rows
            for (int j = 0; j < dim_m; j++) { // Columns
                disp_val = matrix[i * dim_m + j];
                if (abs(disp_val) < __FEMATRIX_ZERO_TOL) {
                    disp_val = 0;
                }
                std::cout << std::noshowpoint << std::setprecision(4) << disp_val;
                if (j < dim_m - 1) std::cout << "\t";
            }
            std::cout << "" << std::endl;
        }
    }
    std::cout << "" << std::endl;
}

/**
 * Display matrix in console.
 */
void FEMatrix::disp() const {
    this->disp_matrix(this->mat, this->n, this->m);
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
 * Updates vector value.
 *
 * @param i Row position
 * @param j Column position
 * @param val Value
 */
void FEMatrix::set(int i, double val) {
    if (this->n == 1) { // Vector is row
        if (i - this->origin >= this->m) {
            throw std::logic_error("[FEMATRIX] Set column vector overflow");
        }
        this->_set(0, i - this->origin, val);
    } else if (this->m == 1) { // Vector is column
        if (i - this->origin >= this->n) {
            throw std::logic_error("[FEMATRIX] Set row vector overflow");
        }
        this->_set(i - this->origin, 0, val);
    } else {
        throw std::logic_error("[FEMATRIX] Matrix must be a vector");
    }
}

/**
 * Updates matrix value, no origin used.
 *
 * @param i Row position
 * @param j Column position
 * @param val Value
 */
void FEMatrix::_set(int i, int j, double val) {
    this->mat[i * this->m + j] = val;
}

/**
 * Returns matrix value, origin used.
 *
 * @param i Row position
 * @param j Column position
 * @return Value at matrix[i][j]
 */
double FEMatrix::get(int i, int j) const {
    if (i >= this->n + this->origin || j >= this->m + this->origin || (i - this->origin) < 0 ||
        (j - this->origin) < 0) {
        throw std::logic_error("[FEMATRIX] Get value from matrix olumn or row position overflow");
    }
    return this->mat[(i - this->origin) * this->m + (j - this->origin)];
}

/**
 * Returns vector value, origin used.
 *
 * @param i Row/Column position
 * @return Value at matrix[i][j]
 */
double FEMatrix::get(int i) const {
    if (this->n == 1) { // Vector is row
        if (i - this->origin >= this->m) {
            throw std::logic_error("[FEMATRIX] Get column vector overflow");
        }
        return this->mat[i - this->origin];
    } else if (this->m == 1) {
        if (i - this->origin >= this->n) {
            throw std::logic_error("[FEMATRIX] Get row vector overflow");
        }
        return this->mat[i - this->origin];
    } else {
        throw std::logic_error("[FEMATRIX] Matrix must be a vector");
    }
}

/**
 * Returns matrix value, no origin is used.
 *
 * @param i Row position
 * @param j Column position
 * @return Value at matrix[i][j]
 */
double FEMatrix::_get(int i, int j) const {
    return this->mat[i * this->m + j];
}

/**
 * Save matrix to file.
 *
 * @param filename Filename.
 */
void FEMatrix::save_to_file(std::string filename) const {
    std::ofstream plik;
    plik.open(filename);
    double save_val; // Save value
    for (int i = 0; i < this->n; i++) {
        for (int j = 0; j < this->m; j++) {
            save_val = this->_get(i, j);
            if (abs(save_val) < __FEMATRIX_ZERO_TOL) {
                save_val = 0;
            }
            plik << save_val;
            if (j < this->m - 1) {
                plik << "\t";
            }
        }

        if (i < this->n - 1) {
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
int *FEMatrix::size() const {
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

    // Create new matrix
    double *newMatrix = new double[matrix.n * matrix.m];
    this->n = matrix.n;
    this->m = matrix.m;

    // Assign values
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            newMatrix[i * this->m + j] = matrix._get(i, j);
        }
    }

    // Delete actual matrix and update
    delete[] this->mat;
    this->mat = newMatrix;

    // Return actual matrix
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
            this->mat[i * this->m + j] += matrix._get(i, j);
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
    FEMatrix newMatrix = this->clone();
    newMatrix += matrix;
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
            this->mat[i * this->m + j] -= matrix._get(i, j);
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
    FEMatrix newMatrix = this->clone();
    newMatrix -= matrix;
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
            newMatrix._set(i, j, -this->_get(i, j));
        }
    }
    return newMatrix;
}

/**
 * Matrix transpose.
 */
void FEMatrix::transpose_self() {

    // Save temporal dimension
    int tempdim = this->n;

    // Create new transposed matrix
    double *newMat = new double[this->n * this->m];
    for (int i = 0; i < this->m; i++) { // Rows
        for (int j = 0; j < this->n; j++) { // Columns
            newMat[i * this->n + j] = this->_get(j, i);
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

    // Delete
    delete[] newMat;

}

/**
 * Matrix transpose and return new.
 *
 * @return New transposed matrix
 */
FEMatrix FEMatrix::transpose() const {
    FEMatrix matrix = this->clone();
    matrix.transpose_self();
    return matrix;
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
    double *auxMatrix = new double[a * b];
    double sum = 0; // Stores partial sum

    // Multiply
    for (int i = 0; i < a; i++) { // Rows of new matrix
        for (int j = 0; j < b; j++) { // Columns of new matrix
            sum = 0;
            for (int k = 0; k < this->m; k++) {
                sum += this->_get(i, k) * matrix._get(k, j);
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

    // Delete matrix
    delete[] auxMatrix;

    // Return self
    return *this;

}

/**
 * Matrix multiplication and return new matrix.
 *
 * @param matrix Matrix to multiply
 * @return
 */
FEMatrix FEMatrix::operator*=(const FEMatrix &matrix) const {
    FEMatrix newMatrix = this->clone();
    newMatrix *= matrix;
    return newMatrix;
}

/**
 * Multiply self matrix by a constant.
 *
 * @param a Constant
 * @return
 */
FEMatrix &FEMatrix::operator*=(double a) {
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            this->_set(i, j, a * this->_get(i, j));
        }
    }
    return *this;
}

/**
 * Multiply matrix by a constant and return new matrix
 *
 * @param a Constant
 * @return
 */
FEMatrix FEMatrix::operator*=(double a) const {
    FEMatrix newMatrix = this->clone();
    newMatrix *= a;
    return newMatrix;
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
            newMatrix._set(i, j, this->_get(i, j));
        }
    }
    newMatrix.set_origin(this->origin_temp);
    return newMatrix;
}

/**
 * Set matrix origin.
 *
 * @param o New origin
 */
void FEMatrix::set_origin(int o) {
    if (o < 0) {
        throw std::logic_error("[FEMATRIX] Invalid origin");
    }
    this->origin = o;
    this->origin_temp = o;
}

/**
 * Return max value of matrix.
 *
 * @return
 */
double FEMatrix::max() const {
    double max = this->_get(0, 0);
    double num;
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            num = this->_get(i, j);
            if (num > max) {
                max = num;
            }
        }
    }
    return max;
}

/**
 * Return min value of matrix.
 *
 * @return
 */
double FEMatrix::min() const {
    double min = this->_get(0, 0);
    double num;
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            num = this->_get(i, j);
            if (num < min) {
                min = num;
            }
        }
    }
    return min;
}

/**
 * Disables origin.
 */
void FEMatrix::disable_origin() {
    this->origin = 0;
}

/**
 * Enables origin
 */
void FEMatrix::enable_origin() {
    this->origin = this->origin_temp;
}

/**
 * Check if matrix is identity.
 *
 * @return
 */
bool FEMatrix::is_identity() const {

    // If matrix is not square
    if (!this->is_square()) return false;

    // Check matrix
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            if (i == j) {
                if (fabs(this->_get(i, j) - 1) > __FEMATRIX_ZERO_TOL) {
                    return false;
                }
            } else {
                if (fabs(this->_get(i, j)) > __FEMATRIX_ZERO_TOL) {
                    return false;
                }
            }
        }
    }

    // Matrix is identity
    return true;

}

/**
 * Check if matrix is symmetric.
 *
 * @return
 */
bool FEMatrix::is_symmetric() const {

    // If matrix is not square
    if (!this->is_square()) return false;

    // Check matrix
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            if (i < j) {
                if (fabs(this->_get(i, j) - this->_get(j, i)) > __FEMATRIX_ZERO_TOL) {
                    return false;
                }
            }
        }
    }

    // Matrix is symmetric
    return true;

}

/**
 * Make matrix symmetric.
 *
 * @param upper Uses upper matrix to make symmetric matrix.
 * @return
 */
void FEMatrix::make_symmetric(bool upper) {

    // If matrix is not square
    if (!this->is_square()) {
        throw std::logic_error("[FEMATRIX] Cannot make symmetric a non-square matrix");
    }

    // Make symmetric
    if (upper) {
        for (int i = 0; i < this->n; i++) { // Rows
            for (int j = 0; j < this->m; j++) { // Columns
                if (i > j) {
                    this->_set(i, j, this->_get(j, i));
                }
            }
        }
    } else {
        for (int i = 0; i < this->n; i++) { // Rows
            for (int j = 0; j < this->m; j++) { // Columns
                if (i < j) {
                    this->_set(i, j, this->_get(j, i));
                }
            }
        }
    }

}

/**
 * Make matrix symmetric, uses upper by default.
 */
void FEMatrix::make_symmetric() {
    this->make_symmetric(true);
}

/**
 * Sum all matrix.
 *
 * @return
 */
double FEMatrix::sum() const {
    double st = 0;
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            st += this->_get(i, j);
        }
    }
    return st;
}

/**
 * Get matrix row.
 *
 * @param i Row number
 * @param from Init column position
 * @param to Final column position
 * @return
 */
FEMatrix FEMatrix::get_row(int i, int from, int to) const {

    // Check row consistency
    i -= this->origin;
    from -= this->origin;
    to -= this->origin;
    if (i > this->n) {
        throw std::logic_error("[FEMATRIX] Row position overflow");
    }
    if (from < 0 || from > this->m || to < 0 || to > this->m || to < from) {
        throw std::logic_error("[FEMATRIX] Column position overflow");
    }

    // Create new vector
    FEMatrix row = FEMatrix(1, to - from + 1);
    for (int j = from; j <= to; j++) { // Columns
        row.set(0, j - from, this->_get(i, j));
    }
    row.set_origin(this->origin_temp);

    // Return row vector
    return row;

}

/**
 * Get matrix full row.
 *
 * @param i Row number
 * @return
 */
FEMatrix FEMatrix::get_row(int i) const {
    return this->get_row(i, this->origin_temp, this->m);
}

/**
 * Get matrix column.
 *
 * @param j Column number
 * @param from Init row position
 * @param to Final row position
 * @return
 */
FEMatrix FEMatrix::get_column(int j, int from, int to) const {

    // Check column consistency
    j -= this->origin;
    from -= this->origin;
    to -= this->origin;
    if (j > this->m) {
        throw std::logic_error("[FEMATRIX] Column position overflow");
    }
    if (from < 0 || from > this->n || to < 0 || to > this->n || to < from) {
        throw std::logic_error("[FEMATRIX] Row position overflow");
    }

    // Create new vector
    FEMatrix column = FEMatrix(to - from + 1, 1);
    for (int i = from; i <= to; i++) { // Columns
        column.set(i - from, 0, this->_get(i, j));
    }
    column.set_origin(this->origin_temp);

    // Return column vector
    return column;

}

/**
 * Get matrix full column.
 *
 * @param j Column number
 * @return
 */
FEMatrix FEMatrix::get_column(int j) const {
    return this->get_column(j, this->origin_temp, this->n);
}

/**
 * Return max dimension of matrix, usefull for vectors.
 *
 * @return
 */
int FEMatrix::length() const {
    if (this->n > this->m) return this->n;
    return this->m;
}

/**
 * Equal operator.
 *
 * @param matrix Matrix to compare with
 * @return
 */
bool FEMatrix::operator==(const FEMatrix &matrix) const {
    if (this->n != matrix.n || this->m != matrix.m) return false;
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            if (fabs(this->_get(i, j) - matrix._get(i, j)) > __FEMATRIX_ZERO_TOL) {
                return false;
            }
        }
    }
    return true;
}

/**
 * Not equal operator.
 *
 * @param matrix Matrix to compare with
 * @return
 */
bool FEMatrix::operator!=(const FEMatrix &matrix) const {
    if (this->n != matrix.n || this->m != matrix.m) return true;
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            if (fabs(this->_get(i, j) - matrix._get(i, j)) < __FEMATRIX_ZERO_TOL) {
                return false;
            }
        }
    }
    return true;
}

/**
 * Calculates determinant of the matrix.
 *
 * @return
 */
double FEMatrix::det() const {
    if (!this->is_square()) {
        throw std::logic_error("[FEMATRIX] Cannot calculate determinant for a non-square matrix");
    }
    if (this->n == 1) {
        return this->_get(0, 0);
    }
    return this->_det_recursive(this->mat, this->n);
}

/**
 * Calculate determinant of the matrix.
 *
 * @param matrix Submatrix to calcualte determinant
 * @param d Actual dimension of the matrix
 * @return
 */
double FEMatrix::_det_recursive(double *matrix, int d) const {

    // If dimension is 2
    if (d == 2) {
        return matrix[0] * matrix[3] - matrix[1] * matrix[2];
    }

    // Dimension is greather than 2, split matrix into sub-matrices and evaluate
    int nd = d - 1;
    int i, j, k, p;
    double *submat = new double[nd * nd];

    // Iterates submatrices
    double dsum = 0; // Saves actual sum of the determinant
    int sign = 1; // Sign
    for (k = 0; k < d; k++) { // k stores column number

        // Fill submatrix
        for (i = 0; i < nd; i++) { // Rows
            for (j = 0; j < nd; j++) { // Column
                if (j >= k) p = j + 1; // Calculates position in global matrix
                else p = j;
                submat[i * nd + j] = matrix[(i + 1) * d + p];
            }
        }

        // Calculates recursively and multiply by column value at first row
        dsum += matrix[k] * this->_det_recursive(submat, nd) * sign;
        sign *= -1; // Flip sign

    }

    // Destroy memory
    delete[] submat;

    // Return value
    return dsum;

}

/**
 * Get norm of vector.
 *
 * @return
 */
double FEMatrix::norm() const {
    if (!this->is_vector()) {
        throw std::logic_error("[FEMATRIX] Matrix must be a vector");
    }
    double sumnorm = 0;
    for (int i = 0; i < this->length(); i++) {
        sumnorm += this->mat[i] * this->mat[i];
    }
    return sqrt(sumnorm);
}

/**
 * Check if matrix is vector.
 *
 * @return
 */
bool FEMatrix::is_vector() const {
    return this->n == 1 || this->m == 1;
}

/**
 * Check if matrix is diagonal.
 *
 * @return
 */
bool FEMatrix::is_diag() const {
    if (!this->is_square()) return false;
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            if (i != j) { // Out of diagonal, must be zero
                if (fabs(this->_get(i, j)) > __FEMATRIX_ZERO_TOL) return false;
            } else { // Diagonal, must be different than zero
                if (fabs(this->_get(i, j)) < __FEMATRIX_ZERO_TOL) return false;
            }
        }
    }
    return true;
}

/**
 * Check if all elements of matrix are same value.
 *
 * @param a Value to compare
 * @return
 */
bool FEMatrix::is_equal_double(double a) const {
    for (int i = 0; i < this->n; i++) { // Rows
        for (int j = 0; j < this->m; j++) { // Columns
            if (fabs(this->_get(i, j) - a) > __FEMATRIX_ZERO_TOL) return false;
        }
    }
    return true;
}

/**
 * Check if matrix only contains zeros.
 *
 * @return
 */
bool FEMatrix::is_zeros() const {
    return this->is_equal_double(0);
}

/**
 * Check if matrix only contains ones.
 * @return
 */
bool FEMatrix::is_ones() const {
    return this->is_equal_double(1);
}