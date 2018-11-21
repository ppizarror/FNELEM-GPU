/**
FNELEM-GPU CPU MATRIX INVERSION
Performs matrix inversion using Gauss Jordan algorithm.

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

// Common class definition
class FEMatrix;

// Constant definition
#define _FEMATRIX_MIN_INVERSION_VALUE 0.0005

/**
 * Performs CPU matrix inversion using Gauss-Jordan elimination algorithm.
 *
 * @param matrix Matrix to inverse
 * @return Inverse matrix
 */
FEMatrix matrix_inverse_cpu(FEMatrix *matrix) {

    // Check matrix is square
    if (!matrix->is_square()) {
        throw std::logic_error("Matrix not square, cannot be inverted using Gauss-Jordan");
    }

    // Creates variables
    int dimension = matrix->get_square_dimension();
    double temporary, r;
    int i, j, k, temp;

    // Create augmented matrix
    matrix->disable_origin();
    int aug_n = dimension;
    int aug_m = 2 * dimension;
    double *augmentedMatrix = new double[aug_n * aug_m];
    for (int row = 0; row < dimension; row++) {
        for (int col = 0; col < 2 * dimension; col++) {
            if (col < dimension) {
                augmentedMatrix[row * aug_m + col] = matrix->get(row, col);
            } else {
                if (row == col % dimension) {
                    augmentedMatrix[row * aug_m + col] = 1;
                } else {
                    augmentedMatrix[row * aug_m + col] = 0;
                }
            }
        }
    }
    matrix->enable_origin();

    for (j = 0; j < dimension; j++) {
        temp = j;

        // Finding maximum jth column element in last (dimension-j) rows
        for (i = j + 1; i < dimension; i++)
            if (augmentedMatrix[i * aug_m + j] > augmentedMatrix[temp * aug_m + j])
                temp = i;

        if (fabs(augmentedMatrix[temp * aug_m + j]) < _FEMATRIX_MIN_INVERSION_VALUE) {
            std::cout << "[FEMatrix] Element are too small to deal with" << std::endl;
            break;
        }

        // Swapping row which has maximum jth column element
        if (temp != j) {
            for (k = 0; k < 2 * dimension; k++) {
                temporary = augmentedMatrix[j * aug_m + k];
                augmentedMatrix[j * aug_m + k] = augmentedMatrix[temp * aug_m + k];
                augmentedMatrix[temp * aug_m + k] = temporary;
            }
        }

        // Performing row operations to form required identity matrix out of the input matrix
        for (i = 0; i < dimension; i++)
            if (i != j) {
                r = augmentedMatrix[i * aug_m + j];
                for (k = 0; k < 2 * dimension; k++)
                    augmentedMatrix[i * aug_m + k] -=
                            (augmentedMatrix[j * aug_m + k] / augmentedMatrix[j * aug_m + j]) * r;
            } else {
                r = augmentedMatrix[i * aug_m + j];
                for (k = 0; k < 2 * dimension; k++)
                    augmentedMatrix[i * aug_m + k] /= r;
            }

    }

    // Stores inverse matrix
    double *invMatrix = new double[dimension * dimension];
    k = 0;
    for (i = 0; i < dimension; i++) {
        for (j = dimension; j < 2 * dimension; j++) {
            invMatrix[k] = augmentedMatrix[i * aug_m + j];
            k += 1;
        }
    }

    // Create matrix
    FEMatrix returnMatrix = FEMatrix(dimension, dimension, invMatrix);

    // Variable deletion
    delete[] augmentedMatrix;
    delete[] invMatrix;

    // Return matrix
    return returnMatrix;

}