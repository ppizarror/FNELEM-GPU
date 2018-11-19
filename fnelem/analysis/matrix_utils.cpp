/**
FNELEM-GPU MATRIX UTILITARY FUNCTIONS
Performs matrix inversion using Gauss Jordan algorithm.
Based on: https://github.com/ZhengzhongSun/Matrix-Inversion-with-CUDA

@package fnelem.analysis
@author ppizarror
@date 19/11/2018
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
#include <string>
#include <iostream>

/**
 * Save matrix to file.
 *
 * @param A Matrix
 * @param s File name
 * @param n Number of rows
 * @param h Number of columns
 */
void save_matrix_to_file(double *A, std::string s, int n, int h) {
    std::ofstream plik;
    plik.open(s);
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < h; i++) {
            plik << A[j * n + i] << "\t";
        }
        plik << std::endl;
    }
    plik.close();
}