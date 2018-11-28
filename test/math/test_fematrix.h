/**
FNELEM-GPU - FEMATRIX TEST
Test finite element matrix (FEMATRIX).

@package test.math
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

// Include sources
#include "../test_utils.h"
#include "../../fnelem/math/fematrix.h"
#include "../../fnelem/math/matrix_inversion_cpu.h"

void __test_fematrix_init() {
    test_print_title("FEMATRIX", "test_fematrix_init");
    FEMatrix *matrix = new FEMatrix(3, 5);
    matrix->fill_zeros();
    assert(matrix->is_equal());
    int *dim = matrix->size();
    assert(dim[0] == 3);
    assert(dim[1] == 5);
    matrix->set(0, 0, 10);
    assert(matrix->get(0, 0) == 10);
    assert(!matrix->is_vector());
    delete[] dim;
    delete matrix;
}

void __test_fematrix_disp() {
    test_print_title("FEMATRIX", "test_fematrix_disp");
    FEMatrix *matrix = new FEMatrix(3, 3);
    matrix->fill_ones();
    matrix->set_disp_exponent(4);
    matrix->set_disp_precision(4);
    matrix->set_disp_identation(1);
    matrix->disp();
    assert(matrix->is_square());
    delete matrix;
}

void __test_fematrix_array() {
    test_print_title("FEMATRIX", "test_fematrix_array");
    FEMatrix *mat = new FEMatrix(2, 2);
    double *arr = mat->get_array();
    assert(arr[0] == 0);
    delete[] arr;
    delete mat;
}

void __test_fematrix_cpu_inversion() {
    test_print_title("FEMATRIX", "test_fematrix_cpu_inversion");
    int n = 3;
    auto *L = new double[n * n];
    L[0 * 3 + 0] = 1;
    L[0 * 3 + 1] = 2;
    L[0 * 3 + 2] = 3;
    L[1 * 3 + 0] = 5;
    L[1 * 3 + 1] = 2;
    L[1 * 3 + 2] = 1;
    L[2 * 3 + 0] = 2;
    L[2 * 3 + 1] = 2;
    L[2 * 3 + 2] = 3;
    FEMatrix *mat = new FEMatrix(n, n, L);
    mat->disp();
    mat->save_to_file("test-matrix-cpu.txt");
    FEMatrix *matInverse = matrix_inverse_cpu(mat);
    matInverse->save_to_file("test-matrix-cpu-inversion.txt");
    matInverse->disp();
    std::cout << matInverse->to_string(true) << std::endl;
    std::cout << matInverse->to_string(false) << std::endl;
    delete[] L;
    delete mat;
    delete matInverse;
}

void __test_fematrix_add() {
    test_print_title("FEMATRIX", "test_fematrix_add");
    FEMatrix *m1 = new FEMatrix(3, 3);
    m1->fill_ones();
    FEMatrix *m2 = new FEMatrix(3, 3);
    m2->fill_ones();
    m2->set(0, 0, 3);
    m2->set(1, 1, 5);
    m2->disp();

    // Adds
    *m2 += m1;
    assert(m2->get(0, 0) == 4);
    assert(m2->get(1, 1) == 6);

    // Assignation
    *m1 = m2;
    assert(m1->get(0, 0) == 4);
    assert(m1->get(1, 1) == 6);

    // Destroy
    delete m1;
    delete m2;
}

void __test_fematrix_substract() {
    test_print_title("FEMATRIX", "test_fematrix_substract");
    FEMatrix *m1 = new FEMatrix(3, 3);
    m1->fill_ones();
    FEMatrix *m2 = new FEMatrix(3, 3);
    m2->fill_ones();
    m2->set(0, 0, 3);
    m2->set(1, 1, 5);
    FEMatrix *m3 = *m2 - *m1;
    FEMatrix *m4 = -*m3;
    assert(m3->get(0, 0) == 2);
    assert(m4->get(0, 0) == -2);
    delete m1;
    delete m2;
    delete m3;
    delete m4;
}

void __test_fematrix_transpose() {
    test_print_title("FEMATRIX", "test_fematrix_transpose");
    FEMatrix *m1 = new FEMatrix(2, 3);
    m1->set(0, 0, 1);
    m1->set(0, 1, 1);
    m1->set(0, 2, 5);
    m1->set(1, 0, 9);
    m1->set(1, 1, 10);
    m1->set(1, 1, -1);
    m1->disp();
    m1->transpose_self();
    m1->disp();
    delete m1;
}

void __test_fematrix_multiplication() {
    test_print_title("FEMATRIX", "test_fematrix_multiplication");
    FEMatrix *m1 = new FEMatrix(2, 3);
    m1->set(0, 0, 1);
    m1->set(0, 1, 2);
    m1->set(0, 2, 3);
    m1->set(1, 0, 4);
    m1->set(1, 1, 5);
    m1->set(1, 2, 6);
    m1->disp();
    FEMatrix *m2 = m1->clone();
    m2->transpose_self();
    *m1 *= *m2;
    m1->disp();
    assert(m1->max() == 77);
    assert(m1->min() == 14);
    delete m1;
    delete m2;
}

void __test_fematrix_identity() {
    test_print_title("FEMATRIX", "test_fematrix_identity");
    FEMatrix *m = new FEMatrix(4, 4);
    m->set_origin(1);
    m->set(1, 1, 1);
    m->set(1, 2, 5);
    m->set(1, 3, 7);
    m->set(1, 4, 8);
    m->set(2, 1, 123);
    m->set(2, 2, 2432);
    m->set(2, 3, 22);
    m->set(2, 4, 2);
    m->set(3, 1, 45);
    m->set(3, 2, 345);
    m->set(3, 3, 65);
    m->set(3, 4, 7);
    m->set(4, 1, 2);
    m->set(4, 2, 8);
    m->set(4, 3, 22);
    m->set(4, 4, 34);
    m->disp();

    // Multiply by inverse
    FEMatrix *im = matrix_inverse_cpu(m);
    im->disp();
    *m *= *im;
    m->disp();

    // Check identity
    assert(m->is_identity());
    delete im;
    delete m;
}

void __test_fematrix_symmetric() {
    test_print_title("FEMATRIX", "test_fematrix_symmetric");
    FEMatrix *m = new FEMatrix(3, 3);
    m->set_origin(1);
    m->set(1, 1, 3);
    m->set(2, 2, 6);
    m->set(3, 3, 7);
    m->set(1, 2, 5);
    m->set(2, 1, 5);
    m->disp();
    assert(m->is_symmetric());
    m->transpose_self();
    assert(m->is_symmetric());
    delete m;
}

void __test_fematrix_make_symmetric() {
    test_print_title("FEMATRIX", "test_fematrix_make_symmetric");
    FEMatrix *m = new FEMatrix(3, 3);
    m->set_origin(1);
    m->set(1, 1, 3);
    m->set(2, 2, 6);
    m->set(3, 3, 7);
    m->set(1, 2, 5);
    m->set(2, 3, 8);
    m->disp();
    assert(!m->is_symmetric());
    m->make_symmetric();
    m->disp();
    assert(m->is_symmetric());
    delete m;
}

void __test_fematrix_constant_multiplication() {
    test_print_title("FEMATRIX", "test_fematrix_constant_multiplication");
    FEMatrix *m = new FEMatrix(5, 7);
    m->fill_ones();
    *m *= 5;
    assert(m->sum() == 5 * 7 * 5);
    delete m;
}

void __test_fematrix_row_column() {
    test_print_title("FEMATRIX", "test_fematrix_row_column");
    FEMatrix *m = new FEMatrix(4, 4);
    m->set_origin(1);
    m->set(1, 1, 1);
    m->set(1, 2, 2);
    m->set(1, 3, 3);
    m->set(1, 4, 4);
    m->set(2, 1, 5);
    m->set(2, 2, 6);
    m->set(2, 3, 7);
    m->set(2, 4, 8);
    m->set(3, 1, 9);
    m->set(3, 2, 10);
    m->set(3, 3, 11);
    m->set(3, 4, 12);
    m->set(4, 1, 13);
    m->set(4, 2, 14);
    m->set(4, 3, 15);
    m->set(4, 4, 16);
    m->disp();

    // Test row
    FEMatrix *row1 = m->get_row(1, 1, 4); // [1, 2, 3, 4]
    row1->disp();
    assert(row1->get(1) == 1 && row1->get(2) == 2 && row1->get(3) == 3 && row1->get(4) == 4);
    FEMatrix *row4 = m->get_row(4); // [13, 14, 15, 16]
    row4->disp();
    assert(row4->get(1) == 13 && row4->get(2) == 14 && row4->get(3) == 15 && row4->get(4) == 16);
    FEMatrix *r = m->get_row(3, 2, 2); // [9, 10, 11, 12]
    assert(r->get(1) == 10);

    FEMatrix *r1 = m->get_row(3, 1, 3); // [9, 10, 11, 12]
    r1->disp();
    assert(r1->length() == 3);

    // Test column
    FEMatrix *col1 = m->get_column(1); // [1, 5, 9, 13]
    col1->disp();
    assert(col1->length() == 4);
    assert(col1->is_vector());

    // Test multipication
    FEMatrix *colt = col1->transpose();
    colt->disp();
    *colt *= *col1;
    assert(colt->get(1) == 1 + 5 * 5 + 9 * 9 + 13 * 13);

    // Variable deletion
    delete m;
    delete col1;
    delete row1;
    delete row4;
    delete colt;
    delete r;
    delete r1;
}

void __test_fematrix_equal() {
    test_print_title("FEMATRIX", "test_fematrix_equal");
    FEMatrix *a = new FEMatrix(12, 4);
    a->fill(2);
    assert(a->is_equal());
    *a *= 0.5;
    FEMatrix *b = new FEMatrix(12, 4);
    b->fill_ones();
    assert(a->equals(b));
    *b *= 0.3;
    assert(!a->equals(b));
    delete a;
    delete b;
}

void __test_fematrix_determinant() {
    test_print_title("FEMATRIX", "test_fematrix_determinant");
    FEMatrix *mat1 = new FEMatrix(1, 1);
    mat1->set(0, 0, 3);
    assert(mat1->det() == 3); // Determinant 1x1

    // Determinant 2x2
    FEMatrix *mat2 = new FEMatrix(2, 2);
    mat2->set_origin(1);
    mat2->set(1, 1, 2);
    mat2->set(1, 2, 4);
    mat2->set(2, 1, 7);
    mat2->set(2, 2, 3);
    assert(is_num_equal(mat2->det(), -22));

    // Determinant 3x3
    FEMatrix *mat3 = new FEMatrix(3, 3);
    mat3->set(0, 0, 1);
    mat3->set(0, 1, 2);
    mat3->set(0, 2, 3);
    mat3->set(1, 0, 5);
    mat3->set(1, 1, 2);
    mat3->set(1, 2, 1);
    mat3->set(2, 0, 2);
    mat3->set(2, 1, 2);
    mat3->set(2, 2, 3);
    mat3->disp();
    assert(is_num_equal(mat3->det(), -4));

    // Determinant 4x4
    FEMatrix *mat4 = new FEMatrix(4, 4);
    mat4->set_origin(1);
    mat4->set(1, 1, 2);
    mat4->set(1, 2, 4);
    mat4->set(1, 3, 7);
    mat4->set(1, 4, 8);
    mat4->set(2, 1, 7);
    mat4->set(2, 2, 3);
    mat4->set(2, 3, 3);
    mat4->set(2, 4, 5);
    mat4->set(3, 1, 9);
    mat4->set(3, 2, 7);
    mat4->set(3, 3, 2);
    mat4->set(3, 4, 1);
    mat4->set(4, 1, 0);
    mat4->set(4, 2, 5);
    mat4->set(4, 3, 7);
    mat4->set(4, 4, 3);
    FEMatrix *mat4t = mat4->transpose();
    assert(is_num_equal(mat4->det(), -580));
    assert(is_num_equal(mat4t->det(), -580));

    // Determinant of inverse [4x4]
    FEMatrix *imat4 = matrix_inverse_cpu(mat4);
    assert(is_num_equal(imat4->det(), -0.00172413793103448208));

    // Determinant of sum
    FEMatrix *mat4s = *mat4 + *imat4;
    assert(is_num_equal(mat4s->det(), -1451.98793103448247165943));

    // Determinant zero for a NxN ones matrix
    FEMatrix *mat_ones = new FEMatrix(10, 10);
    mat_ones->fill_ones();
    assert(is_num_equal(mat_ones->det(), 0));

    // Destroy variables
    delete mat1;
    delete mat2;
    delete mat3;
    delete mat4;
    delete mat4t;
    delete imat4;
    delete mat4s;
    delete mat_ones;
}

void __test_fematrix_norm() {
    test_print_title("FEMATRIX", "test_fematrix_norm");
    FEMatrix *vector = new FEMatrix(6, 1);
    vector->set(0, 3);
    vector->set(1, 4);
    vector->set(2, 5);
    vector->set(3, 6);
    vector->set(4, 7);
    vector->set(5, 8);
    vector->disp();
    assert(is_num_equal(vector->norm(), 14.10673597966588488362));
    delete vector;
}

void __test_fematrix_diagonal() {
    test_print_title("FEMATRIX", "test_fematrix_diagonal");
    FEMatrix *nodiagonal = new FEMatrix(4, 2);
    assert(!nodiagonal->is_diag());
    FEMatrix *diagonal = new FEMatrix(4, 4);
    diagonal->set_origin(1);
    diagonal->set(1, 1, 1);
    diagonal->set(2, 2, 2);
    diagonal->set(3, 3, 3);
    diagonal->set(4, 4, 4);
    assert(diagonal->is_diag());
    diagonal->set(1, 2, -1);
    assert(!diagonal->is_diag());
    delete nodiagonal;
    delete diagonal;
}

void __test_fematrix_double_equal() {
    test_print_title("FEMATRIX", "test_fematrix_double_equal");
    FEMatrix *mat = new FEMatrix(5, 5);
    mat->fill_ones();
    assert(mat->is_double(1));
    assert(!mat->is_double(0));
    assert(mat->is_ones());
    mat->fill_zeros();
    assert(!mat->is_double(1));
    assert(mat->is_double(0));
    assert(mat->is_zeros());
    mat->set_name("MAT TEST");
    assert(mat->get_name() == "MAT TEST");
    delete mat;
}

void test_fematrix_suite() {
    __test_fematrix_init();
    __test_fematrix_disp();
    __test_fematrix_array();
    __test_fematrix_cpu_inversion();
    __test_fematrix_add();
    __test_fematrix_substract();
    __test_fematrix_transpose();
    __test_fematrix_multiplication();
    __test_fematrix_identity();
    __test_fematrix_symmetric();
    __test_fematrix_make_symmetric();
    __test_fematrix_constant_multiplication();
    __test_fematrix_row_column();
    __test_fematrix_equal();
    __test_fematrix_determinant();
    __test_fematrix_norm();
    __test_fematrix_diagonal();
    __test_fematrix_double_equal();
}