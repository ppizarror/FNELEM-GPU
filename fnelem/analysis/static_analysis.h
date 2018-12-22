/**
FNELEM-GPU ANALYSIS - STATIC ANALYSIS
Performs static analysis calculation.

@package fnelem.model.base
@author ppizarror
@date 21/12/2018
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
#ifndef __FNELEM_ANALYSIS_STATIC_ANALYSIS_H
#define __FNELEM_ANALYSIS_STATIC_ANALYSIS_H

// Include headers
#include "../model/base/model.h"
#include "../math/matrix_inversion_cpu.h"
#include "../math/matrix_inversion_cuda.h"

// Library imports
#include <chrono>
#include <iostream>
#include <vector>

class StaticAnalysis {
private:

    // Stores model object
    Model *model;

    // Number of degrees of freedom
    int ndof = 0;

    // Matrix stiffness
    FEMatrix *Kt = nullptr;

    // Displacement vector
    FEMatrix *u = nullptr;

    // Force vector
    FEMatrix *F = nullptr;

    // Start dof numeration
    void define_dof();

    // Build stiffness matrix
    void build_stiffness_matrix();

    // Build force vector
    void build_force_vector();

    // Return yes/no
    std::string yes_no(bool v) const;

public:

    // Constructor
    explicit StaticAnalysis(Model *model);

    // Destructor
    ~StaticAnalysis();

    // Start analysis
    void analyze(bool use_gpu);

    // Return stiffness matrix
    FEMatrix *get_stiffness_matrix() const;

    // Return displacement vector
    FEMatrix *get_displacements_vector() const;

    // Return force vector
    FEMatrix *get_force_vector() const;

    // Get number of degrees of freedom
    int get_ndof() const;

    // Display analysis information to console
    void disp() const;

    // Clear data on demand
    void clear();

};

#endif // __FNELEM_ANALYSIS_STATIC_ANALYSIS_H