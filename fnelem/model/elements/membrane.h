/**
FNELEM-GPU MEMBRANE ELEMENT
Bidimensional membrane element composed by 4 nodes.

@package fnelem.model.elements
@author ppizarror
@date 26/11/2018
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


#ifndef FNELEM_GPU_MEMBRANE_H
#define FNELEM_GPU_MEMBRANE_H

// Libray imports
#include "element.h"

// Constant definition
#define __MEMBRANE_INTEGRATION_NPOINTS 15

class Membrane : public Element {
private:

    // Elasticity modulus
    double E = 0;

    // Poisson modulus
    double poisson = 0;

    // Thickness of the section
    double t = 0;

    // Width of the section
    double b = 0;

    // Height of the section
    double h = 0;

    // Equivalent node forces
    FEMatrix *Feq;

    // Calculate local stiffness matrix
    void generate_local_stiffness();

    // Calculate global stiffness matrix
    void generate_global_stiffness();

    // Calculate Aij stiffness value
    double k_aij(FEMatrix *A, int i, int j) const;

    // Calculate Bij stiffness value
    double k_bij(FEMatrix *A, int i, int j) const;

    // Calculate Cij stiffness value
    double k_cij(FEMatrix *A, int i, int j) const;

    // Validate (x,y) point to perform stress/deformation analysis
    void validate_xy(double x, double y) const;

    // Generate stress vector
    FEMatrix *generate_stress_npoints_matrix() const;

public:

    // Constructor
    Membrane(std::string tag, Node *n1, Node *n2, Node *n3, Node *n4,
             double E, double poisson, double thickness);

    // Destructor
    ~Membrane();

    // Return membrane width
    double get_width() const;

    // Return membrane height
    double get_height() const;

    // Display membrane information
    void disp() const override;

    // Set degree of freedom ID from nodes
    void set_dofid() override;

    // Get displacement vector from (x,y) point inside membrane
    FEMatrix *get_displacement(double x, double y) const;

    // Get deformation/strain vector from (x,y) point inside membrane
    FEMatrix *get_deformation(double x, double y) const;

    // Calculate strain vector from (x,y) point inside membrane
    FEMatrix *get_stress(double x, double y) const;

    // Get local resistant force
    FEMatrix *get_force_local() const override;

    // Get global resistant force
    FEMatrix *get_force_global() const override;

    // Add equivalent force to internal forces
    void add_equivalent_force_node(int nodenum, FEMatrix *f);

    // Add resistant force to reaction
    void add_force_to_reaction() override;

    // Save properties to file
    void save_properties(std::ofstream &file) const override;

    // Save internal stress to file
    void save_internal_stress(std::ofstream &file) const override;

};

#endif // FNELEM_GPU_MEMBRANE_H