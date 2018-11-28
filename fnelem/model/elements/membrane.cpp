#include <utility>

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

// Include header
#include "membrane.h"

/**
 * Membrane constructor.
 *
 * @param tag Membrane tag
 * @param n1 Node 1
 * @param n2 Node 2
 * @param n3 Node 3
 * @param n4 Node 4
 * @param E Elasticity modulus of the section
 * @param poisson Poisson constant of the section
 * @param thickness Thickness of the section
 */
Membrane::Membrane(std::string tag, Node *n1, Node *n2, Node *n3, Node *n4, double E, double poisson, double thickness)
        : Element(std::move(tag)) {

    // Stores material
    this->ndof = 8;
    this->E = E;
    this->poisson = poisson;

    // Stores nodes
    this->nodes->push_back(n1);
    this->nodes->push_back(n2);
    this->nodes->push_back(n3);
    this->nodes->push_back(n4);
    this->nnodes = 4;

    // Generate constitutive matrix
    this->constitutive = new FEMatrix(3, 3);
    this->constitutive->set(0, 0, 1 / (1 - poisson * poisson));
    this->constitutive->set(0, 1, poisson / (1 - poisson * poisson));
    this->constitutive->set(1, 0, poisson / (1 - poisson * poisson));
    this->constitutive->set(1, 1, 1 / (1 - poisson * poisson));
    this->constitutive->set(2, 2, 1 / (2 + 2 * poisson));
    (*this->constitutive) *= this->E;

    // Calculates dimension
    //     4 ------------- 3
    //     |       y       |
    //  2h |       #x      |
    //     |               |
    //     1 ------------- 2
    //            2b
    double db1, db2, dh1, dh2;
    db1 = fabs(n1->get_pos_x() - n2->get_pos_x()) / 2;
    db2 = fabs(n3->get_pos_x() - n4->get_pos_x()) / 2;
    dh1 = fabs(n1->get_pos_y() - n4->get_pos_y()) / 2;
    dh2 = fabs(n2->get_pos_y() - n3->get_pos_y()) / 2;

    if (fabs(db1 - db2) > __FEMATRIX_ZERO_TOL) {
        throw std::logic_error("[MEMBRANE] Invalid node dimension at @x");
    }
    if (fabs(dh1 - dh2) > __FEMATRIX_ZERO_TOL) {
        throw std::logic_error("[MEMBRANE] Invalid node dimension at @y");
    }

    // Stores geometry
    this->t = thickness;
    this->b = db1;
    this->h = dh1;

    // Init matrices
    this->Feq = FEMatrix_vector(8);
    this->stiffness_local = new FEMatrix(8, 8);
    this->stiffness_global = new FEMatrix(8, 8);
    this->dofid = FEMatrix_vector(8);

    // Set as initialized
    this->initialize();

    // Generates local matrix
    this->generate_local_stiffness();
    this->generate_global_stiffness();

    this->stiffness_local->set_disp_precision(4);
    this->stiffness_global->set_disp_precision(4);

}

/**
 * Destructor.
 */
Membrane::~Membrane() {
    delete this->Feq;
}

/**
 * Return membrane width.
 *
 * @return
 */
double Membrane::get_width() const {
    return 2 * this->b;
}

/**
 * Return membrane height.
 *
 * @return
 */
double Membrane::get_height() const {
    return 2 * this->h;
}

/**
 * Calculate local stiffness matrix
 */
void Membrane::generate_local_stiffness() {

    // Init A-vector
    FEMatrix *A = FEMatrix_vector(6);
    A->set_origin(1);

    // Calculate constitutive relationship
    A->set(1, (this->t * this->h * this->constitutive->get(0, 0)) / (6 * this->b));
    A->set(2, (this->t * this->b * this->constitutive->get(1, 1)) / (6 * this->h));
    A->set(3, (this->t * this->constitutive->get(0, 1)) / 4);
    A->set(4, (this->t * this->b * this->constitutive->get(2, 2)) / (6 * this->h));
    A->set(5, (this->t * this->h * this->constitutive->get(2, 2)) / (6 * this->b));
    A->set(6, (this->t * this->constitutive->get(2, 2)) / 4);

    // Generates local stiffness matrix
    this->stiffness_local->set_origin(1);
    this->stiffness_local->set(1, 1, 2 * this->k_aij(A, 1, 4));
    this->stiffness_local->set(1, 2, this->k_aij(A, 3, 6));
    this->stiffness_local->set(1, 3, this->k_cij(A, 4, 1));
    this->stiffness_local->set(1, 4, this->k_bij(A, 3, 6));
    this->stiffness_local->set(1, 5, -this->k_aij(A, 1, 4));
    this->stiffness_local->set(1, 6, -this->k_aij(A, 3, 6));
    this->stiffness_local->set(1, 7, this->k_cij(A, 1, 4));
    this->stiffness_local->set(1, 8, this->k_bij(A, 6, 3));

    this->stiffness_local->set(2, 2, 2 * this->k_aij(A, 2, 4));
    this->stiffness_local->set(2, 3, this->k_bij(A, 6, 3));
    this->stiffness_local->set(2, 4, this->k_cij(A, 2, 5));
    this->stiffness_local->set(2, 5, -this->k_aij(A, 3, 6));
    this->stiffness_local->set(2, 6, -this->k_aij(A, 2, 5));
    this->stiffness_local->set(2, 7, this->k_bij(A, 3, 6));
    this->stiffness_local->set(2, 8, this->k_cij(A, 5, 2));

    this->stiffness_local->set(3, 3, 2 * this->k_aij(A, 1, 4));
    this->stiffness_local->set(3, 4, -this->k_aij(A, 3, 6));
    this->stiffness_local->set(3, 5, this->k_cij(A, 1, 4));
    this->stiffness_local->set(3, 6, this->k_bij(A, 3, 6));
    this->stiffness_local->set(3, 7, -this->k_aij(A, 1, 4));
    this->stiffness_local->set(3, 8, this->k_aij(A, 6, 3));

    this->stiffness_local->set(4, 4, 2 * this->k_aij(A, 2, 4));
    this->stiffness_local->set(4, 5, this->k_bij(A, 6, 3));
    this->stiffness_local->set(4, 6, this->k_cij(A, 5, 2));
    this->stiffness_local->set(4, 7, this->k_aij(A, 3, 6));
    this->stiffness_local->set(4, 8, -this->k_aij(A, 5, 2));

    this->stiffness_local->set(5, 5, 2 * this->k_aij(A, 1, 4));
    this->stiffness_local->set(5, 6, this->k_aij(A, 3, 6));
    this->stiffness_local->set(5, 7, this->k_cij(A, 4, 1));
    this->stiffness_local->set(5, 8, -this->k_bij(A, 6, 3));

    this->stiffness_local->set(6, 6, 2 * this->k_aij(A, 2, 4));
    this->stiffness_local->set(6, 7, this->k_bij(A, 6, 3));
    this->stiffness_local->set(6, 8, this->k_cij(A, 2, 5));

    this->stiffness_local->set(7, 7, 2 * this->k_aij(A, 1, 4));
    this->stiffness_local->set(7, 8, -this->k_bij(A, 3, 6));

    this->stiffness_local->set(8, 8, 2 * this->k_aij(A, 2, 4));

    // Extends matrix transpose
    this->stiffness_local->make_symmetric();
    this->stiffness_local->set_origin(0);

    // Variable deletion
    delete A;

}

/**
 * Generate global stiffness matrix.
 */
void Membrane::generate_global_stiffness() {

    // As membrane does not have any rotation local and global matrices are the same
    (*this->stiffness_global) = this->stiffness_local;

}

/**
 * Return Aij stiffness modifier.
 *
 * @param A Matrix
 * @param i Position-i
 * @param j Position-j
 * @return
 */
double Membrane::k_aij(FEMatrix *A, int i, int j) const {
    return A->get(i) + A->get(j);
}

/**
 * Return Aij stiffness modifier.
 *
 * @param A Matrix
 * @param i Position-i
 * @param j Position-j
 * @return
 */
double Membrane::k_bij(FEMatrix *A, int i, int j) const {
    return A->get(i) - A->get(j);
}

/**
 * Return Cij stiffness modifier.
 *
 * @param A Matrix
 * @param i Position-i
 * @param j Position-j
 * @return
 */
double Membrane::k_cij(FEMatrix *A, int i, int j) const {
    return A->get(i) - 2 * A->get(j);
}

/**
 * Display membrane information.
 */
void Membrane::disp() const {
    std::cout << "Membrane information" << std::endl;
    Element::disp();

    // Display geometric / constitutive values
    std::cout << "\n\tWidth:\t\t\t\t" << 2 * this->b << std::endl;
    std::cout << "\tHeight:\t\t\t\t" << 2 * this->h << std::endl;
    std::cout << "\tElastic modulus:\t" << this->E << std::endl;
    std::cout << "\tPoisson modulus:\t" << this->poisson << std::endl;

    // Node information
    std::string nodetag;
    Node *n;
    for (unsigned long i = 0; i < this->nnodes; i++) {
        n = this->nodes->at(i);
        nodetag += n->get_model_tag();
        if (i < this->nnodes - 1) {
            nodetag += ", ";
        }
    }
    std::cout << "\tElement nodes:\t\t" << nodetag << std::endl;

    // Display local stiffnesss matrix
    std::cout << "\tLocal stiffness matrix (8x8):" << std::endl;
    this->stiffness_local->set_disp_identation(2);
    this->stiffness_local->disp();
    this->stiffness_local->set_disp_identation(0);

}

/**
 * Set ID degrees of freedom from node definition.
 */
void Membrane::set_dofid() {

    // Get nodes
    Node *n1 = this->nodes->at(0);
    Node *n2 = this->nodes->at(1);
    Node *n3 = this->nodes->at(2);
    Node *n4 = this->nodes->at(3);

    // Set dofid
    this->dofid->set(0, n1->get_dof(1));
    this->dofid->set(1, n1->get_dof(2));
    this->dofid->set(2, n2->get_dof(1));
    this->dofid->set(3, n2->get_dof(2));
    this->dofid->set(4, n3->get_dof(1));
    this->dofid->set(5, n3->get_dof(2));
    this->dofid->set(6, n4->get_dof(1));
    this->dofid->set(7, n4->get_dof(2));

}

/**
 * Validate (x,y) internal point to calculate stress/deformation.
 *
 * @param x X-position in membrane
 * @param y Y-position in membrane
 * @return
 */
void Membrane::validate_xy(double x, double y) const {
    if (fabs(x) > this->b || fabs(y) > this->h) {
        throw std::logic_error("[MEMBRANE] Position (x,y) out of membrane");
    }
}

/**
 * Get displacement vector from (x,y) point inside membrane.
 *
 * @param x X-position
 * @param y Y-position
 * @return
 */
FEMatrix *Membrane::get_displacement(double x, double y) const {

    // Check point is valid
    this->validate_xy(x, y);

    // Generate N matrix
    double N1 = (this->b - x) * (this->h - y) / (4 * this->b * this->h);
    double N2 = (this->b + x) * (this->h - y) / (4 * this->b * this->h);
    double N3 = (this->b + x) * (this->h + y) / (4 * this->b * this->h);
    double N4 = (this->b - x) * (this->h + y) / (4 * this->b * this->h);

    FEMatrix *N = new FEMatrix(2, 8);
    N->set(0, 0, N1);
    N->set(0, 2, N2);
    N->set(0, 4, N3);
    N->set(0, 6, N4);
    N->set(1, 1, N1);
    N->set(1, 3, N2);
    N->set(1, 5, N3);
    N->set(1, 7, N4);

    // Find node displacements
    Node *n1 = this->nodes->at(0);
    Node *n2 = this->nodes->at(1);
    Node *n3 = this->nodes->at(2);
    Node *n4 = this->nodes->at(3);

    FEMatrix *d = FEMatrix_vector(8);
    d->set(0, n1->get_displacement(1));
    d->set(1, n1->get_displacement(2));
    d->set(2, n2->get_displacement(1));
    d->set(3, n2->get_displacement(2));
    d->set(4, n3->get_displacement(1));
    d->set(5, n3->get_displacement(2));
    d->set(6, n4->get_displacement(1));
    d->set(7, n4->get_displacement(2));

    // Calculates deformation
    FEMatrix *dsp = (*N * *d);

    // Delete data
    delete N;
    delete d;

    // Return matrix
    return dsp;

}