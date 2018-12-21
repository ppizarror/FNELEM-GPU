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
Membrane::Membrane(std::string tag, Node *n1, Node *n2, Node *n3, Node *n4,
                   double E, double poisson, double thickness)
        : Element(std::move(tag)) {

    // Stores material
    this->ndof = 8;
    this->E = E;
    this->poisson = poisson;

    // Check nodes only are 2-D
    if (n1->get_ndof() != 2 || n2->get_ndof() != 2 || n3->get_ndof() != 2 || n4->get_ndof() != 2) {
        throw std::logic_error("[MEMBRANE] Membrane element only works with 2D nodes");
    }

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

    if (fabs(db1 - db2) > FNELEM_CONST_ZERO_TOLERANCE) {
        throw std::logic_error("[MEMBRANE] Invalid node dimension at @x");
    }
    if (fabs(dh1 - dh2) > FNELEM_CONST_ZERO_TOLERANCE) {
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
    std::cout << "Membrane information:" << std::endl;
    Element::disp();

    // Display geometric / constitutive values
    std::cout << "\n\tWidth:\t\t\t\t" << 2 * this->b << std::endl;
    std::cout << "\tHeight:\t\t\t\t" << 2 * this->h << std::endl;
    std::cout << "\tElastic modulus:\t" << this->E << std::endl;
    std::cout << "\tPoisson modulus:\t" << this->poisson << std::endl;

    // Node information
    std::string nodetag;
    Node *n;
    for (std::size_t i = 0; i < this->nnodes; i++) {
        n = this->nodes->at(i);
        nodetag += n->get_model_tag();
        if (i < this->nnodes - 1) {
            nodetag += ", ";
        }
    }
    std::cout << "\tElement nodes:\t\t" << nodetag << std::endl;

    // Display constitutive matrix
    std::cout << "\tConstitutive matrix (3x3):" << std::endl;
    this->constitutive->set_disp_identation(2);
    this->constitutive->disp();
    this->constitutive->set_disp_identation(0);

    // Display local stiffnesss matrix
    std::cout << "\tLocal stiffness matrix (8x8):" << std::endl;
    this->stiffness_local->set_disp_identation(2);
    this->stiffness_local->disp();
    this->stiffness_local->set_disp_identation(0);

    // Display equivalent force
    std::cout << "\tEquivalent force (1x8):" << std::endl;
    std::cout << "\t\t";
    for (int i = 0; i < 4; i++) {
        std::cout << std::setw(7) << "(" << (i + 1) << ")";
        if (i < 4) std::cout << "\t";
    }
    std::cout << std::endl;
    for (int i = 0; i < 2; i++) {
        std::cout << "\t\t";
        for (int j = 0; j < 4; j++) {
            std::cout << std::setw(8) << this->Feq->get(i + 2 * j);
            if (j < 3) std::cout << "\t";
        }
        std::cout << std::endl;
    }

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

    // Get node displacements
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

    // Calculates displacement; [2x8]x[8x1] = [2x1]
    FEMatrix *dsp = (*N * *d);

    // Delete data
    delete N;
    delete d;

    // Return matrix
    return dsp;

}

/**
 * Get deformation vector from (x,y) point inside membrane.
 *
 * @param x X-position
 * @param y Y-position
 * @return
 */
FEMatrix *Membrane::get_deformation(double x, double y) const {

    // Check point is valid
    this->validate_xy(x, y);

    // Generate B matrix
    double a1 = (this->b + x) / (4 * this->b * this->h);
    double a2 = (this->b - x) / (4 * this->b * this->h);
    double a3 = (this->h + y) / (4 * this->b * this->h);
    double a4 = (this->h - y) / (4 * this->b * this->h);

    FEMatrix *B = new FEMatrix(3, 8);
    B->set(0, 0, -a4);
    B->set(0, 2, a4);
    B->set(0, 4, a3);
    B->set(0, 6, -a3);
    B->set(1, 1, -a2);
    B->set(1, 3, -a1);
    B->set(1, 5, a1);
    B->set(1, 7, a2);
    B->set(2, 0, -a2);
    B->set(2, 1, -a4);
    B->set(2, 2, -a1);
    B->set(2, 3, a4);
    B->set(2, 4, a1);
    B->set(2, 5, a3);
    B->set(2, 6, a2);
    B->set(2, 7, -a3);

    // Get node displacements
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

    // Calculates deformation; [3x8]x[8x1] = [3x1]
    FEMatrix *def = (*B * *d);

    // Delete data
    delete B;
    delete d;

    // Return matrix
    return def;

}

/**
 * Calculate strain vector from (x,y) point inside membrane.
 *
 * @param x X-position
 * @param y Y-position
 * @return
 */
FEMatrix *Membrane::get_stress(double x, double y) const {
    FEMatrix *def = this->get_deformation(x, y);
    FEMatrix *stress = *this->constitutive * *def; // [3x3]x[3x1] = [3x1]
    delete def;
    return stress;
}

/**
 * Calculate local resistant force.
 *
 * @return
 */
FEMatrix *Membrane::get_force_local() const {

    // Get node displacements
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

    // Calculate force by multiplication with local stiffness matrix
    FEMatrix *fr_local = (*this->stiffness_local * *d);
    delete d;

    // Return matrix
    return fr_local;

}

/**
 * Calculate global resistant force.
 *
 * @return
 */
FEMatrix *Membrane::get_force_global() const {
    FEMatrix *fr_local = this->get_force_local();
    *fr_local -= this->Feq;
    return fr_local;
}

/**
 * Add equivalent force to nodes.
 *
 * @param nodenum Number of node (1-4)
 * @param f Force vector
 */
void Membrane::add_equivalent_force_node(int nodenum, FEMatrix *f) {
    int fnode = f->length();
    if (nodenum < 1 || nodenum > 4) {
        throw std::logic_error("[MEMBRANE] Invalid node number");
    }
    if (!f->is_vector()) {
        throw std::logic_error("[MEMBRANE] Load must be a vector");
    }
    int pos = fnode * nodenum - 1; // Force position
    f->disable_origin();
    for (int i = 0; i < fnode; i++) {
        this->Feq->set(pos + i - 1, this->Feq->get(pos + i - 1) + f->get(i));
    }
    f->enable_origin();
}

/**
 * Add resistant force to reaction.
 */
void Membrane::add_force_to_reaction() {

    // Calculate global resistant force
    FEMatrix *fr_global = this->get_force_global();

    // Get nodes
    Node *n1 = this->nodes->at(0);
    Node *n2 = this->nodes->at(1);
    Node *n3 = this->nodes->at(2);
    Node *n4 = this->nodes->at(3);

    // Generate load vector
    FEMatrix *load = FEMatrix_vector(2);

    load->set(0, -this->Feq->get(0));
    load->set(1, -this->Feq->get(1));
    n1->apply_load(load);

    load->set(0, -this->Feq->get(2));
    load->set(1, -this->Feq->get(3));
    n2->apply_load(load);

    load->set(0, -this->Feq->get(4));
    load->set(1, -this->Feq->get(5));
    n3->apply_load(load);

    load->set(0, -this->Feq->get(6));
    load->set(1, -this->Feq->get(7));
    n4->apply_load(load);

    // Add resistant force as load
    load->set(0, fr_global->get(0));
    load->set(1, fr_global->get(1));
    n1->apply_element_stress(load);

    load->set(0, fr_global->get(2));
    load->set(1, fr_global->get(3));
    n2->apply_element_stress(load);

    load->set(0, fr_global->get(4));
    load->set(1, fr_global->get(5));
    n3->apply_element_stress(load);

    load->set(0, fr_global->get(6));
    load->set(1, fr_global->get(7));
    n4->apply_element_stress(load);

    // Delete variables
    delete fr_global;
    delete load;

}

/**
 * Save properties to file.
 *
 * @param file File handler
 */
void Membrane::save_properties(std::ofstream &file) const {
    file << "\tMembrane " << this->get_model_tag() << ":";
    file << "\n\t\tWidth (2b):\t\t" << 2 * this->b;
    file << "\n\t\tHeight (2h):\t" << 2 * this->h;
    file << "\n\t\tThickness:\t\t" << this->t;
    file << "\n\t\tElastic mod:\t" << this->E;
    file << "\n\t\tPoisson mod:\t" << this->poisson;

    // Node information
    std::string nodetag;
    Node *n;
    for (std::size_t i = 0; i < this->nnodes; i++) {
        n = this->nodes->at(i);
        nodetag += n->get_model_tag();
        if (i < this->nnodes - 1) {
            nodetag += ", ";
        }
    }
    file << "\n\t\tElement nodes:\t" << nodetag;

}

/**
 * Generate tension matrix of integration points.
 *
 * @return
 */
FEMatrix *Membrane::generate_stress_npoints_matrix() const {

    // Calculate total elements
    int el = static_cast<int>(pow(FNELEM_CONST_MEMBRANE_INTEGRATION_NPOINTS + 1, 2));

    // Generate vector
    FEMatrix *tvec = new FEMatrix(el, 7);

    // Calculate integration point differential evaluation
    double dx = (2 * this->b) / (FNELEM_CONST_MEMBRANE_INTEGRATION_NPOINTS + 1);
    double dy = (2 * this->h) / (FNELEM_CONST_MEMBRANE_INTEGRATION_NPOINTS + 1);

    // Get first node coordinates
    double cglobx = this->nodes->at(0)->get_pos_x();
    double cgloby = this->nodes->at(0)->get_pos_y();

    int k = 0;
    double x = 0, y = 0; // Stores each integration point coordinates
    for (int i = 1; i < FNELEM_CONST_MEMBRANE_INTEGRATION_NPOINTS + 2; i++) {
        for (int j = 1; j < FNELEM_CONST_MEMBRANE_INTEGRATION_NPOINTS + 2; j++) {

            // Calculate integration point coordinates
            x = -this->b + (i - 1) * dx;
            y = -this->h + (j - 1) * dy;

            // Calculate tension
            FEMatrix *stress = this->get_stress(x, y);

            // Stores to matrix
            tvec->set(k, 0, cglobx + x + this->b);
            tvec->set(k, 1, cgloby + y + this->h);
            tvec->set(k, 2, x);
            tvec->set(k, 3, y);
            tvec->set(k, 4, stress->get(0));
            tvec->set(k, 5, stress->get(1));
            tvec->set(k, 6, stress->get(2));

            // Deletes matrix
            delete stress;
            k += 1;

        }
    }

    // Return matrix
    return tvec;

}

/**
 * Save internal stress to file.
 *
 * @param file
 */
void Membrane::save_internal_stress(std::ofstream &file) const {

    // Stores forces
    FEMatrix *fr = this->get_force_global();

    // Saves forces to each node
    file << "\tMembrane " << this->get_model_tag() << ":";
    file << "\n\t\tNode " << this->nodes->at(0)->get_model_tag() << " (-b, -h):\t" << fr->get(0) << ",\t" << fr->get(1);
    file << "\n\t\tNode " << this->nodes->at(1)->get_model_tag() << " (+b, -h):\t" << fr->get(2) << ",\t" << fr->get(3);
    file << "\n\t\tNode " << this->nodes->at(2)->get_model_tag() << " (+b, +h):\t" << fr->get(4) << ",\t" << fr->get(5);
    file << "\n\t\tNode " << this->nodes->at(3)->get_model_tag() << " (-b, +h):\t" << fr->get(6) << ",\t" << fr->get(7);

    // Writes tension
    FEMatrix *tm = this->generate_stress_npoints_matrix();
    int tmlen = tm->length();
    file << "\n\t\tStress " << this->get_model_tag() << " [GLOBALX GLOBALY X Y SIGMAX SIGMAY SIGMAXY DISPLX DISPLY]";
    for (int i = 0; i < tmlen; i++) {
        file << "\n\t\t\t" << tm->get(i, 0) << "\t" << tm->get(i, 1) << "\t" << tm->get(i, 2) << "\t";
        file << tm->get(i, 3) << "\t" << tm->get(i, 4) << "\t" << tm->get(i, 5) << "\t" << tm->get(i, 6);
    }

    // Deletes variables
    delete fr;
    delete tm;

}