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

// Include source
#include "static_analysis.h"
#include "../model/base/constants.h"

/**
 * Constructor.
 *
 * @param model Model definition
 */
StaticAnalysis::StaticAnalysis(Model *model) {
    this->model = model;
    this->ndof = 0;
}

/**
 * Destructor.
 */
StaticAnalysis::~StaticAnalysis() {
    if (this->ndof != 0) {
        delete this->Kt;
        delete this->u;
        delete this->F;
    }
}

/**
 * Start static analysis.
 *
 * @param use_gpu Use GPU inversion
 */
void StaticAnalysis::analyze(bool use_gpu) {

    // Define DOFID
    this->define_dof();
    return;

    // Apply load patterns
    this->model->apply_load_patterns();

    // Build analysis matrix data
    this->build_stiffness_matrix();
    this->build_force_vector();

    // Inverse matrix
    FEMatrix *invKt;
    if (!use_gpu) {
        invKt = matrix_inverse_cpu(this->Kt);
    } else {
        invKt = matrix_inverse_cpu(this->Kt);
    }

    // Solve matrix system
    this->u = *invKt * *this->F;

    // Update model
    this->model->update(this->u);

}

/**
 * Return matrix stiffness.
 *
 * @return
 */
FEMatrix *StaticAnalysis::get_stiffness_matrix() const {
    if (this->ndof == 0) {
        return nullptr;
    } else {
        return this->Kt->clone();
    }
}

/**
 * Return displacement vector.
 *
 * @return
 */
FEMatrix *StaticAnalysis::get_displacements_vector() const {
    if (this->ndof == 0) {
        return nullptr;
    } else {
        return this->u->clone();
    }
}

/**
 * Return force vector.
 *
 * @return
 */
FEMatrix *StaticAnalysis::get_force_vector() const {
    if (this->ndof == 0) {
        return nullptr;
    } else {
        return this->F->clone();
    }
}

/**
 * Return number of degrees of freedom.
 *
 * @return
 */
int StaticAnalysis::get_ndof() const {
    return this->ndof;
}

/**
 * Return yes/no.
 *
 * @param v Boolean variable
 * @return
 */
std::string StaticAnalysis::yes_no(bool v) const {
    if (v) {
        return "yes";
    } else {
        return "no";
    }
}

/**
 * Display static analysis on console.
 */
void StaticAnalysis::disp() const {
    std::cout << "Static analysis information:" << std::endl;
    if (this->ndof == 0) {
        std::cout << "\tAnalysis is not yet initialized" << std::endl;
        return;
    }

    std::cout << "\tStiffness matrix:" << std::endl;
    this->Kt->set_disp_identation(2);
    this->Kt->disp();
    this->Kt->set_disp_identation(0);
    std::cout << "\tStiffness determinant: " << this->Kt->det() << std::endl;
    std::cout << "\tStiffness symmetric: " << this->yes_no(this->Kt->is_symmetric()) << std::endl;

    std::cout << "\tForce vector:" << std::endl;
    this->F->set_disp_identation(2);
    this->F->disp();
    this->F->set_disp_identation(0);

    std::cout << "\tDisplacements vector:" << std::endl;
    this->u->set_disp_identation(2);
    this->u->disp();
    this->u->set_disp_identation(0);
}

/**
 * Start dof numbering.
 */
void StaticAnalysis::define_dof() {

    // Apply model restraints
    this->model->apply_restraints();

    // Check each node, then assign node DOFID
    std::vector<Node *> *nodes = this->model->get_nodes();
    FEMatrix *dof;
    int dof_count = 0;
    for (auto &node : *nodes) {
        dof = node->get_dof();
        for (int i = 0; i < dof->length(); i++) {
            if (fabs(dof->get(i) + 1) > FNELEM_CONST_ZERO_TOLERANCE) {
                dof_count += 1;
                node->set_dof(i + 1, dof_count);
            }
        }
    }
    this->ndof = dof_count;

    // Update element dofid
    std::vector<Element *> *elements = this->model->get_elements();
    for (auto &element : *elements) {
        element->set_dofid();
    }

}

/**
 * Build stiffness matrix.
 */
void StaticAnalysis::build_stiffness_matrix() {

    // Create stiffness matrix
    this->Kt = new FEMatrix(this->ndof, this->ndof);

    std::vector<Element *> *elements = this->model->get_elements();
    FEMatrix *dof;
    FEMatrix *Ktelem;
    int ndof, i, j;
    for (auto &element : *elements) {

        dof = element->get_dofid();
        dof->set_origin(1);
        ndof = element->get_ndof();
        Ktelem = element->get_stiffness_global();

        // Performs index method
        for (int r = 1; r <= ndof; r++) {
            for (int s = 1; s <= ndof; s++) {
                i = static_cast<int>(dof->get(r));
                j = static_cast<int>(dof->get(s));

                if (i != 0 && j != 0) {
                    this->Kt->set(i, j, this->Kt->get(i, j) + Ktelem->get(r, s));
                }
            }
        }

        delete Ktelem;
        delete dof;
    }

}

/**
 * Build force vector.
 */
void StaticAnalysis::build_force_vector() {

    // Create force
    this->F = FEMatrix_vector(this->ndof);

    FEMatrix *dof;
    int ndof;
    std::vector<Node *> *nodes = this->model->get_nodes();
    for (auto &node : *nodes) {
        dof = node->get_dof();
        ndof = node->get_ndof();
    }

}

/**
 * Clear data on demand.
 */
void StaticAnalysis::clear() {
    this->model->clear();
}