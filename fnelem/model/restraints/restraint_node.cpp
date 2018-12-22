/**
FNELEM-GPU RESTRAINTS - NODE RESTRAINT.
Node restraint definition.

@package fnelem.model.restraints
@author ppizarror
@date 30/11/2018
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
#include "restraint_node.h"

/**
 * Constructor.
 *
 * @param tag Node restraint tag
 * @param n Node to apply restraints
 */
RestraintNode::RestraintNode(std::string tag, Node *n) : Restraint(std::move(tag)) {

    // Generate inner restraints vector
    this->dofid = FEMatrix_vector(n->get_ndof());
    this->dofid->fill(-1); // IF -1 no dofid has been filled

    // Stores node reference
    this->node = n;
    this->ndof = n->get_ndof();

}

/**
 * Destructor.
 */
RestraintNode::~RestraintNode() {
    delete this->dofid;
}

/**
 * Add DOFID restraint.
 *
 * @param id
 */
void RestraintNode::add_dofid(int id) {

    // Check if id is valid
    if (id < 1 || id > this->ndof) {
        throw std::logic_error("[NODE-RESTRAINT] Local DOFID restraint greather than number of Node NDOF");
    }

    // Stores data
    this->dofid->set(id - 1, id);

}

/**
 * Apply node restraints.
 */
void RestraintNode::apply() {
    for (int i = 0; i < this->ndof; i++) {
        if (fabs(this->dofid->get(i) + 1) > FNELEM_CONST_ZERO_TOLERANCE) {
            this->node->set_dof(static_cast<int>(this->dofid->get(i)), -1);
        }
    }
}

/**
 * Display node restraint information.
 */
void RestraintNode::disp() const {
    std::cout << "Node restraint information:" << std::endl;
    Restraint::disp();
    std::cout << "\n\tRestrained node:\t" << this->node->get_model_tag();

    // Generate restrained DOFID
    std::cout << "\n\tRestrained DOFID:\t";
    if (this->dofid->is_double(-1)) {
        std::cout << "NONE";
    } else {
        std::string resdof;
        for (int i = 0; i < this->ndof; i++) {
            if (fabs(this->dofid->get(i) + 1) > FNELEM_CONST_ZERO_TOLERANCE) {
                resdof += std::to_string(static_cast<int>(this->dofid->get(i)));
                if (i < this->ndof - 1) {
                    resdof += "\t";
                }
            }
        }
        std::cout << resdof;
    }
    std::cout << std::endl;

}

/**
 * Add all dof.
 */
void RestraintNode::add_all() {
    for (int i = 1; i <= this->ndof; i++) {
        this->add_dofid(i);
    }
}