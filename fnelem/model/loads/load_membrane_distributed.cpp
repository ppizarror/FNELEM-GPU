/**
FNELEM-GPU LOAD - MEMBRANE DISTRIBUTED LOAD.
Distributed load in membrane element.

@package fnelem.model.load
@author ppizarror
@date 01/12/2018
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
#include "load_membrane_distributed.h"

/**
 * Constructor, need membrane object, node to apply load, load factors and distances from nodes,
 * node enumeration follows:
 *
 *  4 ------------- 3
 *  |               |
 *  |               |
 *  1 ------------- 2
 *
 * Only node combination 1-2, 2-3, 3-4, 4-1 are allowed.
 *
 * @param tag Load tag
 * @param membrane Membrane object reference
 * @param node1 Node position to apply load1
 * @param node2 Node position to apply load2
 * @param load1 Load 1
 * @param dist1 Distance from node 1 to apply load 1
 * @param load2 Load 2
 * @param dist2 Distance from node 1 to apply load 2
 */
LoadMembraneDistributed::LoadMembraneDistributed(std::string tag, Membrane *membrane, int node1,
                                                 int node2, double load1, double dist1, double load2, double dist2)
        : Load(std::move(tag)) {

    // Check nodes are well defined
    if (fabs(node2 - node1) == 2) {
        throw std::logic_error("[LOAD-MEMBRANE-DISTRIBUTED] Diagonal node definition is not allowed");
    }
    if (node1 < 1 || node1 > 4 || node2 < 1 || node2 > 4) {
        throw std::logic_error("[LOAD-MEMBRANE-DISTRIBUTED] Node position must be between 1 and 4");
    }
    if (node1 == node2) {
        throw std::logic_error("[LOAD-MEMBRANE-DISTRIBUTED] Node position cannot be the same");
    }

    // Check distances are well defined
    if (dist1 < 0 || dist1 >= 1 || dist2 <= 0 || dist2 > 1) {
        throw std::logic_error(
                "[LOAD-MEMBRANE-DISTRIBUTED] Load distances are not well defined, both must be between 0 and 1");
    }
    if (fabs(dist1 - dist2) < FNELEM_CONST_ZERO_TOLERANCE) {
        throw std::logic_error("[LOAD-MEMBRANE-DISTRIBUTED] Load distances cannot be the same");
    }

    // Get nodes from element
    std::vector<Node *> *nodes = membrane->get_nodes();
    Node *n1 = nodes->at(static_cast<unsigned long>(node1 - 1));
    Node *n2 = nodes->at(static_cast<unsigned long>(node2 - 1));

    // Check nodes only are 2-D
    if (n1->get_ndof() != 2 || n2->get_ndof() != 2) {
        throw std::logic_error("[LOAD-MEMBRANE] Membrane distributed load only works with 2D nodes");
    }

    // Calculates length
    double dx = n2->get_pos_x() - n1->get_pos_x();
    double dy = n2->get_pos_y() - n1->get_pos_y();
    this->L = sqrt(pow(dx, 2) + pow(dy, 2));

    // Calculates angle
    this->theta = atan(dy / dx);

    // Stores distances, loads and nodes
    this->membrane = membrane;
    this->nnode1 = node1;
    this->nnode2 = node2;
    this->node1 = n1;
    this->node2 = n2;
    this->load1 = load1;
    this->load2 = load2;
    this->dist1 = dist1 * this->L;
    this->dist2 = dist2 * this->L;

}

/**
 * Destructor.
 */
LoadMembraneDistributed::~LoadMembraneDistributed() = default;

/**
 * Apply load.
 *
 * @param factor Load factor
 */
void LoadMembraneDistributed::apply(double factor) {

    // Calculates integrals
    double steps = (this->dist2 - this->dist1) / FNELEM_CONST_GAUSS_INTEGRAL_POINTS;
    double v1 = 0.0, v2 = 0.0;
    for (int i = 0; i < FNELEM_CONST_GAUSS_INTEGRAL_POINTS; i++) {
        v1 += this->v1_int(this->dist1 + (i + 0.5) * steps) * steps;
        v2 += this->v2_int(this->dist1 + (i + 0.5) * steps) * steps;
    }

    // Rotate forces
    double v1_x = v1 * sin(this->theta);
    double v1_y = v1 * cos(this->theta);
    double v2_x = v2 * sin(this->theta);
    double v2_y = v2 * cos(this->theta);

    // Generate load vectors
    FEMatrix *loadvector1 = FEMatrix_vector(this->node1->get_ndof());
    FEMatrix *loadvector2 = FEMatrix_vector(this->node2->get_ndof());
    loadvector1->set(0, factor * v1_x);
    loadvector1->set(1, factor * v1_y);
    loadvector2->set(0, factor * v2_x);
    loadvector2->set(1, factor * v2_y);

    // Apply equivalent forces
    this->membrane->add_equivalent_force_node(this->nnode1, loadvector1);
    this->membrane->add_equivalent_force_node(this->nnode2, loadvector2);

    // Apply force to nodes
    this->node1->apply_load(loadvector1);
    this->node2->apply_load(loadvector2);

    // Variable deletion
    delete loadvector1;
    delete loadvector2;

}

/**
 * Display load information.
 */
void LoadMembraneDistributed::disp() const {
    std::cout << "Load membrane distribuited information:" << std::endl;
    Load::disp();
}

/**
 * Distributed load function.
 *
 * @param x Distance evaluation
 * @return
 */
double LoadMembraneDistributed::rho(double x) const {
    return this->load1 + (this->load2 - this->load1) * (x - this->dist1) / (this->dist2);
}

/**
 * Calculates N1 interpolation.
 *
 * @param x Distance
 * @return
 */
double LoadMembraneDistributed::N1(double x) const {
    return 1 - 3 * pow(x / this->L, 2) + 2 * pow(x / this->L, 3);
}

/**
 * Calculates N3 interpolation.
 *
 * @param x Distance
 * @return
 */
double LoadMembraneDistributed::N3(double x) const {
    return 3 * pow(x / this->L, 2) - 2 * pow(x / this->L, 3);
}

/**
 * X-displacement interpolation function.
 *
 * @param x Distance
 * @return
 */
double LoadMembraneDistributed::v1_int(double x) const {
    return this->rho(x) * this->N1(x);
}

/**
 * Y-displacement interpolation function.
 *
 * @param x Distance
 * @return
 */
double LoadMembraneDistributed::v2_int(double x) const {
    return this->rho(x) * this->N3(x);
}