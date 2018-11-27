/**
FNELEM-GPU NODE ELEMENT - NODE DEFINITION
Structural nodes.

@package fnelem.model.node
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

// Include header
#include "node.h"

/**
 * Destroy inner variables.
 */
void Node::destroy() {
    delete this->coords;
    delete this->displ;
    delete this->gdlid;
    delete this->loads;
    delete this->reaction;
}

/**
 * Destroy node.
 */
Node::~Node() {
    this->destroy();
};

/**
 * Creates a 2D node.
 *
 * @param tag Tag of the node
 * @param posx Position x
 * @param posy Position y
 */
Node::Node(std::string tag, double posx, double posy) : ModelComponent(std::move(tag)) {
    this->ngdl = 2;
    this->coords = FEMatrix_vector(this->ngdl);
    this->coords->set(0, posx);
    this->coords->set(1, posy);
    this->init();
}

/**
 * Creates a 2D node.
 *
 * @param tag Tag of the node
 * @param posx Position x
 * @param posy Position y
 */
Node::Node(std::string tag, double posx, double posy, double posz) : ModelComponent(std::move(tag)) {
    this->ngdl = 3;
    this->coords = FEMatrix_vector(this->ngdl);
    this->coords->set(0, posx);
    this->coords->set(1, posy);
    this->coords->set(2, posz);
    this->init();
}

/**
 * Init internal variables.
 */
void Node::init() {

    // Init ID of degrees of freedom
    this->gdlid = FEMatrix_vector(this->ngdl);

    // Init displacements
    this->displ = FEMatrix_vector(this->ngdl);

    // Init reactions
    this->reaction = FEMatrix_vector(this->ngdl);

    // Init node loads
    this->loads = FEMatrix_vector(this->ngdl);

    // Save initial values
    for (int i = 0; i < this->ngdl; i++) {
        this->gdlid->set(i, -1);
        this->displ->set(i, 0);
        this->reaction->set(i, 0);
        this->loads->set(i, 0);
    }

}

/**
 * Return number of degrees of freedom.
 *
 * @return Ngdl
 */
int Node::get_ngdl() const {
    return this->ngdl;
}

/**
 * Return a copy of node coordinates.
 *
 * @return Coordinates
 */
FEMatrix *Node::get_coordinates() const {
    return this->coords->clone();
}

/**
 * Get GDLID vector.
 *
 * @return
 */
FEMatrix *Node::get_gdlid() const {
    return this->gdlid->clone();
}

/**
 * Assign operator.
 *
 * @param node
 * @return
 */
Node &Node::operator=(const Node &node) {

    // Destroy inner vectors
    this->destroy();

    // Assign operator
    this->ngdl = node.ngdl;

    this->coords = node.coords;
    this->displ = node.displ;
    this->gdlid = node.gdlid;
    this->loads = node.loads;
    this->reaction = node.reaction;

    // Call ModelComponent assign
    ModelComponent::operator=(node);
    return *this;

}

/**
 * Get load results.
 *
 * @return
 */
FEMatrix *Node::get_load_results() const {
    return this->loads->clone();
}

/**
 * Get node displacements.
 *
 * @return
 */
FEMatrix *Node::get_displacements() const {
    return this->displ->clone();
}

/**
 * Get node reactions.
 *
 * @return
 */
FEMatrix *Node::get_reactions() const {
    return this->reaction->clone();
}

/**
 * Check vector.
 *
 * @param mat External matrix
 * @param vector_name Vector name
 */
void Node::check_vector(const FEMatrix *mat, std::string vector_name) {
    if (!mat->is_vector()) {
        throw std::logic_error("[NODE] " + vector_name + " must be a vector");
    }
    if (mat->length() != this->ngdl) {
        throw std::logic_error("[NODE] Number of degrees of freedom does not match");
    }
}

/**
 * Set GDLID of local degree of freedom.
 *
 * @param local_id ID of local freedom of liberty (1...n)
 * @param global_id ID of global freedom of liberty
 */
void Node::set_gdlid(int local_id, int global_id) {
    if (local_id < 1 || local_id > this->ngdl) {
        throw std::logic_error("[NODE] Local GDLID greather than number of Node NGDL");
    }
    this->gdlid->set(local_id - 1, global_id);
}

/**
 * Set GDLID from vector.
 *
 * @param gdlid Vector of GDLID vector
 */
void Node::set_gdlid(FEMatrix *gdl) {
    this->check_vector(gdl, "Node GDLID");
    (*this->gdlid) = gdl;
}

/**
 * Set node displacements.
 *
 * @param local_id ID of local freedom of liberty (1...n)
 * @param global_id Displacement
 */
void Node::set_displacement(int local_id, double d) {
    if (local_id < 1 || local_id > this->ngdl) {
        throw std::logic_error("[NODE] Local GDLID greather than number of Node NGDL");
    }
    this->displ->set(local_id - 1, d);
}

/**
 * Set node displacements from vector.
 *
 * @param d Vector of node displacements
 */
void Node::set_displacement(FEMatrix *d) {
    this->check_vector(d, "Node displacements");
    (*this->displ) = d;
}

/**
 * Apply load to node.
 *
 * @param load Node load
 */
void Node::apply_load(FEMatrix *load) {
    this->check_vector(load, "Node loads");
    (*this->reaction) -= *load;
}

/**
 * Apply element inner stress to node reactions.
 *
 * @param sigma Element stress
 */
void Node::apply_element_stress(FEMatrix *sigma) {
    this->check_vector(sigma, "Element stress");
    (*this->reaction) += *sigma;
}

/**
 * Display node information.
 */
void Node::disp() const {
    std::cout << "Node information" << std::endl;
    ModelComponent::disp();
    std::cout << "\n\tNumber degrees of freedom:\t" << this->ngdl << std::endl;
    std::cout << "\tCoordinates:\t" << this->coords->to_string_line() << std::endl;
    std::cout << "\tGLOBAL ID:\t\t" << this->gdlid->to_string_line(true) << std::endl;
    std::cout << "\tDisplacements:\t" << this->displ->to_string_line() << std::endl;
    std::cout << "\tReactions:\t\t" << this->reaction->to_string_line() << std::endl;
}

/**
 * Save node properties to file.
 *
 * @param file
 */
void Node::save_properties(std::ofstream &file) const {
    file << "\t Node " << this->get_model_tag() << ": ";
    file << this->coords->to_string_line() << std::endl;
}

/**
 * Save node displacements to file.
 *
 * @param file
 */
void Node::save_displacements(std::ofstream &file) const {
    file << "\t Node " << this->get_model_tag() << ": ";
    file << this->displ->to_string_line() << std::endl;
}

/**
 * Save node displacements to file.
 *
 * @param file
 */
void Node::save_reactions(std::ofstream &file) const {
    file << "\t Node " << this->get_model_tag() << ": ";
    file << this->reaction->to_string_line() << std::endl;
}

/**
 * Return x position.
 *
 * @return
 */
double Node::get_pos_x() const {
    return this->coords->get(0);
}

/**
 * Return y position.
 *
 * @return
 */
double Node::get_pos_y() const {
    return this->coords->get(1);
}

/**
 * Return z position.
 *
 * @return
 */
double Node::get_pos_z() const {
    if (this->ngdl == 2) {
        throw std::logic_error("[NODE] z-coordinate does not exist in a 2D node");
    }
    return this->coords->get(2);
}