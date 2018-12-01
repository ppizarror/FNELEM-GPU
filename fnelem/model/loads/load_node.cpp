/**
FNELEM-GPU LOAD - NODE LOAD.
Load applied to node.

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
#include "load_node.h"

/**
 * Constructor.
 *
 * @param tag Node load tag
 * @param n Target node
 * @param l Target load
 */
LoadNode::LoadNode(std::string tag, Node *n, FEMatrix *l) : Load(std::move(tag)) {

    // Stores node ndof
    this->ndof = n->get_ndof();

    // Check if load and node have same ndof
    if (l->length() != this->ndof || !l->is_vector()) {
        throw std::logic_error("[LOAD-NODE] Load vector invalid length, must be a vector, with node NODF components");
    }

    // Clones load
    this->load = l->clone();
    this->node = n;
}

/**
 * Destructor.
 */
LoadNode::~LoadNode() {
    delete this->load;
}

/**
 * Apply load.
 *
 * @param factor Load factor
 */
void LoadNode::apply(double factor) {

    // Generate new modified load
    FEMatrix *loadmod = this->load->clone();
    *loadmod *= factor;
    this->node->apply_load(loadmod);

    // Deletes created load
    delete loadmod;

}

/**
 * Display node load information.
 */
void LoadNode::disp() const {
    std::cout << "Load node information:" << std::endl;
    Load::disp();
    std::cout << "\n\tLoads:\t" << this->load->to_string_line() << std::endl;
}