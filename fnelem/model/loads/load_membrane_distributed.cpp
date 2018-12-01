#include <utility>

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
    Load::apply(factor);
}

/**
 * Display load information.
 */
void LoadMembraneDistributed::disp() const {
    Load::disp();
}