/**
FNELEM-GPU BASE - MODEL
Model main class, integrates all components.

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
#ifndef __FNELEM_MODEL_BASE_MODEL_H
#define __FNELEM_MODEL_BASE_MODEL_H

// Include headers
#include "../nodes/node.h"
#include "../elements/element.h"
#include "../restraints/restraint.h"
#include "../loads/load_pattern.h"

// Library imports
#include <iostream>
#include <vector>

/**
 * Model class definition.
 */
class Model {
private:

    // Number of dimensions
    int ndim = 0;

    // Number of degrees of freedom
    int ndof = 0;

    // Nodes vector
    std::vector<Node *> *nodes = nullptr;

    // Elements vector
    std::vector<Element *> *elements = nullptr;

    // Restraints vector
    std::vector<Restraint *> *restraints = nullptr;

    // Load pattern vector
    std::vector<LoadPattern *> *loadpatterns = nullptr;

    // Write title header to file
    void write_file_title(std::ofstream &plik, std::string title) const;

    // Check if model is defined
    void check_defined(FEMatrix *u) const;

    // Check all variables are not null
    void check_non_null() const;

public:

    // Constructor
    Model(int ndim, int ndof);

    // Destructor
    ~Model();

    // Init model
    void initialize();

    // Set nodes
    void set_nodes(std::vector<Node *> *node);

    // Get nodes
    std::vector<Node *> *get_nodes() const;

    // Set elements
    void set_elements(std::vector<Element *> *element);

    // Get elements
    std::vector<Element *> *get_elements() const;

    // Set restraints
    void set_restraints(std::vector<Restraint *> *restraint);

    // Get restraints
    std::vector<Restraint *> *get_restrants() const;

    // Set load patterns
    void set_load_patterns(std::vector<LoadPattern *> *loadpattern);

    // Get load patterns
    std::vector<LoadPattern *> *get_load_patterns() const;

    // Apply model restraints
    void apply_restraints() const;

    // Apply load patterns
    void apply_load_patterns() const;

    // Update model components after solve method is done
    void update(FEMatrix *u);

    // Save results
    void save_results(std::string filename) const;

    // Display model information to console
    void disp() const;

    // Clear all elements of model
    void clear();

};

#endif // __FNELEM_MODEL_BASE_MODEL_H