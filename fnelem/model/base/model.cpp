/**
FNELEM-GPU BASE - MODEL
Model main class, integrates all components.

@package fnelem.model.base
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

// Include source
#include "model.h"
#include "constants.h"

/**
 * Constructor.
 *
 * @param ndim Number of dimensions
 * @param ndof Number of degrees of freedom
 */
Model::Model(int ndim, int ndof) {
    if (ndim < 1 || ndim > 3) {
        throw std::logic_error("[MODEL] Dimension number must be greater than one and lesser than 4");
    }
    if (ndof < 1) {
        throw std::logic_error("[MODEL] DOF must be greater than one");
    }
    this->ndim = ndim;
    this->ndof = ndof;
}

/**
 * Destructor.
 */
Model::~Model() = default;

/**
 * Init model.
 */
void Model::initialize() {
    this->check_non_null();
    for (auto &node : *this->nodes) {
        node->initialize();
    }
    for (auto &element : *this->elements) {
        element->initialize();
    }
}

/**
 * Set model nodes.
 *
 * @param node Node vector
 */
void Model::set_nodes(std::vector<Node *> *node) {
    this->nodes = node;
}

/**
 * Get model nodes.
 *
 * @return
 */
std::vector<Node *> *Model::get_nodes() const {
    return this->nodes;
}

/**
 * Set model elements.
 *
 * @param element Elements vector
 */
void Model::set_elements(std::vector<Element *> *element) {
    this->elements = element;
}

/**
 * Return model elements.
 * @return
 */
std::vector<Element *> *Model::get_elements() const {
    return this->elements;
}

/**
 * Set model restraints.
 *
 * @param restraint Restraints vector.
 */
void Model::set_restraints(std::vector<Restraint *> *restraint) {
    this->restraints = restraint;
}

/**
 * Return model restraints.
 *
 * @return
 */
std::vector<Restraint *> *Model::get_restrants() const {
    return this->restraints;
}

/**
 * Set model load patterns.
 *
 * @param loadpattern Load patterns vector
 */
void Model::set_load_patterns(std::vector<LoadPattern *> *loadpattern) {
    this->loadpatterns = loadpattern;
}

/**
 * Return model load pattern.
 *
 * @return
 */
std::vector<LoadPattern *> *Model::get_load_patterns() const {
    return this->loadpatterns;
}

/**
 * Apply model restraints.
 */
void Model::apply_restraints() const {
    for (auto &restraint : *this->restraints) {
        restraint->apply();
    }
}

/**
 * Apply model load patterns.
 */
void Model::apply_load_patterns() const {
    for (auto &loadpattern : *this->loadpatterns) {
        loadpattern->apply();
    }
}

/**
 * Update model after solve is done, needs node displacements vector.
 *
 * @param u
 */
void Model::update(FEMatrix *u) {

    // Check if model is full defined, if not throws exception
    this->check_defined(u);
    u->disable_origin();

    // Define displacement vector to each node
    FEMatrix *dof, *d;
    int node_ndof;
    for (Node *&node : *this->nodes) {
        dof = node->get_dofid();
        node_ndof = node->get_ndof();

        // Build displacement vector for each node
        d = FEMatrix_vector(node_ndof);
        for (int j = 0; j < node_ndof; j++) {
            if (dof->get(j) > 0) {
                d->set(j, u->get(static_cast<int>(dof->get(j)) - 1));
            }
        }

        // Save displacements
        node->set_displacement(d);

        // Delete vars
        delete d;
        delete dof;
    }
    u->enable_origin();

    // Add resistant forces to reactions
    for (Element *&element:*this->elements) {
        element->add_force_to_reaction();
    }

}

/**
 * Save results to file.
 *
 * @param filename Results filename
 */
void Model::save_results(std::string filename) const {
    std::ofstream plik;
    plik.open(filename);

    // Write software header
    plik << "FNELEM-GPU -  Finite element structural analysis using CUDA and GPU.\n";
    plik << "              v" << FNELEM_ABOUT_VERSION_V << " (";
    plik << FNELEM_ABOUT_VERSION_DATE << ") @ " << FNELEM_ABOUT_VERSION_AUTHOR << "\n";

    // Write model properties
    this->write_file_title(plik, "Input model properties:");

    // Write node information
    plik << "\nNodes:\n";
    plik << "\tNode count:\t" << this->nodes->size() << "\n";
    for (auto &node : *this->nodes) {
        node->save_properties(plik);
    }

    // Element properties
    plik << "\nElements:\n";
    plik << "\tElement count:\t\t" << this->elements->size() << "\n";
    for (auto &element : *this->elements) {
        element->save_properties(plik);
    }

    // Element results
    this->write_file_title(plik, "Analysis results:");

    // Save node results
    plik << "\nNode displacements:\n";
    for (auto &node : *this->nodes) {
        node->save_displacements(plik);
    }
    plik << "\nNode reactions:\n";
    for (auto &node : *this->nodes) {
        node->save_reactions(plik);
    }

    // Save element results
    plik << "\nElement stresses:";
    for (auto &element : *this->elements) {
        plik << "\n";
        element->save_internal_stress(plik);
    }

    // Close file
    plik.close();
}

/**
 * Write title header to file.
 *
 * @param plik File handle
 * @param title String title
 */
void Model::write_file_title(std::ofstream &plik, std::string title) const {
    plik << "\n--------------------------------------------------------------------\n";
    plik << title << "\n";
    plik << "--------------------------------------------------------------------\n";
}

/**
 * Display model information to console.
 */
void Model::disp() const {
    std::cout << "Model information:" << std::endl;
    std::cout << "\tNDIM Number of dimensions: " << this->ndim << std::endl;
    std::cout << "\tNDOF Number of degrees of freedom: " << this->ndof << std::endl;
}

/**
 * Check if model is defined.
 *
 * @param u Displacement vector
 */
void Model::check_defined(FEMatrix *u) const {
    this->check_non_null();
    if (!u->is_vector()) {
        throw std::logic_error("[MODEL] Displacement must be a vector");
    }
    if (!u->is_vector() || u->length() != this->ndof) {
        throw std::logic_error("[MODEL] Displacement NDOF must be the same as model NDOF");
    }
}

/**
 * Check all variables are non null.
 */
void Model::check_non_null() const {
    if (this->nodes == nullptr) {
        throw std::logic_error("[MODEL] Node vector must be defined");
    }
    if (this->elements == nullptr) {
        throw std::logic_error("[MODEL] Elements vector must be defined");
    }
    if (this->restraints == nullptr) {
        throw std::logic_error("[MODEL] Restraints vector must be defined");
    }
    if (this->loadpatterns == nullptr) {
        throw std::logic_error("[MODEL] LoadPatterns vector must be defined");
    }
}

/**
 * Clear model elements.
 */
void Model::clear() {
    for (auto &node : *this->nodes) {
        delete node;
    }
    for (auto &element : *this->elements) {
        delete element;
    }
    for (auto &restraint : *this->restraints) {
        delete restraint;
    }
    for (auto &loadpattern : *this->loadpatterns) {
        loadpattern->clear();
        delete loadpattern;
    }
}