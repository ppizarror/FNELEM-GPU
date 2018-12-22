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

// Init header file
#ifndef __FNELEM_MODEL_NODES_NODE_H
#define __FNELEM_MODEL_NODES_NODE_H

// Include headers
#include "../base/model_component.h"
#include "../../math/fematrix.h"
#include "../../math/fematrix_utils.h"

// Library imports
#include <iostream>

/**
 * Node element.
 */
class Node : public ModelComponent {
private:

    // Number of degres of freedom
    int ndof = 0;

    // ID of degrees of freedom
    FEMatrix *dofid;

    // Coordinates of the node
    FEMatrix *coords;

    // Displacements vector
    FEMatrix *displ;

    // Loads vector
    FEMatrix *loads;

    // Reaction vector
    FEMatrix *reaction;

    // Init internal variables
    void init();

    // Destroy malloc inner variables
    void destroy();

    // Check external vector
    void check_vector(const FEMatrix *mat, std::string vector_name);

public:

    // Destroy node
    ~Node();

    // 2D node
    Node(std::string tag, double posx, double posy);

    // 3D node
    Node(std::string tag, double posx, double posy, double posz);

    // Assign
    Node &operator=(const Node &node);

    // Return number of degrees of freedom
    int get_ndof() const;

    // Return node coordinates
    FEMatrix *get_coordinates() const;

    // Get degrees of freedom of node
    FEMatrix *get_dof() const;

    // Get results of loads
    FEMatrix *get_load_results() const;

    // Get node displacements
    FEMatrix *get_displacements() const;

    // Get node reactions
    FEMatrix *get_reactions() const;

    // Set node degrees of freedom
    void set_dof(int local_id, int global_id);

    // Get node degree of freedom
    int get_dof(int local_id);

    // Set node degrees of freedom from vector/matrix
    void set_dof(FEMatrix *gdl);

    // Set node displacements
    void set_displacement(int local_id, double d);

    // Get node displacement at local dof
    double get_displacement(int local_id);

    // Set node displacements from vector/matrix
    void set_displacement(FEMatrix *d);

    // Apply load to node
    void apply_load(FEMatrix *load);

    // Adds element inner stress to reactions
    void apply_element_stress(FEMatrix *sigma);

    // Get node reaction ad local dof
    double get_reaction(int local_id);

    // Display node information to console
    void disp() const override;

    // Save properties to file
    void save_properties(std::ofstream &file) const;

    // Save properties to file
    void save_displacements(std::ofstream &file) const;

    // Save properties to file
    void save_reactions(std::ofstream &file) const;

    // Get position x
    double get_pos_x() const;

    // Get position x
    double get_pos_y() const;

    // Get position x
    double get_pos_z() const;

    // Initialize element
    virtual void initialize();

};

#endif // __FNELEM_MODEL_NODES_NODE_H