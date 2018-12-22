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

// Init header file
#ifndef __FNELEM_MODEL_RESTRAINTS_RESTRAINT_NODE_H
#define __FNELEM_MODEL_RESTRAINTS_RESTRAINT_NODE_H

// Include headers
#include "restraint.h"

// Library imports
#include <iostream>

class RestraintNode : public Restraint {
private:

    // Stores vector of DOFID restraints
    FEMatrix *dofid;

    // Stores node reference
    Node *node;

    // Stores NDOF of node
    int ndof = 0;

public:

    // Constructor
    RestraintNode(std::string tag, Node *n);

    // Destructor
    ~RestraintNode() override;

    // Add restraint
    void add_dofid(int id);

    // Add all DOF
    void add_all();

    // Apply restraints
    void apply() override;

    // Display restraint information to console
    void disp() const override;

};

#endif // __FNELEM_MODEL_RESTRAINTS_RESTRAINT_NODE_H