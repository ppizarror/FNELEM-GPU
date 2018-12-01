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

// Init header file
#ifndef __FNELEM_MODEL_LOADS_LOAD_MEMBRANE_DISTRIBUTED_H
#define __FNELEM_MODEL_LOADS_LOAD_MEMBRANE_DISTRIBUTED_H

// Include headers
#include "load.h"
#include "../elements/membrane.h"

// Library imports
#include <iostream>

class LoadMembraneDistributed : public Load {
private:

    // Loads
    double load1 = 0;
    double load2 = 0;

    // Load distance
    double dist1 = 0;
    double dist2 = 0;

    // Distance between loads
    double L = 0;

    // Angle of application
    double theta = 0;

    // Node references
    int nnode1 = 0;
    int nnode2 = 0;
    Node *node1;
    Node *node2;

    // Membrane reference
    Membrane *membrane;

    // Distributed load function
    double rho(double x) const;

    // Interpolation N1
    double N1(double x) const;

    // Interpolation N2
    double N3(double x) const;

    // X-displacement interpolation function
    double v1_int(double x) const;

    // Y-displacement interpolation function
    double v2_int(double x) const;

public:

    // Constructor
    LoadMembraneDistributed(std::string tag, Membrane *membrane, int node1,
                            int node2, double load1, double dist1, double load2, double dist2);

    // Destructor
    ~LoadMembraneDistributed();

    // Apply load
    void apply(double factor) override;

    // Display load information
    void disp() const override;

};

#endif // __FNELEM_MODEL_LOADS_LOAD_MEMBRANE_DISTRIBUTED_H