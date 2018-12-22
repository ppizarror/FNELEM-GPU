/**
FNELEM-GPU - TEST STATIC ANALYSIS
Test static analysis class.

@package test.analysis
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

// Include sources
#include "../test_utils.h"
#include "../../fnelem/analysis/static_analysis.h"
#include "../../fnelem/model/base/model.h"
#include "../../fnelem/model/elements/membrane.h"
#include "../../fnelem/model/restraints/restraint_node.h"
#include "../../fnelem/model/loads/load_pattern_constant.h"
#include "../../fnelem/model/loads/load_node.h"

void __test_static_analysis_test1() {
    test_print_title("STATIC-ANALYSIS", "test_static_analysis_test1");

    // Simple membrane test
    // Example from book
    //   INTRODUCCION AL ANALISIS ESTRUCTURAL POR ELEMENTOS FINITOS
    //   Author: JORGE EDUARDO HURTADO GÓMEZ
    //   http://bdigital.unal.edu.co/10002/6/958932276X.2002.pdf
    //   Page 92
    //
    //           1000 kN
    //             |
    //             v
    //    2 ------ 4 ------ 6
    //    |        |        |
    //    |   (1)  |   (2)  |  2m
    //    |        |        |
    //    1 ------ 3 ------ 5
    //    ^   2m       2m   ^
    // ==========================

    // Create model
    Model *model = new Model(2, 8);

    // Create nodes
    std::vector<Node *> *nodes = new std::vector<Node *>();
    nodes->push_back(new Node("N1", 0, 0));
    nodes->push_back(new Node("N2", 0, 2));
    nodes->push_back(new Node("N3", 2, 0));
    nodes->push_back(new Node("N4", 2, 2));
    nodes->push_back(new Node("N5", 4, 0));
    nodes->push_back(new Node("N6", 4, 2));

    // Creates membrane
    std::vector<Element *> *elements = new std::vector<Element *>();
    elements->push_back(new Membrane("MEM1", nodes->at(0), nodes->at(2), nodes->at(3),
                                     nodes->at(1), 2000, 0.2, 1.0));
    elements->push_back(new Membrane("MEM2", nodes->at(2), nodes->at(4), nodes->at(5),
                                     nodes->at(3), 2000, 0.2, 1.0));

    // Create restraints
    std::vector<Restraint *> *restraints = new std::vector<Restraint *>();
    RestraintNode *r1 = new RestraintNode("R1", nodes->at(0));
    RestraintNode *r2 = new RestraintNode("R2", nodes->at(4));
    r1->add_all();
    r2->add_all();
    restraints->push_back(r1);
    restraints->push_back(r2);

    // Create loads
    std::vector<Load *> *loads = new std::vector<Load *>();
    FEMatrix *loadv = FEMatrix_vector(2);
    loadv->set(1, -1000);
    loads->push_back(new LoadNode("P", nodes->at(3), loadv));

    // Create load pattern
    std::vector<LoadPattern *> *loadpattern = new std::vector<LoadPattern *>();
    loadpattern->push_back(new LoadPatternConstant("LoadConstant", loads));

    // Define elements
    model->set_nodes(nodes);
    model->set_elements(elements);
    model->set_restraints(restraints);
    model->set_load_patterns(loadpattern);

    // Create analysis
    StaticAnalysis *analysis = new StaticAnalysis(model);
    analysis->disp();

    // Check all results
    assert(analysis->get_ndof() == 0);
    assert(analysis->get_displacements_vector() == nullptr);
    assert(analysis->get_force_vector() == nullptr);
    assert(analysis->get_stiffness_matrix() == nullptr);

    // Run analysis
    analysis->analyze(false);
    analysis->disp();

    // Save results to file
    model->save_results("out/test-static-analysis-1.txt");

    // Delete data
    analysis->clear();
    delete elements;
    delete loadpattern;
    delete loadv;
    delete loads;
    delete restraints;
    delete nodes;
    delete model;
    delete analysis;
}

/**
 * Performs TEST-STATIC-ANALYSIS suite.
 */
void test_static_analysis_suite() {
    __test_static_analysis_test1();
}