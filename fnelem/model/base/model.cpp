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

/**
 * Constructor.
 *
 * @param tag Model tag
 * @param ndim Number of dimensions
 * @param ndof Number of degrees of freedom
 */
Model::Model(std::string tag, int ndim, int ndof) {
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