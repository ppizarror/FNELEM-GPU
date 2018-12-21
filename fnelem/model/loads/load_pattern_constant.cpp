#include <utility>

/**
FNELEM-GPU LOAD - CONSTANT LOAD PATTERN.
Constant load pattern class.

@package fnelem.model.load
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

// Include header
#include "load_pattern_constant.h"

/**
 * Load pattern constant constructor.
 *
 * @param tag Load pattern constant tag
 * @param loads Load vector
 */
LoadPatternConstant::LoadPatternConstant(std::string tag, std::vector<Load *> *loads) :
        LoadPattern(std::move(tag)) {
    this->load_array = loads;
}

/**
 * Destructor.
 */
LoadPatternConstant::~LoadPatternConstant() = default;

/**
 * Apply load pattern.
 */
void LoadPatternConstant::apply() {
    for (auto &i : *this->load_array) {
        i->apply(1);
    }
}

/**
 * Display load pattern information.
 */
void LoadPatternConstant::disp() const {
    std::cout << "Constant load pattern information:" << std::endl;
    LoadPattern::disp();
    std::cout << "\n\tTotal loads at vector:\t" << this->load_array->size() << std::endl;
}