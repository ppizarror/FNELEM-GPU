/**
FNELEM-GPU LOAD - LOAD PATTERN ABSTRACT CLASS.
General load pattern definition.

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

// Init header file
#ifndef __FNELEM_MODEL_LOADS_LOAD_PATTERN_H
#define __FNELEM_MODEL_LOADS_LOAD_PATTERN_H

// Include headers
#include "../base/model_component.h"

// Library imports
#include <iostream>

class LoadPattern : public ModelComponent {
public:

    // Constructor
    explicit LoadPattern(std::string tag);

    // Destructor
    ~LoadPattern();

    // Apply load pattern
    virtual void apply();

    // Display load pattern information to console
    void disp() const override;

};

#endif // __FNELEM_MODEL_LOADS_LOAD_PATTERN_H