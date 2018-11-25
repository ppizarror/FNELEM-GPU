/**
FNELEM-GPU BASE ELEMENTS - MODEL COMPONENT
Base element of the platform.

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

// Init header file
#ifndef BASE_MODEL_COMPONENT_HEADER_FILE
#define BASE_MODEL_COMPONENT_HEADER_FILE

// Library import
#include <iostream>
#include <string>

/**
 * Base model component.
 */
class ModelComponent {
private:

    // Element tag
    std::string tagID = "";

public:

    // Init model
    ModelComponent();

    // Init model with tag
    explicit ModelComponent(std::string tag);

    // Destroy object
    ~ModelComponent();

    // Returns tag
    std::string get_model_tag() const;

    // Assign operator
    ModelComponent &operator=(const ModelComponent &model);

    // Display model
    void disp() const;

};

#endif // BASE_MODEL_COMPONENT_HEADER_FILE