#include <utility>

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

// Include source
#include "model_component.h"

/**
 * Init model.
 */
ModelComponent::ModelComponent() = default;

/**
 * Init model with tag.
 *
 * @param tag
 */
ModelComponent::ModelComponent(std::string tag) {
    this->tagID = std::move(tag);
}

/**
 * Returns model tag.
 *
 * @return Model tag.
 */
std::string ModelComponent::get_model_tag() {
    return this->tagID;
}

/**
 * Object destruction.
 */
ModelComponent::~ModelComponent() = default;