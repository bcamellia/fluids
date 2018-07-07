#ifndef SHADERSCB_H
#define SHADERSCB_H

#include <iostream>
#include <string>
#include "Renderer.h"

#define STRING(A) std::string(#A);
#define PARTICLE_SHADER "particleShader"

namespace ShaderResource{

	void initShaders();

};

#endif