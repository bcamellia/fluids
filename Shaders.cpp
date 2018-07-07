#include "Shaders.h"

namespace ShaderResource{

	namespace{
	
	const std::string particleVertShader = STRING(
	attribute vec3 inVertex;
	attribute vec4 inColor;
	varying vec4 outColor;
	uniform mat4 mPMVMatrix;
	uniform vec3 lightPos;
	void main(){
		gl_Position = mPMVMatrix * vec4(inVertex.xyz, 1.0);
		gl_PointSize = 1.0;
		outColor = inColor;
	}
	);

	const std::string particleFragShader = STRING(
		varying vec4 outColor;
		void main(){
			gl_FragColor = outColor;
	}
	);

	const std::string textVertShader = STRING(
	uniform sampler2D tex1;
	uniform vec2 pixel;
	attribute float agamma;
	attribute float ashift;
	varying float vshift;
	varying float vgamma;
	void main()
	{
		gl_FrontColor = gl_Color;
		gl_TexCoord[0].xy = gl_MultiTexCoord0.xy;
		gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
		vshift = ashift;
		vgamma = agamma;
	}
	);

	const std::string textFragShader = STRING(
		uniform sampler2D tex1;
	uniform vec2 pixel;
	varying float vgamma;
	varying float vshift;

	void main() {
		vec2 uv      = gl_TexCoord[0].xy;
		vec4 current = texture2D(tex1, uv);
		vec4 previous= texture2D(tex1, uv+vec2(-1,0)*pixel);
		vec4 next    = texture2D(tex1, uv+vec2(+1,0)*pixel);

		float r = current.r;
		float g = current.g;
		float b = current.b;

		if( vshift <= 0.333 )
		{
			float z = vshift/0.333;
			r = mix(current.r, previous.b, z);
			g = mix(current.g, current.r,  z);
			b = mix(current.b, current.g,  z);
		} 
		else if( vshift <= 0.666 )
		{
			float z = (vshift-0.33)/0.333;
			r = mix(previous.b, previous.g, z);
			g = mix(current.r,  previous.b, z);
			b = mix(current.g,  current.r,  z);
		}
		else if( vshift < 1.0 )
		{
			float z = (vshift-0.66)/0.334;
			r = mix(previous.g, previous.r, z);
			g = mix(previous.b, previous.g, z);
			b = mix(current.r,  previous.b, z);
		}

		vec3 color = pow( vec3(r,g,b), vec3(1.0/vgamma));
		gl_FragColor.rgb = color*gl_Color.rgb;
		gl_FragColor.a = (color.r+color.g+color.b)/3.0 * gl_Color.a;
	}
	);

	}

	void initShaders(){

		int attribs[10];
		int uniforms[10];
		
		attribs[0] = VERTEX_ARRAY;
		attribs[1] = COLOR_ARRAY;
		uniforms[0] = eMatrix;
		uniforms[1] = eLightPos;

		shaderConstructionInfo sci(attribs,2,uniforms,2);
		sci.vertShaderSource = particleVertShader;
		sci.fragShaderSource = particleFragShader;

		Renderer::loadShader(sci, PARTICLE_SHADER);


		/*uniforms[0] = eTex1;
		uniforms[1] = ePixel;

		sci = shaderConstructionInfo(0,0,uniforms,2);
		sci.vertShaderSource = textVertShader;
		sci.fragShaderSource = textFragShader;

		Renderer::loadShader(sci, TEXT_SHADER);*/

	}

};