#ifndef RENDERER_H

#define RENDERER_H

#define _USE_MATH_DEFINES
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#include <math.h>
#include <Windows.h>
#include <vector>
#include <cassert>
#include <string>
#include <time.h>
#include <GL/glew.h>
#include <GL/wglew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <boost/unordered_map.hpp>

#include <direct.h>
#define GetCurrentDir _getcwd

#define WINDOW_TITLE "Fluid Sim"
#define WINDOW_HEIGHT 480.0f
#define WINDOW_WIDTH 640.0f

static const char* globalAttributeNames[] = {
"inVertex", "inUV", "inColor", "pointSize"
};
enum{VERTEX_ARRAY, TEXCOORD_ARRAY, COLOR_ARRAY, POINT_SIZE, numAttributes};

static const char* globalUniformNames[] = {
	  "mPMVMatrix","color", "tex1", "tex2", "tex3", "lightPos", "pixel"};
enum {eMatrix, eColor, eTex1, eTex2, eTex3, eLightPos, ePixel, numUniform};

typedef std::vector<std::pair<unsigned int, const char*> > attribute_pair_list;

static double diffclock(clock_t clock1, clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffs=(diffticks)/CLOCKS_PER_SEC;
	return diffs;
}

struct Shader
{
	GLuint id;
	GLuint vertShader, fragShader;
	GLuint* uniLoc;
	std::vector<unsigned int> attribs;
	Shader(){
		vertShader = -1; fragShader = -1; id = -1;
		uniLoc = new GLuint[numUniform];
	}
	void addUniforms(const std::vector<unsigned int>& uniforms){
		for(unsigned int i=0; i<uniforms.size(); i++){
			uniLoc[uniforms[i]] = glGetUniformLocation(id,globalUniformNames[uniforms[i]]);
		}
	}
	void addAttributes(const attribute_pair_list& attribs){
		for(unsigned int i=0; i<attribs.size(); i++){
			this->attribs.push_back(attribs[i].first);
		}
	}
	void bind(){
		glUseProgram(id);
		for(unsigned int i=0; i<attribs.size(); i++){
			glEnableVertexAttribArray(attribs[i]);
		}
	}
	void unbind(){
		glUseProgram(0);
		for(unsigned int i=0; i<attribs.size(); i++){
			glDisableVertexAttribArray(attribs[i]);
		}
	}
	~Shader(){
		delete [] uniLoc;
		GLenum ErrorCheckValue = glGetError();
		glUseProgram(0);
		if(vertShader!=-1){
			glDetachShader(id, vertShader);
			glDeleteShader(vertShader);
		}
		if(fragShader!=-1){
			glDetachShader(id, fragShader);
			glDeleteShader(fragShader);
		}
		if(id!=-1){
			glDeleteProgram(id);
		}
	}
};

struct shaderConstructionInfo{
	std::string vertShaderSource;
	std::string fragShaderSource;
	attribute_pair_list attributes;
	std::vector<unsigned int> uniforms;
	shaderConstructionInfo(int attributes[], unsigned int numAttributes, int uniforms[], unsigned int numUniforms){
		for(unsigned int i=0; i<numAttributes; i++){
			this->attributes.push_back(std::make_pair(attributes[i],globalAttributeNames[attributes[i]]));
		}
		for(unsigned int i=0; i<numUniforms; i++){
			this->uniforms.push_back(uniforms[i]);
		}
	}
};

namespace Renderer{

	void init(HINSTANCE hInstance, HINSTANCE hPrevInstance, TCHAR *lpCmdLine, int nCmdShow);
	void loadShader(const shaderConstructionInfo& info, std::string shaderName);
	Shader* getShader(std::string shaderName);
	void loadModel(std::string fileName, std::string modelName);
	//Model* getModel(std::string modelName);
	void loadTexture(std::string fileName, std::string texName);

	static glm::mat4 camera;
	static glm::mat4 perspective;

	void swapBuffers();
	void saveScreenshot(std::string filename, int w, int h);
	void destroy();

};

#endif