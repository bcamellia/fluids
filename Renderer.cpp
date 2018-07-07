#include "Renderer.h"
#include <TCHAR.h>

namespace Renderer{

	void showErrorMessage(std::string msg);
	Shader* compileShader(const shaderConstructionInfo& info);

	namespace{
		HWND hWnd;
		HDC hDC;
		HGLRC hRC;

		boost::unordered_map<std::string, Shader*> shaders;
	}

	void showErrorMessage(std::string msg){
		MessageBox(hWnd,msg.c_str(),"",MB_OK);
		exit(1);
	}

	void destroy(){	
		shaders.clear();
		wglMakeCurrent(hDC, 0); // Remove the rendering context from our device context  
		wglDeleteContext(hRC); // Delete our rendering context
	}

	bool createContext(){

		hDC = GetDC(hWnd);
		PIXELFORMATDESCRIPTOR pfd; // Create a new PIXELFORMATDESCRIPTOR (PFD)  
		memset(&pfd, 0, sizeof(PIXELFORMATDESCRIPTOR)); // Clear our  PFD  
		pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR); // Set the size of the PFD to the size of the class  
		pfd.dwFlags = PFD_DOUBLEBUFFER | PFD_SUPPORT_OPENGL | PFD_DRAW_TO_WINDOW; // Enable double buffering, opengl support and drawing to a window  
		pfd.iPixelType = PFD_TYPE_RGBA; // Set our application to use RGBA pixels  
		pfd.cColorBits = 32; // Give us 32 bits of color information (the higher, the more colors)  
		pfd.cDepthBits = 32; // Give us 32 bits of depth information (the higher, the more depth levels)  
		pfd.iLayerType = PFD_MAIN_PLANE; // Set the layer of the PFD  
		int nPixelFormat = ChoosePixelFormat(hDC, &pfd); // Check if our PFD is valid and get a pixel format back  
		if (nPixelFormat == 0){ 
			return false;  
		}

		bool bResult = SetPixelFormat(hDC, nPixelFormat, &pfd); // Try and set the pixel format based on our PFD  
		if (!bResult){ // If it fails  
			return false;  
		}
		HGLRC tempOpenGLContext = wglCreateContext(hDC); // Create an OpenGL 2.1 context for our device context  
		wglMakeCurrent(hDC, tempOpenGLContext); // Make the OpenGL 2.1 context current and active 
		return true;

	}

	LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
	{
		//switch (message)
		{
			return DefWindowProc(hWnd, message, wParam, lParam);
		}
	}

	void init(HINSTANCE hInstance, HINSTANCE hPrevInstance, TCHAR *lpCmdLine, int nCmdShow){

		hDC = GetDC(hWnd);
		PIXELFORMATDESCRIPTOR pfd; // Create a new PIXELFORMATDESCRIPTOR (PFD)  
		memset(&pfd, 0, sizeof(PIXELFORMATDESCRIPTOR)); // Clear our  PFD  
		pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR); // Set the size of the PFD to the size of the class  
		pfd.dwFlags = PFD_DOUBLEBUFFER | PFD_SUPPORT_OPENGL | PFD_DRAW_TO_WINDOW; // Enable double buffering, opengl support and drawing to a window  
		pfd.iPixelType = PFD_TYPE_RGBA; // Set our application to use RGBA pixels  
		pfd.cColorBits = 32; // Give us 32 bits of color information (the higher, the more colors)  
		pfd.cDepthBits = 32; // Give us 32 bits of depth information (the higher, the more depth levels)  
		pfd.iLayerType = PFD_MAIN_PLANE; // Set the layer of the PFD  
		int nPixelFormat = ChoosePixelFormat(hDC, &pfd); // Check if our PFD is valid and get a pixel format back  
		if (nPixelFormat == 0){ 

		}

		bool bResult = SetPixelFormat(hDC, nPixelFormat, &pfd); // Try and set the pixel format based on our PFD  
		if (!bResult){ // If it fails  
			//return false;  
		}
		HGLRC tempOpenGLContext = wglCreateContext(hDC); // Create an OpenGL 2.1 context for our device context  
		wglMakeCurrent(hDC, tempOpenGLContext); // Make the OpenGL 2.1 context current and active 

		WNDCLASS sWC;
		sWC.style = CS_HREDRAW | CS_VREDRAW;
		sWC.lpfnWndProc = WndProc;
		sWC.cbClsExtra = 0;
		sWC.cbWndExtra = 0;
		sWC.hInstance = hInstance;
		sWC.hIcon = 0;
		sWC.hCursor = 0;
		sWC.lpszMenuName = 0;
		sWC.hbrBackground = (HBRUSH) GetStockObject(WHITE_BRUSH);
		sWC.lpszClassName = "DEMO";

		ATOM registerClass = RegisterClass(&sWC);
		if (!registerClass)
		{
			MessageBox(0, ("Failed to register the window class"), ("Error"), MB_OK | MB_ICONEXCLAMATION);
		}

		RECT	sRect;
		SetRect(&sRect, 0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);
		AdjustWindowRectEx(&sRect, WS_CAPTION | WS_SYSMENU, false, 0);
		hWnd = CreateWindow( sWC.lpszClassName, _T(WINDOW_TITLE), WS_VISIBLE | WS_SYSMENU,
			0, 0, WINDOW_WIDTH, WINDOW_HEIGHT, NULL, NULL, hInstance, NULL);


		if(!createContext()){
			showErrorMessage("Couldn't initialize GLEW");
		}

		GLenum error = glewInit(); 

		if (error != GLEW_OK || !glewIsSupported("GL_VERSION_2_0")){
			showErrorMessage("Couldn't initialize GLEW");
		}

		glClearColor(.8f, 0.8f, 0.8f, 1.0f);
		glEnable(GL_DEPTH_TEST);
		glDepthMask(GL_TRUE);
		glClearDepth(1.f);
		glDepthFunc(GL_LEQUAL);
		glViewport(0,0,WINDOW_WIDTH,WINDOW_HEIGHT);

	}

	void loadModel(std::string fileName, std::string modelName){
	}

	/*Model* getModel(std::string modelName){
	auto it = &models.find(modelName)->second;
	if(it == nullptr){
	showErrorMessage("Internal issue locating model.."); 
	}
	return it;
	}*/

	Shader* getShader(std::string shaderName){
		int n = shaders.size();

		auto it = shaders.find(shaderName);
		if(it == shaders.end()){
			showErrorMessage("Internal issue locating shader.."); 
		}
		return it->second;
	}

	/*sf::Texture* getTexture(std::string texName){
	auto it = &textures.find(texName)->second;
	if(it == nullptr){
	showErrorMessage("Internal issue locating texture.."); 
	}
	return it;
	}*/

	void loadShader(const shaderConstructionInfo& info, std::string shaderName){
		Shader* shader = compileShader(info);
		shader->addAttributes(info.attributes);
		shader->addUniforms(info.uniforms);
		shaders.insert(std::make_pair(shaderName,shader));
	}

	/*void loadTexture(std::string fileName, std::string texName){
	sf::Texture tex;
	tex.loadFromFile(texName);
	textures.insert(std::make_pair(texName,tex));
	}*/

	Shader* compileShader(const shaderConstructionInfo& info){

		Shader* shader = new Shader();

		shader->fragShader = glCreateShader(GL_FRAGMENT_SHADER);

		const char* fragstr = info.fragShaderSource.c_str();
		const char* vertstr = info.vertShaderSource.c_str();
		glShaderSource(shader->fragShader, 1, (const char**)&fragstr, NULL);

		// Compile the source code
		glCompileShader(shader->fragShader);

		// Check if compilation succeeded
		GLint bShaderCompiled;
		glGetShaderiv(shader->fragShader, GL_COMPILE_STATUS, &bShaderCompiled);
		if (!bShaderCompiled)
		{
			// An error happened, first retrieve the length of the log message
			int i32InfoLogLength, i32CharsWritten;
			glGetShaderiv(shader->fragShader, GL_INFO_LOG_LENGTH, &i32InfoLogLength);

			// Allocate enough space for the message and retrieve it
			char* pszInfoLog = new char[i32InfoLogLength];
			glGetShaderInfoLog(shader->fragShader, i32InfoLogLength, &i32CharsWritten, pszInfoLog);

			/*
			Displays the message in a dialog box when the application quits
			using the shell PVRShellSet function with first parameter prefExitMessage.
			*/
			char* pszMsg = new char[i32InfoLogLength+256];
			strcpy(pszMsg, "Failed to compile fragment shader: ");
			strcat(pszMsg, pszInfoLog);

			showErrorMessage("Error creating fragment shader.");

			delete [] pszMsg;
			delete [] pszInfoLog;

		}

		// Loads the vertex shader in the same way
		shader->vertShader = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(shader->vertShader, 1, (const char**)&vertstr, NULL);
		glCompileShader(shader->vertShader);
		glGetShaderiv(shader->vertShader, GL_COMPILE_STATUS, &bShaderCompiled);
		if (!bShaderCompiled)
		{
			int i32InfoLogLength, i32CharsWritten;
			glGetShaderiv(shader->vertShader, GL_INFO_LOG_LENGTH, &i32InfoLogLength);
			char* pszInfoLog = new char[i32InfoLogLength];
			glGetShaderInfoLog(shader->vertShader, i32InfoLogLength, &i32CharsWritten, pszInfoLog);
			char* pszMsg = new char[i32InfoLogLength+256];
			strcpy(pszMsg, "Failed to compile vertex shader: ");
			strcat(pszMsg, pszInfoLog);

			showErrorMessage("Error creating vertex shader.");

			delete [] pszMsg;
			delete [] pszInfoLog;
			assert(false);
		}
		shader->id = glCreateProgram();

		// Attach the fragment and vertex shaders to it
		glAttachShader(shader->id, shader->fragShader);
		glAttachShader(shader->id, shader->vertShader);

		// Bind the custom attributes
		for(unsigned int i=0; i<info.attributes.size(); i++){
			std::pair<int,const char*> att = info.attributes[i];
			glBindAttribLocation(shader->id, att.first, att.second);
		}
		// Link the program
		glLinkProgram(shader->id);

		// Check if linking succeeded in the same way we checked for compilation success
		GLint bLinked;
		glGetProgramiv(shader->id, GL_LINK_STATUS, &bLinked);
		if (!bLinked)
		{
			int i32InfoLogLength, i32CharsWritten;
			glGetProgramiv(shader->id, GL_INFO_LOG_LENGTH, &i32InfoLogLength);
			char* pszInfoLog = new char[i32InfoLogLength];
			glGetProgramInfoLog(shader->id, i32InfoLogLength, &i32CharsWritten, pszInfoLog);
			char* pszMsg = new char[i32InfoLogLength+256];
			strcpy(pszMsg, "Failed to link program: ");
			strcat(pszMsg, pszInfoLog);

			delete [] pszMsg;
			delete [] pszInfoLog;
			assert(false);
		}

		return shader;

	}

	void swapBuffers(){
		SwapBuffers(hDC);
	}

	void saveScreenshot(std::string filename, int w, int h)
	{	
		//This prevents the images getting padded 
		// when the width multiplied by 3 is not a multiple of 4
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

		int nSize = w*h*3;
		// First let's create our buffer, 3 channels per Pixel
		char* dataBuffer = (char*)malloc(nSize*sizeof(char));

		// Let's fetch them from the backbuffer	
		// We request the pixels in GL_BGR format, thanks to Berzeger for the tip
		glReadPixels((GLint)0, (GLint)0,
			(GLint)w, (GLint)h,
			GL_BGR, GL_UNSIGNED_BYTE, dataBuffer);

		//Now the file creation
		FILE *filePtr = fopen(filename.c_str(), "wb");
			unsigned char TGAheader[12]={0,0,2,0,0,0,0,0,0,0,0,0};
			unsigned char header[6] = { w%256,w/256,
			h%256,h/256,
			24,0};
			// We write the headers
			fwrite(TGAheader,	sizeof(unsigned char),	12,	filePtr);
			fwrite(header,	sizeof(unsigned char),	6,	filePtr);
			// And finally our image data
			fwrite(dataBuffer,	sizeof(GLubyte),	nSize,	filePtr);
			fclose(filePtr);
			
	}
}