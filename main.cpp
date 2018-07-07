#include <iostream>
#include "Renderer.h"
#include "FluidSimulator.h"
#include "octree.h"
#include "Shaders.h"
#include <iomanip>
#include <glm/glm.hpp>

clock_t oldClock, nextClock;
int numFrames = 0, fps = 60.0;
float currentTime = 0;
int nFrames = 0;
std::vector<int> fpsData;
int WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, TCHAR * pCmdLine, int nCmdShow);
void calcFPS();
void writeFPS();


int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, TCHAR* pCmdLine, int nCmdShow){
	Renderer::init(hInstance, hPrevInstance, pCmdLine, nCmdShow);
	Renderer::camera = IDENTITY_MAT4;


	Renderer::camera = glm::translate(Renderer::camera, glm::vec3(.5,.8,5));
	Renderer::perspective = glm::perspective(45.0f, WINDOW_WIDTH/WINDOW_HEIGHT , 0.1f, 1000.0f);
	ShaderResource::initShaders();

	FluidSimulator fluid;
	fluid.setShaderProgram(PARTICLE_SHADER);

	char buffer[MAX_PATH];
	GetModuleFileName( NULL, buffer, MAX_PATH );
	std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
	pos = std::string( buffer ).find_last_of( "\\/" ,pos-1);
	std::string dir = std::string( buffer ).substr( 0, pos);
	char num[10];
	itoa(NUM_PARTICLES,num,10);
	dir = dir + "\\frames_" + std::string(num) + "\\";
	mkdir(dir.c_str());
	o = clock();
	fluid.addParticles(NUM_PARTICLES,GRID,NUM_PROCESS);
		while(1){
			glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
			fluid.stepParticleSimulation();
			glm::mat4 view = glm::inverse(Renderer::camera);
			fluid.renderParticles(view, Renderer::perspective);
			Renderer::swapBuffers();
			//calcFPS();
			if(fluid.mDone){
#ifdef SAVE_FPS_DATA
				writeFPS();
#endif
				break;
			}
			_RPT1(_CRT_WARN,"FPS: %i\n",fps);
#ifdef COMPILE_MOVIE_FRAMES
			std::ostringstream ss;
			ss << std::setw(4) << std::setfill('0') << nFrames;
			std::string result = ss.str();

			Renderer::saveScreenshot(dir + "img_" + result + ".tga", WINDOW_WIDTH, WINDOW_HEIGHT);
			nFrames++;
#endif
		}
		//fluid.reset();
		exit(1);
};

void calcFPS(){
	numFrames++;
	nextClock = clock();
	currentTime = diffclock(nextClock, oldClock);
	if(currentTime >= 1.0){
		fps = numFrames;
		numFrames = 0;
		currentTime = 0;
		oldClock = clock();
#ifdef SAVE_FPS_DATA
		fpsData.push_back(fps);
#endif
	}
}

void writeFPS(){
	char buffer[MAX_PATH];
	GetModuleFileName( NULL, buffer, MAX_PATH );
	std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
	pos = std::string( buffer ).find_last_of( "\\/" ,pos-1);
	std::string dir = std::string( buffer ).substr( 0, pos);
	char num[10];
	itoa(NUM_PARTICLES,num,10);
	dir = dir + "\\frames_" + std::string(num) + "\\";
	std::ofstream file1;
	file1.open(dir + "mFPS.txt");
	for(int i=0; i<fpsData.size(); i++){
		file1 << boost::lexical_cast<std::string>(fpsData[i]) << std::endl;
	}
	file1.close();
}