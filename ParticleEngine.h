#ifndef PARTICLEENGINE_H
#define PARTICLEENGINE_H

#include <glm/glm.hpp>
#include "Renderer.h"

struct Particle{
	glm::vec3 vel;
	glm::vec3 pos;

	Particle(){
		pos = glm::vec3(0,0,0);
		vel = glm::vec3(0,0,0);
	}
	Particle(float _x, float _y, float _z) : pos(_x,_y,_z){
		vel = glm::vec3(0,0,0);
	}
	Particle(const glm::vec3& position) : pos(position){
		vel = glm::vec3(0,0,0);
	}
};

struct FluidParticle : Particle{
	glm::vec3 sph_force;
	glm::vec4 col;
	float density;

	FluidParticle() : Particle(), sph_force(0,0,0), density(0.0){
	}

	FluidParticle(float _x, float _y, float _z) : Particle(_x,_y,_z), sph_force(0,0,0), density(0.0){
	}

	FluidParticle(const glm::vec3& position) : Particle(position) , sph_force(0,0,0), density(0.0){
	}

	void clearForces(){
		sph_force = glm::vec3(0,0,0);
	}

};

class ParticleEngine{
public:
	ParticleEngine(bool enableParticleCollisions = false);
	~ParticleEngine();
	virtual void addParticles(int numParticles);
	virtual void renderParticles(const glm::mat4& view, const glm::mat4& proj);

	void setShaderProgram(std::string shaderName);
	virtual void stepParticleSimulation();

protected:

	std::vector<FluidParticle> particles;

	//Collisions
	bool part2partCollisionEnabled;
	std::function<void(Particle* part, int nProcess)> part2partCollisionCallback;

	Shader* shader;

	int numParticles;
	glm::vec3 volMin, volMax;
	glm::vec3 gravity;

};

#endif