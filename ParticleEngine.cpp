#include "ParticleEngine.h"

using namespace glm;

//#define round(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))

ParticleEngine::ParticleEngine(bool enableParticleCollisions){
	part2partCollisionEnabled = enableParticleCollisions;
}

ParticleEngine::~ParticleEngine(){
	particles.resize(0);
}

void ParticleEngine::renderParticles(const glm::mat4& view, const glm::mat4& proj){
	glm::mat4 world = glm::mat4(1,0,0,0,
								 0,1,0,0,
								 0,0,1,0,
								 0,0,0,1);

	shader->bind();
	
	glUniformMatrix4fvARB(shader->uniLoc[eMatrix], 1, GL_FALSE, glm::value_ptr(world));
    glUniform3fARB(shader->uniLoc[eLightPos], 0, 0, 0);
	glVertexAttribPointerARB(VERTEX_ARRAY, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), &particles[0].pos.x);
	glDrawArrays(GL_POINTS, 0, numParticles);
	
	shader->unbind();

}

void ParticleEngine::stepParticleSimulation(){
	
	//Concurrency::combinable<boost::unordered_multimap<vec3, Particle*> > new_part_lookup;
	//new_part_lookup.local() = boost::unordered_multimap<vec3, Particle*> (numParticles);

	//Concurrency::parallel_for(0,NUM_nProcess,[&,this](int nProcess){
	for(int i=0; i<numParticles; i++){
		Particle* currentPart = &particles[i];
		vec3 pos;
		pos = currentPart->pos;
	}

}

void ParticleEngine::addParticles(int numParticles){
	this->numParticles = numParticles;
}

void ParticleEngine::setShaderProgram(std::string shaderName){
	shader = Renderer::getShader(shaderName);
}