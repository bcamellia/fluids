#ifndef FLUIDSIMULATOR_H
#define FLUIDSIMULATOR_H

#include "ParticleEngine.h"
#include <time.h>
#include "SpatialIndexer.h"
#include <boost\lexical_cast.hpp>
#include <boost\assign.hpp>
#include <boost\iterator\counting_iterator.hpp>
#include <fstream>
#include <functional>
#include "FluidSimControl.h"

class FluidSimulator : public ParticleEngine{
public:
	FluidSimulator();
	~FluidSimulator();
	void initParameters(float pDist);
	void reset();
	void addParticles(int numParticles, int spatial_container_type, int numProcess);
	void renderParticles(const glm::mat4& view, const glm::mat4& proj);
	void stepParticleSimulation();
	bool mDone;

private:
	float mDT; //s
	float mTime; //s
	float mMass; //kg
	double mVolume; //m^3
	float mRestDensity; //kg / m^3
	float mInitialDensity; //kg / m^3
	float mGasConstant; //Pa(m^3) / K*mol
	float mViscosity; //kg / ms
	float mSurfTenThreshold2; //kg / m
	float mSurfTenCoefficient; //kg / m
	int numProcess;
	Concurrency::combinable<float> mKineticEnergy; //J

	float mRes;

	float CELL_SIZE; //length of a side in meters

	float poly6Kern, poly6GradientKern, poly6LaplacianKern, spikyGradientKern, viscLaplacianKern;
	float mSmoothRadius, mSmoothRadius2;

	clock_t prev;

	IndexContainer* map;

	int* numNeighbors;
	std::vector<glm::vec2> kinEnergy;
	int** mNeighbors;
	int** mNeighborTemp;
	float** nDists;

	glm::vec3 calcPoly6Grad(const glm::vec3& R, float h2r2);
	float calcPoly6Lap(float r2, float h2r2);
	glm::vec3 calcSpikyGrad(const glm::vec3& R, float r);
	float calcViscLap(float r);

	void computeDensity(FluidParticle* fluidPart, int partIndex, int nProcess);
	void computeForces(FluidParticle* fluidPart, int nProcess);
	void advanceParticles(FluidParticle* fluidpart, int nProcess);

	void initFluidVisualization();
	void buildFluidMeshMarchingCubes();

	std::vector<glm::vec3>* mCubeVertices;
	std::vector<unsigned short>* mCubeIndices;
	Concurrency::combinable<int>** pointField;
	int maxPointFieldX;
	int maxPointFieldY;
	Concurrency::combinable<std::vector<glm::ivec2> > cellsToCheck;
	boost::unordered_set<glm::ivec2,boost::hash<glm::ivec2> > accumCells;
	std::vector<glm::vec3> fluidMesh;
	std::vector<glm::vec4> fluidMeshCol;
	std::vector<unsigned short> fluidMeshIndices;
	float meshResolution;

};

#endif