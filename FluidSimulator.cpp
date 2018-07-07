#include "FluidSimulator.h"

using namespace glm;

float interval = 100/1000.0;
float incr = interval;
clock_t timer1, timer2;
float avgDens = 0;
float avgForce = 0;
float avgUpdate = 0;
float avgNeighbor = 0;
int frame = 0;

float inline __declspec (naked) __fastcall fsqrt(float n)
{
	_asm fld qword ptr [esp+4]
    _asm fsqrt
	_asm ret 8
}

FluidSimulator::FluidSimulator() : ParticleEngine(true) // collisions must be enabled
{
}

void FluidSimulator::initParameters(float pDist){

	mDT = .002;

	part2partCollisionEnabled = true;

	mGasConstant = 1.0; // J/Kg-K for water - http://www.engineeringtoolbox.com/individual-universal-gas-constant-d_588.html (compromise to maintain larger time step) 
	mSurfTenCoefficient = 0.0712; //water at 30 C - http://www.engineeringtoolbox.com/water-surface-tension-d_597.html
	mViscosity = .000798; // water at 30 C - http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html 
	mKineticEnergy.clear();
	mDone = false;
	mRes = 1.0;
	mRestDensity = 1000.0;
	int avgNumP = 20;

	mSmoothRadius = pow((3*avgNumP*mVolume)/(4*M_PI*numParticles), 1/3.0); //copenhagen article, setup volume(1 cm^3)
	CELL_SIZE = 2.0 * mSmoothRadius;
	mMass = mVolume * mRestDensity / numParticles;

	mTime = 0;

	mSmoothRadius2 = mSmoothRadius*mSmoothRadius;
	poly6GradientKern = -945.0f / (32.0f * M_PI * pow(mSmoothRadius, 9) );	
	poly6LaplacianKern = 945.0f / (8.0f * M_PI * pow(mSmoothRadius, 9) );	
	poly6Kern = 315.0f / (64.0f * M_PI * pow(mSmoothRadius, 9) );
	spikyGradientKern = -45.0f / (M_PI * pow(mSmoothRadius, 6) );
	viscLaplacianKern = 45.0f / (M_PI * pow(mSmoothRadius, 5) );

	mInitialDensity = 2 * mMass * (pow(mSmoothRadius2,3)) * poly6Kern; //density of initial particle (mass * W), no distance to self

	mSurfTenThreshold2 = mRestDensity / avgNumP;
	gravity = vec3(0,-9.8,0);

#ifdef SAVE_KIN_DATA
	kinEnergy.resize((int)(MAX_TIME/mDT)+2);
#endif

}

glm::vec3 FluidSimulator::calcPoly6Grad(const glm::vec3& R, float h2r2){
	float m = poly6GradientKern * (h2r2*h2r2);
	return R * m;
}

float FluidSimulator::calcPoly6Lap(float r2, float h2r2){
	float b = (r2 - .75*h2r2);
	return poly6LaplacianKern * h2r2 * b;
}

glm::vec3 FluidSimulator::calcSpikyGrad(const glm::vec3& R, float r){
	float a = mSmoothRadius - r;
	float m = (spikyGradientKern/r) * (a*a);
	return R*m;
}

float FluidSimulator::calcViscLap(float r){
	return viscLaplacianKern * (1.0 - (r / mSmoothRadius) );
}

void FluidSimulator::reset(){

	kinEnergy.resize(0);
	kinEnergy.resize(0);
	frame = 0;

	for(int i=0; i<numParticles; i++){
		delete [] nDists[i];
		delete [] mNeighbors[i];
	}
	for(int i=0; i<numProcess; i++){
		delete [] mNeighborTemp[i];
	}
	delete [] mNeighborTemp;
	delete [] nDists;
	delete [] mNeighbors;
	delete [] numNeighbors;
	particles.resize(0);
	delete map;
}

FluidSimulator::~FluidSimulator(){
	reset();
}

void FluidSimulator::stepParticleSimulation(){

	/*if(GetAsyncKeyState(VK_UP)){
		gravity = vec3(0,9.8,0);
	}
	if(GetAsyncKeyState(VK_LEFT)){
		gravity = vec3(-9.8,0,0);
	}
	if(GetAsyncKeyState(VK_DOWN)){
		gravity = vec3(0,-9.8,0);
	}
	if(GetAsyncKeyState(VK_RIGHT)){
		gravity = vec3(9.8,0,0);
	}*/

	timer1 = clock();
	switch(map->type){
	case GRID: //multiple processes possible in grid mode 
		Concurrency::parallel_for(0, numProcess, [this](int index){
			int begin = numParticles/numProcess * index;
			int end = begin + (numParticles/numProcess);
			for(int i=begin; i < end; i++){
				FluidParticle* currentPart = &particles[i];
				computeDensity(currentPart, i, index);
			}
		});
		break;
	default: //linear
		for(int i=0; i<numParticles; i++){
			FluidParticle* currentPart = &particles[i];
			computeDensity(currentPart, i, 0);	
		}
		break;
	}
	timer2 = clock();
	avgDens += diffclock(timer2, timer1);


	timer1 = clock();
	Concurrency::parallel_for(0, numProcess, [this](int index){
		int begin = numParticles/numProcess * index;
		int end = begin + (numParticles/numProcess);
		for(int i=begin; i < end; i++){
			FluidParticle* currentPart = &particles[i];
			vec3 pos = currentPart->pos;
			computeForces(currentPart,i);
		}
	});
	timer2 = clock();
	avgForce += diffclock(timer2, timer1);

	timer1 = clock();
	Concurrency::parallel_for(0, numProcess, [this](int index){
		int begin = numParticles/numProcess * index;
		int end = begin + (numParticles/numProcess);
		for(int i=begin; i < end; i++){
		FluidParticle* currentPart = &particles[i];
			advanceParticles(currentPart,i);
		}
	});
	timer2 = clock();
	avgUpdate += diffclock(timer2, timer1);

#ifdef BUILD_MESH
	buildFluidMeshMarchingCubes();
#endif

	timer1 = clock();
	map->clear();

	for(int i=0; i<numParticles; i++){
		FluidParticle* currentPart = &particles[i];
		vec3 pos = currentPart->pos;
		map->insert(pos, i);
	}
	timer2 = clock();
	avgNeighbor += diffclock(timer2, timer1);

	_RPT1(_CRT_WARN,"Force Speed: %f\n",avgForce);
	_RPT1(_CRT_WARN,"Update Speed: %f\n",avgUpdate);
	_RPT1(_CRT_WARN,"Density Speed: %f\n",avgDens);
	_RPT1(_CRT_WARN,"Grid Reconstr. Speed: %f\n",avgNeighbor);

#ifdef SAVE_KIN_DATA

	float kinTotal = 0.5 * mMass * mKineticEnergy.combine([](float l, float r){
		return l+r;
	});

	mKineticEnergy.clear();

	kinEnergy[frame] = (vec2(mTime, kinTotal));
	
	frame++;


	if(mTime>=MAX_TIME-mDT){
		char buffer[MAX_PATH];
		GetModuleFileName( NULL, buffer, MAX_PATH );
		std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
		pos = std::string( buffer ).find_last_of( "\\/" ,pos-1);
		std::string dir = std::string( buffer ).substr( 0, pos);
		char num[10];
		itoa(numParticles,num,10);
		dir = dir + "\\frames_" + std::string(num) + "\\";
		mkdir(dir.c_str());
		std::ofstream file1, file2, file3;
		file1.open(dir + "mTime.txt");
		file2.open(dir + "mKinEnergy.txt");
		file3.open(dir + "mSpeed.txt");
		for(int i=0; i<kinEnergy.size(); i++){
			file1 << boost::lexical_cast<std::string>(kinEnergy[i].x) << std::endl;
			file2 << boost::lexical_cast<std::string>(kinEnergy[i].y) << std::endl;
		}
		file3 << "Force Speed: " << boost::lexical_cast<std::string>(avgForce) << std::endl;
		file3 << "Update Speed: " << boost::lexical_cast<std::string>(avgUpdate) << std::endl;
		file3 << "Density Speed: " << boost::lexical_cast<std::string>(avgDens) << std::endl;
		file3 << "Grid Reconstr. Speed: " << boost::lexical_cast<std::string>(avgNeighbor) << std::endl;

		file1.close();
		file2.close();
		file3.close();

		mDone = true;
	}
#endif

	mTime += mDT;

}


void FluidSimulator::computeDensity(FluidParticle* fluidPart, int partIndex, int nProcess){
	int n=0;
	float r, c, r2;
	float sum = 0.0;
	int pos = 0;
	float maxDist = 0.0;

	int count = map->find(fluidPart->pos, mSmoothRadius, mNeighborTemp[nProcess], MAX_NEIGHBORS, partIndex);

	numNeighbors[partIndex] = 0;
	FluidParticle* fluidNeighbor = 0;

	//particle - to - particle
	for(int i=0; i<count; i++){
		pos = mNeighborTemp[nProcess][i];
		if(pos == partIndex){
			continue;
		}
		fluidNeighbor = &particles[pos];

		r2 = glm::length2((fluidPart->pos) - (fluidNeighbor->pos));

		if(r2 < mSmoothRadius2){ //used to avoid sqrt in distance calculation
			r = sqrt(r2);
			if(r > maxDist){
				maxDist = r;
			}
			c = mSmoothRadius2 - r2;
			sum += (c*c*c);
			nDists[partIndex][n] = r;
			mNeighbors[partIndex][n] = pos;
			n++;
			if(n == MAX_NEIGHBORS){
				break;
			}
		}
	}
	
	numNeighbors[partIndex] = n;

	sum *= mMass * poly6Kern;
	fluidPart->density += sum;

}

void FluidSimulator::computeForces(FluidParticle* fluidPart, int nProcess){
	FluidParticle* fluidNeighbor;
	vec3 force(0,0,0);
	vec3 viscosity_force(0,0,0), pressure_force(0,0,0);
	float dp1, dp2;
	float pp1, pp2;
	float r, r2, h2r2;
	float col_lap = 0;
	vec3 col_grad(0,0,0);
	vec3 rVec(0,0,0);

	//particle to particle

	dp1 = fluidPart->density;
	pp1 = mGasConstant * (dp1 - mRestDensity);

	for(int i=0; i<numNeighbors[nProcess]; i++){
		fluidNeighbor = &particles[mNeighbors[nProcess][i]];
		rVec = fluidPart->pos - fluidNeighbor->pos;
		dp2 = fluidNeighbor->density;

		pp2 = mGasConstant * (dp2 - mRestDensity);
		r = nDists[nProcess][i];
		r2 = r*r;
		h2r2 = mSmoothRadius2 - r2;

		pressure_force += (pp1+pp2) / dp2 * calcSpikyGrad(rVec, r);
		viscosity_force += (fluidNeighbor->vel - fluidPart->vel) / dp2 * calcViscLap(r);

		col_grad += ((mMass * calcPoly6Grad(rVec,h2r2)) / dp2);
		col_lap += ((mMass * calcPoly6Lap(r2,h2r2)) / dp2);
	}

	pressure_force *= -0.5f;
	viscosity_force *= mViscosity;

	force += mMass * (pressure_force + viscosity_force);

	float gradLen2 = glm::length2(col_grad);
	if (gradLen2 >= mSurfTenThreshold2){
		force += -mSurfTenCoefficient * col_lap * col_grad / sqrt(gradLen2);	
	}

	fluidPart->sph_force = force;

}

void FluidSimulator::advanceParticles(FluidParticle* fluidPart, int nProcess){

	vec3 new_vel = fluidPart->vel + (fluidPart->sph_force/fluidPart->density + gravity) * mDT;
	fluidPart->pos += 0.5f * (new_vel + fluidPart->vel) * mDT;
	fluidPart->vel = new_vel;

#ifdef SAVE_KIN_DATA
	mKineticEnergy.local() += glm::length2(fluidPart->vel);
#endif

	//particle to wall

	//+-x
	float mag, dist, EPSILON = .000001f; //prevent division by zero if particle not moving
	float velMag = glm::length(fluidPart->vel);
	float d1 = (fluidPart->pos.x - volMin.x);
	float d2 = (volMax.x - fluidPart->pos.x);
	if(d1<d2){
		dist = d1;
		mag = 1.0;
	}else{
		dist = d2;
		mag = -1.0;
	}
	if(dist < 0){
		fluidPart->vel -= (1+mRes*(abs(dist)/(mDT*velMag+EPSILON))) * glm::vec3(mag,0,0) * glm::dot(fluidPart->vel, glm::vec3(mag,0,0));
	}

	//+-y
	d1 = (fluidPart->pos.y - volMin.y);
	d2 = (volMax.y - fluidPart->pos.y);
	if(d1<d2){
		dist = d1;
		mag = 1.0;
	}else{
		dist = d2;
		mag = -1.0;
	}
	if(dist < 0){
		fluidPart->vel -= (1+mRes*(abs(dist)/(mDT*velMag+EPSILON))) * glm::vec3(0,mag,0) * glm::dot(fluidPart->vel, glm::vec3(0,mag,0));
	}

#ifdef IS_3D
	//+-z
	d1 = (fluidPart->pos.z - volMin.z);
	d2 = (volMax.z - fluidPart->pos.z);
	if(d1<d2){
		dist = d1;
		mag = 1.0;
	}else{
		dist = d2;
		mag = -1.0;
	}
	if(dist < 0){
		glm::vec3 p = fluidPart->vel;
		fluidPart->vel -= (1+mRes*(abs(dist)/(mDT*velMag+EPSILON))) * glm::vec3(0,0,mag) * glm::dot(fluidPart->vel, glm::vec3(0,0,mag));
	}
#endif

#ifdef BUILD_MESH
	vec3 cell_pos = map->getMinCellPos(fluidPart->pos, meshResolution);
	int cx = (int)cell_pos.x;
	int cy = (int)cell_pos.y;
	if((cx < maxPointFieldX) && (cy < maxPointFieldY)){
		pointField[cx][cy].local()++;
	}
#endif

	fluidPart->clearForces();
	fluidPart->density = mInitialDensity;
}

void FluidSimulator::addParticles(int numParticles, int spatial_container_type, int numProcess){

	this->numProcess = numProcess;

	__super::addParticles(numParticles);

	mVolume = .01;

	initParameters(0);

#ifdef IS_3D
	int length = pow(numParticles,1/3.0)+1;
	incr = pow(mVolume,(float)1.0/3.0) / length;
#else
	int length = sqrt((float)numParticles);
	incr = mSmoothRadius*0.5;//pow(mVolume,(float)1.0/2.0) / length;
#endif

	

	assert(CELL_SIZE>0);

	assert(numParticles % numProcess == 0);

	nDists = new float*[numParticles];
	mNeighbors = new int*[numParticles];
	mNeighborTemp = new int*[numProcess];
	numNeighbors = new int[numParticles];

	for(int i=0; i<numProcess; i++){
		mNeighborTemp[i] = new int[numParticles];
	}

	for(int i=0; i<numParticles; i++){
		nDists[i] = new float[MAX_NEIGHBORS+3];
		mNeighbors[i] = new int[MAX_NEIGHBORS];
	}

#ifdef IS_3D
	volMin = vec3(0, 0, 0);
	volMax = vec3(1, 1, 1);
#else
	volMin = vec3(0, 0, 0);
	volMax = vec3(1, 1, 0);
#endif

	switch(spatial_container_type){
	case OCTREE:
		map = new octree(CELL_SIZE, volMin, volMax, numParticles, OCTREE);
		break;
	case KD_TREE:
		map = new kdTree(CELL_SIZE, volMin, volMax, numParticles, KD_TREE);
		break;
	case HASH_MAP:
		map = new hashMap(CELL_SIZE, volMin, volMax, numParticles, HASH_MAP);
		break;
	case GRID:
		map = new grid(CELL_SIZE, volMin, volMax, numParticles, GRID);
		break;
	default:
		map = new hashMap(CELL_SIZE, volMin, volMax, numParticles, HASH_MAP);
		break;
	}
	
	int index=0;

	float mOffset = 0.0;

	for(int i=0; i<length; i++){
		for(int j=0; j<length; j++){
#ifdef IS_3D
			for(int k=0; k<length; k++){
			FluidParticle p = FluidParticle(volMin.x+mOffset+(float)i*incr,volMin.y+mOffset+(float)j*incr,volMin.z+mOffset+(float)k*incr);
			p.vel = vec3(0,0,0);
			p.density = mInitialDensity;
			p.col = glm::vec4(0.0,0,1.0,1.0);
			particles.push_back(p);
			map->insert(p.pos, index);
			index++;
			}
#else
			FluidParticle p = FluidParticle(volMin.x+(float)i*incr,volMin.y+(float)j*incr,0);
			//p.pointSize = 1.0;
			p.vel = vec3(0,0,0);
			p.density = mInitialDensity;
			p.col = glm::vec4(0.0,0,1.0,1.0);
			particles.push_back(p);
			map->insert(p.pos, index);
			index++;
#endif
		}
	}

#ifdef BUILD_MESH
	meshResolution = 1.0;
	initFluidVisualization();
#endif

	prev = clock();
}

void FluidSimulator::initFluidVisualization(){

http://en.wikipedia.org/wiki/File:Marchingsquaresalgorithm.png

	/*
	2	3

	0	1
	*/
	using namespace boost::assign;

	maxPointFieldX = (map->max.x * CELL_SIZE) / meshResolution;
	maxPointFieldY = (map->max.y * CELL_SIZE) / meshResolution;
	pointField = new Concurrency::combinable<int>*[maxPointFieldX];
	for(int i=0; i<maxPointFieldX; i++){
		pointField[i] = new Concurrency::combinable<int>[maxPointFieldY];
	}

	float cell_size = meshResolution;

	mCubeVertices = new std::vector<glm::vec3>[16];
	mCubeIndices = new std::vector<unsigned short>[16];
	int index;
	std::vector<glm::vec3> tempV;
	std::vector<unsigned short> tempI;

	//Case 0
	index = (0)|(0<<1)|(0<<2)|(0<<3);
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 1
	index = (1)|(0<<1)|(0<<2)|(0<<3);
	tempV += vec3(0,0,0), vec3(0,cell_size/2,0), vec3(cell_size/2,0,0);
	tempI += 0,1,2;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 2
	index = (0)|(1<<1)|(0<<2)|(0<<3);
	tempV += vec3(cell_size,0,0), vec3(cell_size,cell_size/2,0), vec3(cell_size/2,0,0);
	tempI += 0,1,2;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 3
	index = (1)|(1<<1)|(0<<2)|(0<<3);
	tempV += vec3(0,0,0), vec3(0,cell_size/2,0), vec3(cell_size,0,0), vec3(cell_size,cell_size/2,0);
	tempI += 0,1,2,1,2,3;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 4
	index = (0)|(0<<1)|(0<<2)|(1<<3);
	tempV += vec3(cell_size,cell_size,0), vec3(cell_size/2,cell_size,0), vec3(cell_size,cell_size/2,0);
	tempI += 0,1,2;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 5
	index = (1)|(0<<1)|(0<<2)|(1<<3);
	tempV += vec3(0,0,0), vec3(0,cell_size/2,0), vec3(cell_size/2,0,0), vec3(cell_size/2,cell_size,0), vec3(cell_size,cell_size/2,0), vec3(cell_size,cell_size,0);
	tempI += 0,1,2,1,2,3,2,3,4,3,4,5;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 6
	index = (0)|(1<<1)|(0<<2)|(1<<3);
	tempV += vec3(cell_size,0,0), vec3(cell_size,cell_size,0), vec3(cell_size/2,0,0), vec3(cell_size/2,cell_size,0);
	tempI += 0,1,2,1,2,3;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 7
	index = (1)|(1<<1)|(0<<2)|(1<<3);
	tempV += vec3(0,0,0), vec3(0,cell_size/2,0), vec3(cell_size,0,0), vec3(cell_size/2,cell_size,0), vec3(cell_size,cell_size,0);
	tempI += 0,1,2,1,2,3,2,3,4;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 8
	index = (0)|(0<<1)|(1<<2)|(0<<3);
	tempV += vec3(0,cell_size/2,0), vec3(0,cell_size,0), vec3(cell_size/2,cell_size,0);
	tempI += 0,1,2;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 9
	index = (1)|(0<<1)|(1<<2)|(0<<3);
	tempV += vec3(0,0,0), vec3(cell_size/2,0,0), vec3(0,cell_size,0), vec3(cell_size/2,cell_size,0);
	tempI += 0,1,2,1,2,3;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 10
	index = (0)|(1<<1)|(1<<2)|(0<<3);
	tempV += vec3(cell_size,0,0), vec3(cell_size/2,0,0), vec3(cell_size,cell_size/2,0), vec3(0,cell_size/2,0), vec3(cell_size/2,cell_size,0), vec3(0,cell_size,0);
	tempI += 0,1,2,1,2,3,2,3,4,3,4,5;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 11
	index = (1)|(1<<1)|(1<<2)|(0<<3);
	tempV += vec3(0,cell_size,0), vec3(cell_size/2,cell_size,0), vec3(0,0,0), vec3(cell_size,cell_size/2,0), vec3(cell_size,0,0);
	tempI += 0,1,2,1,2,3,2,3,4;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 12
	index = (0)|(0<<1)|(1<<2)|(1<<3);
	tempV += vec3(0,cell_size,0), vec3(0,cell_size/2,0), vec3(cell_size,cell_size,0), vec3(cell_size,cell_size/2,0);
	tempI += 0,1,2,1,2,3;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 13
	index = (1)|(0<<1)|(1<<2)|(1<<3);
	tempV += vec3(0,0,0), vec3(cell_size/2,0,0), vec3(0,cell_size,0), vec3(cell_size,cell_size/2,0), vec3(cell_size,cell_size,0);
	tempI += 0,1,2,1,2,3,2,3,4;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 14
	index = (0)|(1<<1)|(1<<2)|(1<<3);
	tempV += vec3(0,cell_size,0), vec3(0,cell_size/2,0), vec3(cell_size,cell_size,0), vec3(cell_size/2,0,0), vec3(cell_size,0,0);
	tempI += 0,1,2,1,2,3,2,3,4;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();

	//Case 15
	index = (1)|(1<<1)|(1<<2)|(1<<3);
	tempV += vec3(0,0,0), vec3(0,cell_size,0), vec3(cell_size,0,0), vec3(cell_size,cell_size,0);
	tempI += 0,1,2,1,2,3;
	mCubeVertices[index] = tempV;
	mCubeIndices[index] = tempI;
	tempV.clear(); tempI.clear();
}

void FluidSimulator::buildFluidMeshMarchingCubes(){

	fluidMesh.clear();
	fluidMeshIndices.clear();
	fluidMeshCol.clear();

		for(int i=0; i<maxPointFieldX; i++){
			for(int j=0; j<maxPointFieldY; j++){
				int c1, c2, c3, c4;
				c1 = pointField[i][j].combine(std::plus<int>()) > 0 ? 1 : 0;
				//.combine(std::plus<int>())

				if((i+1 < maxPointFieldX) && (j+1 < maxPointFieldY)){
					c2 = pointField[i+1][j].combine(std::plus<int>()) > 0 ? 1 : 0;
					c3 = pointField[i][j+1].combine(std::plus<int>()) > 0 ? 1 : 0;
					c4 = pointField[i+1][j+1].combine(std::plus<int>()) > 0 ? 1 : 0;
				}else if(j+1 < maxPointFieldY){
					c3 = pointField[i][j+1].combine(std::plus<int>()) > 0 ? 1 : 0;
					c2 = (c1 > 0) ? 1 : 0;
					c4 = ((c1 == 1) && (c3 == 1)) ? 1 : 0;
				}else if(i+1 < maxPointFieldX){
					c2 = pointField[i+1][j].combine(std::plus<int>()) > 0 ? 1 : 0;
					c3 = (c1 > 0) ? 1 : 0;
					c4 = ((c1 == 1) && (c2 == 1)) ? 1 : 0;
				}else{
					c2 = (c1 > 0) ? 1 : 0;
					c3 = c2;
					c4 = c2;
				}

				int index = (c1)|(c2<<1)|(c3<<2)|(c4<<3);
				if(c1>0){
					pointField[i][j].clear();
				}
				if(index>0){
					auto cubeStateVerts = mCubeVertices[index];
					auto cubeStateIndices = mCubeIndices[index];
					vec3 pos = (vec3(i,j,0)*meshResolution) + volMin;
					int nextTriStart;
					if(fluidMeshIndices.size()==0){
						nextTriStart = 0;
					}else{
						nextTriStart = fluidMeshIndices.back()+1;
					}
					std::for_each(cubeStateVerts.begin(), cubeStateVerts.end(), [&](vec3& v) { v += pos;});
					for(int i=0; i<cubeStateIndices.size(); i++){
						fluidMeshCol.push_back(vec4(0.0,0.0,1.0,1.0));
					}
					std::for_each(cubeStateIndices.begin(), cubeStateIndices.end(), [&](unsigned short& v) { v += nextTriStart;});
					fluidMesh.insert(fluidMesh.end(), cubeStateVerts.begin(), cubeStateVerts.end());
					fluidMeshIndices.insert(fluidMeshIndices.end(), cubeStateIndices.begin(), cubeStateIndices.end());
				}
			}
		}
}

void FluidSimulator::renderParticles(const glm::mat4& view, const glm::mat4& proj){

	shader->bind();

	glm::mat4 rot = IDENTITY_MAT4;

#ifdef IS_3D
	rot = rotate(rot, 15.0f, vec3(1,1,0));
#endif

	glm::mat4 world = view * rot;
	world = proj * world;

	glUniformMatrix4fv(shader->uniLoc[eMatrix], 1, GL_FALSE, glm::value_ptr(world));

	float minX = volMin.x;
	float maxX = volMax.x;
	float maxY = volMax.y;
	float minY = volMin.y;
	float minZ = volMin.z;
	float maxZ = volMax.z;

	float lines[] = {
		minX, minY, minZ,
		maxX, minY, minZ,
		minX, maxY, minZ,
		maxX, maxY, minZ,
		minX, minY, minZ,
		minX, maxY, minZ,
		maxX, minY, minZ,
		maxX, maxY, minZ,

		minX, minY, maxZ,
		maxX, minY, maxZ,
		minX, maxY, maxZ,
		maxX, maxY, maxZ,
		minX, minY, maxZ,
		minX, maxY, maxZ,
		maxX, minY, maxZ,
		maxX, maxY, maxZ,

		minX, minY, minZ,
		minX, minY, maxZ,
		minX, maxY, minZ,
		minX, maxY, maxZ,
		minX, minY, minZ,
		minX, maxY, minZ,
		minX, minY, maxZ,
		minX, maxY, maxZ,

		maxX, minY, minZ,
		maxX, minY, maxZ,
		maxX, maxY, minZ,
		maxX, maxY, maxZ,
		maxX, minY, minZ,
		maxX, maxY, minZ,
		maxX, minY, maxZ,
		maxX, maxY, maxZ,
	};

#ifdef BUILD_MESH
	glVertexAttribPointer(VERTEX_ARRAY, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), &fluidMesh[0].x);
	glVertexAttribPointer(COLOR_ARRAY, 4, GL_FLOAT, GL_FALSE, sizeof(vec4), &fluidMeshCol[0].x);
	glVertexAttribPointer(POINT_SIZE, 1, GL_FLOAT, GL_FALSE, sizeof(FluidParticle), &particles[0].pointSize);
	glDrawElements(GL_TRIANGLES, fluidMeshIndices.size(), GL_UNSIGNED_SHORT, &fluidMeshIndices[0]);
#else
	glVertexAttribPointer(VERTEX_ARRAY, 3, GL_FLOAT, GL_FALSE, sizeof(FluidParticle), &particles[0].pos.x);
	glVertexAttribPointer(COLOR_ARRAY, 4, GL_FLOAT, GL_FALSE, sizeof(FluidParticle), &particles[0].col.x);
	glDrawArrays(GL_POINTS, 0, numParticles);
#endif
	/*float mLine[4] = {0,volMax.y+1, mTime, volMax.y+1};

	if(kinEnergy.size()>2){
	glVertexAttribPointer(VERTEX_ARRAY, 2, GL_FLOAT, GL_FALSE, sizeof(vec2), &kinEnergy[0].x);
	glVertexAttribPointer(COLOR_ARRAY, 4, GL_FLOAT, GL_FALSE, sizeof(vec4), &kinEnergyCol[0].x);
	glDrawArrays(GL_LINE_STRIP, 0, kinEnergy.size());
	glVertexAttribPointer(VERTEX_ARRAY, 2, GL_FLOAT, GL_FALSE, 0, &mLine[0]);
	glDrawArrays(GL_LINE_STRIP, 0, 2);
	}
	*/
	glVertexAttribPointer(VERTEX_ARRAY, 3, GL_FLOAT, GL_FALSE, 0, &lines[0]);
	glDrawArrays(GL_LINES, 0, 32);


	shader->unbind();
}