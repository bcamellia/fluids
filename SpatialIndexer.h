#ifndef SPATIAL_INDEX_H
#define SPATIAL_INDEX_H

#include <kdtree.h>
#include <octree.h>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <ppl.h>
#include "FluidSimControl.h"

enum{ OCTREE, KD_TREE, HASH_MAP, GRID};

namespace boost
{
	inline std::size_t hash_value(glm::ivec3 const& v)
    {
      std::size_t seed = 0;//( (unsigned int)(v.x*73856093) ^ (unsigned int)(v.y*19349663) ^ (unsigned int)(v.z*83492791) ) %  MAX_NUM_FLUID_PARTICLES;
	  hash_combine(seed,v.x);
	  hash_combine(seed,v.y);
	  hash_combine(seed,v.z);
	  return seed;
    }

	inline std::size_t hash_value(glm::ivec2 const& v)
    {
      std::size_t seed = 0;//( (unsigned int)(v.x*73856093) ^ (unsigned int)(v.y*19349663) ^ (unsigned int)(v.z*83492791) ) %  MAX_NUM_FLUID_PARTICLES;
	  hash_combine(seed,v.x);
	  hash_combine(seed,v.y);
	  return seed;
    }
}

typedef std::vector<int> NeighborList;
typedef boost::unordered_multimap<glm::ivec3, int, boost::hash<glm::ivec3> > hMap;

struct IndexContainer{
	IndexContainer(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int maxNumParticles, int type){
		this->cell_size = cell_size;
		this->volMin = volMin;
		this->type = type;
		this->maxNumParticles = maxNumParticles;
		max.x = abs(volMax.x-volMin.x)/cell_size;
		if(max.x == 0)
			max.x++;
		max.y = abs(volMax.y-volMin.y)/cell_size;
		if(max.y == 0)
			max.y++;
		max.z = abs(volMax.z-volMin.z)/cell_size;
		if(max.z == 0)
			max.z++;
	};
	~IndexContainer(){
	};

	int type;

	int maxNumParticles;
	float cell_size;
	glm::ivec3 max;
	glm::vec3 volMin;
	virtual void insert(glm::vec3 pos, int index) = 0;
	virtual void clear() = 0;
	virtual int find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index) = 0;

	//int getSpatialIndex(const glm::vec3& pos);
	glm::vec3 getMinCellPos(const glm::vec3& pos);
	glm::vec3 getActualCellPos(const glm::vec3& pos);
	glm::vec3 getMinCellPos(const glm::vec3& pos, float cell_size);
};

struct kdTree : IndexContainer{
	kdTree(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int maxNumParticles, int type);
	~kdTree();
	void* ptree;
	void insert(glm::vec3 pos, int index);
	void clear();
	int find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index);
};

struct octree : IndexContainer{
	octree(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int maxNumParticles, int type);
	~octree();
	Octree<std::vector<int> > otree;
	int maxLen;
	//std::vector<NeighborList> cells;
	void insert(glm::vec3 pos, int index);
	void clear();
	int find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index);
};

struct hashMap : IndexContainer{
	hashMap(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int maxNumParticles, int type);
	~hashMap();
	std::pair<hMap::iterator, hMap::iterator >** NCells;
	hMap map;
	void erase(glm::vec3 pos, int index);
	void insert(glm::vec3 pos, int index);
	void clear();
	int find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index);
};

struct grid : IndexContainer{ //intended for 2-Dimensional use!
	grid(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int maxNumParticles, int type);
	~grid();
	int** cells;
	int* nParticlesInCell;
	int gridSize;
	int getSpatialIndex(const glm::vec3& pos);
	void insert(glm::vec3 pos, int index);
	void clear();
	int find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index);
};
/*
struct grid3D : IndexContainer{ //intended for 3-Dimensional use!
	grid(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int type);
	~grid();
	int**** cells;
	int *** nParticlesInCell;
	hMap map;
	void insert(glm::vec3 pos, int index);
	void clear();
	int find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index);
};
*/

#endif