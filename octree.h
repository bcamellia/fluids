#ifndef OCTREE_H
#define OCTREE_H

#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

template <class T>
struct octreeNode{
private:
	std::vector<T> elements;
	bool isLeaf;
	glm::vec3 pos;
	octreeNode<T>(const glm::vec3& pos, int depth, int maxDepth);
	octreeNode<T>* children[8];
public:
	octreeNode<T>* getChild();
};

template <class T>
class octree{
	octree<T>();
	octreeNode<T> tree;
	octreeNode<T>* find(T param);
	void insert();
};

#endif