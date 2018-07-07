#include "octree.h"

using namespace glm;

template <class T>
octreeNode<T>::octreeNode(const vec3& pos, int depth, int maxDepth){
	if(depth<maxDepth){
		vec3 halfPos = pos/2;
		children[0] = octreeNode<T>(vec3(pos.x+halfPos.x,pos.y+halfPos.y,pos.z+halfPos.z),depth+1,maxDepth);
		children[1] = octreeNode<T>(vec3(pos.x-halfPos.x,pos.y+halfPos.y,pos.z+halfPos.z),depth+1,maxDepth);
		children[2] = octreeNode<T>(vec3(pos.x+halfPos.x,pos.y-halfPos.y,pos.z+halfPos.z),depth+1,maxDepth);
		children[3] = octreeNode<T>(vec3(pos.x+halfPos.x,pos.y+halfPos.y,pos.z-halfPos.z),depth+1,maxDepth);
		children[4] = octreeNode<T>(vec3(pos.x-halfPos.x,pos.y-halfPos.y,pos.z+halfPos.z),depth+1,maxDepth);
		children[5] = octreeNode<T>(vec3(pos.x+halfPos.x,pos.y-halfPos.y,pos.z-halfPos.z),depth+1,maxDepth);
		children[6] = octreeNode<T>(vec3(pos.x-halfPos.x,pos.y+halfPos.y,pos.z-halfPos.z),depth+1,maxDepth);
		children[7] = octreeNode<T>(vec3(pos.x-halfPos.x,pos.y-halfPos.y,pos.z-halfPos.z),depth+1,maxDepth);
		isLeaf = false;
	}else{
		isLeaf = true;
	}
}

template <class T>
octree<T>::octree(){
	
}

template <class T>
octreeNode<T>* octree<T>::find(T param){

}