#include "SpatialIndexer.h"
#include <time.h>

float avgDen = 0;

glm::vec3 IndexContainer::getActualCellPos(const glm::vec3& pos){
	glm::vec3 id = glm::vec3((pos.x-volMin.x)/cell_size, (pos.y-volMin.y)/cell_size, (pos.z-volMin.z)/cell_size);
	return id;
}

glm::vec3 IndexContainer::getMinCellPos(const glm::vec3& pos, float cell_size){
	glm::ivec3 id = glm::ivec3((pos.x-volMin.x)/cell_size, (pos.y-volMin.y)/cell_size, (pos.z-volMin.z)/cell_size);
	return glm::vec3(id.x,id.y,id.z);
}

glm::vec3 IndexContainer::getMinCellPos(const glm::vec3& pos){
	glm::ivec3 id = glm::ivec3((pos.x-volMin.x)/cell_size, (pos.y-volMin.y)/cell_size, (pos.z-volMin.z)/cell_size);
	return glm::vec3(id.x,id.y,id.z);
}

kdTree::kdTree(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int maxNumParticles, int type) : IndexContainer(cell_size, volMin, volMax, maxNumParticles, type){
	ptree = kd_create(3); //3-dimensions
}

kdTree::~kdTree(){
	kd_free(static_cast<kdtree*>(ptree));
}

void kdTree::clear(){
	kd_clear(static_cast<kdtree*>(ptree));
}

void kdTree::insert(glm::vec3 pos, int index){
	intptr_t j = index;
	kd_insert3(static_cast<kdtree*>(ptree), pos.x, pos.y, pos.z, (void*)j);
}

int kdTree::find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index){
	int count = 0;
	kdres* presults = kd_nearest_range3f(static_cast<kdtree*>(ptree), pos.x, pos.y, pos.z, radius_sq);
	while( !kd_res_end( presults ) ) {
		intptr_t index = (int)kd_res_item_data( presults );
		nArray[count] = (int)index;
		count++;
		if(count == rowLen){break;}
		kd_res_next( presults );
	}
	kd_res_free( presults );
	return count;
}

octree::octree(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int maxNumParticles, int type) : IndexContainer(cell_size, volMin, volMax, maxNumParticles, type){
	int x_len = abs(volMax.x - volMin.x);
	int y_len = abs(volMax.y - volMin.y);
	int z_len = abs(volMax.z - volMin.z);

	int maxLen = x_len > y_len ? x_len : y_len;
	maxLen = maxLen > z_len ? maxLen : z_len;

	otree = Octree<std::vector<int> >(maxLen);
}

void octree::clear(){
	otree = Octree<std::vector<int> >(maxLen);
}

octree::~octree(){
	
}

void octree::insert(glm::vec3 pos, int index){
	glm::vec3 p = getMinCellPos(pos);
	otree(p.x,p.y,p.z).push_back(index);
}

int octree::find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index){
	int count = 0;
	glm::vec3 cellMid = getMinCellPos(pos);
	glm::vec3 cellPos = getActualCellPos(pos); //Translate into standardized coordinates

	int cx = cellMid.x, cy = cellMid.y, cz = cellMid.z;
	cellMid += 0.5;

	int nx,ny,nz;

	if(cellPos.x>cellMid.x){
		nx = cx+1;
	}else{
		nx = cx-1;
	}
	if(cellPos.y>cellMid.y){
		ny = cy+1;
	}else{
		ny = cy-1;
	}
	if(cellPos.z>cellMid.z){
		nz = cz+1;
	}else{
		nz = cz-1;
	}
	/*auto it = &otree.at(cx,cy,cz);
	m.insert(m.end(),it->begin(), it->end());
	it = &otree.at(nx,cy,cz);
	m.insert(m.end(),it->begin(), it->end());
	it = &otree.at(nx,ny,cz);
	m.insert(m.end(),it->begin(), it->end());
	it = &otree.at(nx,ny,nz);
	m.insert(m.end(),it->begin(), it->end());
	it = &otree.at(cx,ny,nz);
	m.insert(m.end(),it->begin(), it->end());
	it = &otree.at(cx,cy,nz);
	m.insert(m.end(),it->begin(), it->end());
	it = &otree.at(nx,cy,nz);
	m.insert(m.end(),it->begin(), it->end());
	it = &otree.at(cx,ny,cz);
	m.insert(m.end(),it->begin(), it->end());*/
	return count;
}

hashMap::hashMap(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int maxNumParticles, int type) : IndexContainer(cell_size, volMin, volMax, maxNumParticles, type){
	map = boost::unordered_multimap<glm::ivec3, int, boost::hash<glm::ivec3> > (maxNumParticles);
	NCells = new std::pair<hMap::iterator, hMap::iterator>*[NUM_PROCESS];
	for(unsigned int i=0; i<NUM_PROCESS; i++){
		NCells[i] = new std::pair<hMap::iterator, hMap::iterator>[8];
	}
}

hashMap::~hashMap(){
	for(int i=0; i<8; i++){
		delete NCells[i];
	}
	delete NCells;
}

void hashMap::insert(glm::vec3 pos, int index){
	glm::vec3 key = getMinCellPos(pos);
	map.insert(std::make_pair(glm::ivec3(key.x,key.y,key.z),index));
}

void hashMap::erase(glm::vec3 pos, int index){
	glm::vec3 key = getMinCellPos(pos);
	glm::ivec3 ikey = glm::ivec3(key.x,key.y,key.z);

	auto range = map.equal_range(ikey);
	for(auto it = range.first; it!= range.second; it++){
		if(it->second = index){
			map.erase(it);
			break;
		}
	}
}

int hashMap::find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index){
	glm::vec3 cellMid = getMinCellPos(pos);
	glm::vec3 cellPos = getActualCellPos(pos); //Translate into standardized coordinates

	int cx = cellMid.x, cy = cellMid.y, cz = cellMid.z;
	cellMid += 0.5;

	int nx,ny,nz;

	if(cellPos.x>cellMid.x){
		nx = cx+1;
	}else{
		nx = cx-1;
	}
	if(cellPos.y>cellMid.y){
		ny = cy+1;
	}else{
		ny = cy-1;
	}
	if(cellPos.z>cellMid.z){
		nz = cz+1;
	}else{
		nz = cz-1;
	}

	int nProcess = 0; //not parallel (for now...)

	NCells[nProcess][0] = map.equal_range(glm::ivec3(cx,cy,cz)); //Always using a smoothing length of half the cell size;
	NCells[nProcess][1] = map.equal_range(glm::ivec3(nx,cy,cz)); //Only need to find 8 neighboring cells;
	NCells[nProcess][2] = map.equal_range(glm::ivec3(cx,ny,cz)); //http://www.rchoetzlein.com/eng/graphics/fluids.htm
	NCells[nProcess][3] = map.equal_range(glm::ivec3(cx,cy,nz));
	NCells[nProcess][4] = map.equal_range(glm::ivec3(cx,ny,nz));
	NCells[nProcess][5] = map.equal_range(glm::ivec3(nx,ny,nz));
	NCells[nProcess][3] = map.equal_range(glm::ivec3(nx,ny,cz));
	NCells[nProcess][7] = map.equal_range(glm::ivec3(nx,cy,nz));

	int count = 0;
	//int ind = rowLen * index;

	for(int i=0; i<8; i++){
		for(auto j_it=NCells[nProcess][i].first; j_it!=NCells[nProcess][i].second; j_it++){
			if(count == rowLen){break;}
			nArray[count] = j_it->second;
			count++;
		}
	}
	return count;
}

void hashMap::clear(){
	map.clear();
}

grid::grid(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int maxNumParticles, int type) : IndexContainer(cell_size, volMin, volMax, maxNumParticles, type){
	
	gridSize = (max.x)*(max.y)*(max.z);
	if (maxNumParticles < gridSize){
		gridSize = maxNumParticles;
	}
	cells = new int*[gridSize];
	nParticlesInCell = new int[gridSize];
	for(int i=0; i<gridSize; i++){
		cells[i] = new int[MAX_NEIGHBORS];
		nParticlesInCell[i] = 0;
	}
}

grid::~grid(){
	for(int i=0; i<gridSize; i++){
		delete cells[i];
	}
	delete [] cells;
	delete [] nParticlesInCell;
}

void grid::clear(){
	for(int i=0; i < gridSize; i++){
		nParticlesInCell[i] = 0;
	}
}

void grid::insert(glm::vec3 pos, int index){
	glm::vec3 id = getMinCellPos(pos);
	int sid = getSpatialIndex(id);
	if(sid>gridSize){
		int i=0;
	}
	int ind = nParticlesInCell[sid];
	if(ind==MAX_NEIGHBORS){
		return;
	}
	nParticlesInCell[sid]++;
	cells[sid][ind] = index;
}

int grid::find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index){
	glm::vec3 cellMid = getMinCellPos(pos);
	glm::vec3 cellPos = getActualCellPos(pos); //Translate into standardized coordinates

	int cx = cellMid.x, cy = cellMid.y, cz = cellMid.z;
	cellMid += .5;

	int nx,ny,nz;

	if(cellPos.x>cellMid.x){
		nx = cx+1;
	}else{
		nx = cx-1;
	}
	if(cellPos.y>cellMid.y){
		ny = cy+1;
	}else{
		ny = cy-1;
	}
#ifdef IS_3D
	if(cellPos.z>cellMid.z){
		nz = cz+1;
	}else{
		nz = cz-1;
	}
#endif

	int count = 0;
#ifdef IS_3D
	int maxCell = 8;
	int numNeighbor[8];
	int spatialIndex[8];

	spatialIndex[0] = getSpatialIndex(glm::vec3(cx,cy,cz));
	spatialIndex[1] = getSpatialIndex(glm::vec3(nx,cy,cz));
	spatialIndex[2] = getSpatialIndex(glm::vec3(nx,ny,cz));
	spatialIndex[3] = getSpatialIndex(glm::vec3(cx,ny,cz));
	spatialIndex[4] = getSpatialIndex(glm::vec3(cx,cy,nz));
	spatialIndex[5] = getSpatialIndex(glm::vec3(nx,cy,nz));
	spatialIndex[6] = getSpatialIndex(glm::vec3(cx,ny,nz));
	spatialIndex[7] = getSpatialIndex(glm::vec3(nx,ny,nz));

	numNeighbor[4] = nParticlesInCell[spatialIndex[4]];
	numNeighbor[5] = nParticlesInCell[spatialIndex[5]];
	numNeighbor[6] = nParticlesInCell[spatialIndex[6]];
	numNeighbor[7] = nParticlesInCell[spatialIndex[7]];
#else
	int maxCell = 4;
	int numNeighbor[4];
	int spatialIndex[4];

	spatialIndex[0] = getSpatialIndex(glm::vec3(cx,cy,0));
	spatialIndex[1] = getSpatialIndex(glm::vec3(nx,cy,0));
	spatialIndex[2] = getSpatialIndex(glm::vec3(nx,ny,0));
	spatialIndex[3] = getSpatialIndex(glm::vec3(cx,ny,0));
#endif

	numNeighbor[0] = nParticlesInCell[spatialIndex[0]];
	numNeighbor[1] = nParticlesInCell[spatialIndex[1]];
	numNeighbor[2] = nParticlesInCell[spatialIndex[2]];
	numNeighbor[3] = nParticlesInCell[spatialIndex[3]];

	for(int i=0; i<maxCell; i++){
		for(int j=0; j<numNeighbor[i]; j++){
			if(count == rowLen){
				break;
			}
			nArray[count] = cells[spatialIndex[i]][j];
			count++;
		}
	}
	

	return count;
}

int grid::getSpatialIndex(const glm::vec3& pos){
	return ( (unsigned int)((pos.x)*73856093) ^
		(unsigned int)((pos.y)*19349663) ^ 
		(unsigned int)((pos.z)*83492791) ) %  gridSize;
}

/*
grid::grid(float cell_size, glm::vec3 volMin, glm::vec3 volMax, int type) : IndexContainer(cell_size, volMin, volMax, type){
	cells = new int**[max.x+1];
	nParticlesInCell = new int*[max.x+1];
	for(int x=0; x<max.x+1; x++){
		cells[x] = new int*[max.y+1];
		nParticlesInCell[x] = new int[max.y+1];
		for(int y=0; y<max.y+1; y++){
			nParticlesInCell[x][y] = 0;
			cells[x][y] = new int[40];
		}
	}
}

grid::~grid(){
	for(int x=0; x<max.x; x++){
		delete cells[x];
	}
	delete cells;
}

void grid::clear(){
		for(int i=0; i < max.x; i++){
				for(int j=0; j<max.y; j++){
				nParticlesInCell[i][j] = 0;
			}
		}
}

void grid::insert(glm::vec3 pos, int index){
	glm::vec3 id = getMinCellPos(pos);
	int ind = nParticlesInCell[(int)id.x][(int)id.y];
	if(ind==40){
		return;
		//int y = 0;
	}
	nParticlesInCell[(int)id.x][(int)id.y]++;
	cells[(int)id.x][(int)id.y][ind] = index;
}

int grid3D::find(glm::vec3 pos, float radius_sq, int* nArray, int rowLen, int index){
	glm::vec3 cellMid = getMinCellPos(pos);
	glm::vec3 cellPos = getActualCellPos(pos); //Translate into standardized coordinates

	int cx = cellMid.x, cy = cellMid.y, cz = cellMid.z;
	cellMid += 0.5;

	int nx,ny,nz;

	if(cellPos.x>cellMid.x){
		nx = cx+1;
	}else{
		nx = cx-1;
	}
	if(cellPos.y>cellMid.y){
		ny = cy+1;
	}else{
		ny = cy-1;
	}
	if(cellPos.z>cellMid.z){
		nz = cz+1;
	}else{
		nz = cz-1;
	}

	int count = 0;
	int n1=0,n2=0,n3=0,n4=0,n5=0,n6=0,n7=0,n8=0;
	n1 = nParticlesInCell[cx][cy][cz];

	if(nx>=0&&nx<max.x){
		n2 = nParticlesInCell[nx][cy][cz];
		if(ny>=0&&ny<max.y){
			n3 = nParticlesInCell[nx][ny][cz];
			if(nz>=0&&nz<max.z){
				n4 = nParticlesInCell[nx][ny][nz];
			}
		}
	}
	if(ny>=0&&ny<max.y){
		n5 = nParticlesInCell[cx][ny][cz];
		if(nz>=0&&nz<max.z){
			n6 = nParticlesInCell[cx][ny][nz];
		}
	}
	if(nz>=0&&nz<max.z){
		n7 = nParticlesInCell[cx][cy][nz];
		if(nx>=0&&nx<max.x){
			n7 = nParticlesInCell[nx][cy][nz];
		}
	}

	for(int i=0; i<n1; i++){
		if(count == rowLen){break;}
		nArray[count] = cells[cx][cy][cz][i];
		count++;
	}
	for(int i=0; i<n2; i++){
		if(count == rowLen){break;}
		nArray[count] = cells[nx][cy][i];
		count++;
	}
	for(int i=0; i<n3; i++){
		if(count == rowLen){break;}
		nArray[count] = cells[nx][ny][i];
		count++;
	}
	for(int i=0; i<n4; i++){
		if(count == rowLen){break;}
		nArray[count] = cells[cx][ny][i];
		count++;
	}

	return count;
}*/