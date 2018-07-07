#ifndef FLUIDSIMCONTROL_H
#define FLUIDSIMCONTROL_H

/*TODO
	-color density field for particles
	-add grid visualization
*/

#define NUM_PARTICLES 1024
#define NUM_PROCESS 4
#define MAX_NEIGHBORS 140
//#define COMPILE_MOVIE_FRAMES
//#define SAVE_FPS_DATA
//#define SAVE_KIN_DATA
//#define MAX_TIME 3.0
//#define IS_3D
//#define BUILD_MESH
#define NUM_PARALLEL_MESH_REC 10

const glm::mat4 IDENTITY_MAT4 = glm::mat4(1,0,0,0,
										  0,1,0,0,
										  0,0,1,0,
										  0,0,0,1);

#endif