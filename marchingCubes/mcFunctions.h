#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "triangleTable.h"
#include "vertex.h"
using namespace std;

class McFunctions { // define the functions related to marching cubes

public:
	
	// Original Marching Cubes function (comparing to the octree marching cubes function)
	// only need to call this in OpenGL program for producing marching cubes surface, except generating obj file
	// All input vector (list) parameters should be empty except SDF data list
	// once you call this function, then you can take grid coordinates list and triangle vertex list for OpenGL rendering
	static void marchingCubes(vector<vertex>& gridCoordinates, vector<vertex>& triangleVertexList, const vector<float>& SDFValues, vector<vertex>& triVertexUniqueList, float isoValue, int height, int width, int length);

	// generate the grid coordinates vertices list
	static void generateGridCoordinates(vector<vertex>& gridCoordinates, const vector<float>& SDFValues, int height, int width, int length);

	// calculate the triangle vertex on one edge of the cube, v1 and v2 are two endpoints of the edge
	static vertex calculateIntersection(vertex v1, vertex v2, float isoValue);

	// calculate the triangle vertices for each cube, and all the triangle vertices computed are stored in triangleVertexList
	static void traveseCubes(const vector<vertex>& gridCoordinates, vector<vertex>& triangleVertexList, vector<vertex>& triVertexUniqueList, float isoValue, int height, int width, int length);

	// calculate the normal vector for each grid coordinate
	static void calculateNormalForCubeVertex(vector<vertex>& gridCoordinates, int y, int x, int z, int height, int width, int length);
	
	// eliminate the repeated vertices shared by two cubes on the same edge
	static void postProcessingForTriVertices(vector<vertex>& triVertexUniqueList);
	
	// call this function if you need an output obj file of marching cubes surface 
	static void generateOBJfile(const vector<vertex>& triangleVertexList, const vector<vertex>& triVertexUniqueList, string filename);


};





