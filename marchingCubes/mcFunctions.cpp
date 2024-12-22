/*
 Marching Cubes Research Project, Computer Science, USC
 Original Marching Cubes algorithm part
 Author Name: Xinjie Zhu
 Advised by Professor Oded Stein
*/

#include "mcFunctions.h"


using namespace std;

// generate the grid coordinates vertices list
void McFunctions::generateGridCoordinates(vector<vertex> & gridCoordinates, const vector<float>& SDFValues, int height, int width, int length)
{
   
    // go through each grid vertex
    for(int i = 0 ; i < height; i++) // y axis
            for (int k = 0; k < width; k++) // x axis
                for (int j = 0; j < length; j++) { // z axis, and we use the -z axis

                vertex currentVertex;

                // calculate the position for each grid coordinate vertex
                currentVertex.x = 1.0 * k / (width - 1);
                currentVertex.y = 1.0 * i / (height - 1);
                currentVertex.z = -1.0 * j / (length - 1); // we use the -z axis, so we need to multiply -1.0

                // the index of this vertex in SDF data list and grid coordinates list
                int vertexIndex = i * width * length + k * length + j;
       
                // get the scalar value from SDF data
                currentVertex.scalar = SDFValues[vertexIndex];
                gridCoordinates.push_back(currentVertex);
           
                
            }

    // calculate the normal for each cube vertex
    // because calculating normals needs to be done after every grid vertex are created and scalars are filled into grids
    // so we do this through a new for loop
    for (int i = 0; i < height; i++) // y axis
        for (int k = 0; k < width; k++) // x axis
            for (int j = 0; j < length; j++) { // z axis, and we use -z axis
                calculateNormalForCubeVertex(gridCoordinates, i, k, j, height, width, length);
            }


}

// calculate the intersection point using interpolation
vertex McFunctions::calculateIntersection(vertex v1, vertex v2, float isoValue)
{
    /* Following part of code references Paul Bourke (1994)
       url: https://paulbourke.net/geometry/polygonise/
   ***********************************************/
    // prevent the surface crack due to floating accuracy
    if (abs(isoValue - v1.scalar) < 0.00001)
        return(v1);
    if (abs(isoValue - v2.scalar) < 0.00001)
        return(v2);
    if (abs(v1.scalar - v2.scalar) < 0.00001)
        return(v1);
    /*************************************************/

    vertex intersectPoint;

    float factor = (isoValue - v1.scalar) / (v2.scalar - v1.scalar);

    // interpolate the positions between two cube vertices
    intersectPoint.x = v1.x + factor * (v2.x - v1.x);
    intersectPoint.y = v1.y + factor * (v2.y - v1.y);
    intersectPoint.z = v1.z + factor * (v2.z - v1.z);

    // the scalar of intersect point is useless, so just assign 0.0 to it
    intersectPoint.scalar = 0.0f;

    // interpolate the normals between two cube vertices
    vector <float> normal;
    float nx = v1.normal[0] + factor * (v2.normal[0] - v1.normal[0]);
    float ny = v1.normal[1] + factor * (v2.normal[1] - v1.normal[1]);
    float nz = v1.normal[2] + factor * (v2.normal[2] - v1.normal[2]);
    normal.push_back(nx);
    normal.push_back(ny);
    normal.push_back(nz);

    // fill the normal into intersect point
    intersectPoint.normal = normal;


    return intersectPoint;
}

/** Calculate the triangle vertices for each cube
 All the triangle vertices computed are stored in triangleVertexList
 The difference between triangleVertexList and triVertexUniqueList is that:
 triangleVertexList is stored for OpenGL rendering (like triangle ABC, BCD, there are vertices repeated)
 triVertexUniqueList is stored for generating obj file, and there are no repeated vertex in this list 
**/
void McFunctions::traveseCubes(const vector<vertex>& gridCoordinates, vector<vertex>& triangleVertexList, vector<vertex>& triVertexUniqueList, float isoValue, int height, int width, int length)
{
    // traverse each cube cell in the grid
    for(int y = 0; y < height - 1; y++)
        for(int x = 0; x < width - 1; x++)
            for (int z = 0; z < length - 1; z++) {

               

                //get the index of each vertex (of current cube) in the vertices list
                int v0 = y * width * length + x * length + z;
                int v1 = y * width * length + x * length + z + 1;
                int v2 = y * width * length + (x + 1) * length + z + 1;
                int v3 = y * width * length + (x + 1) * length + z;

                int v4 = (y + 1) * width * length + x * length + z;
                int v5 = (y + 1) * width * length + x * length + z + 1;
                int v6 = (y + 1) * width * length + (x + 1) * length + z + 1;
                int v7 = (y + 1) * width * length + (x + 1) * length + z;

                /* Following part of code references Paul Bourke (1994)
                 url: https://paulbourke.net/geometry/polygonise/
                 ***********************************************/
                
                //calculate case index (from 0 - 255), based on the comparation of isoValue with each cube vertex's SDF value 
                int caseIndex = 0;
                if (gridCoordinates[v0].scalar < isoValue) caseIndex += 1;
                if (gridCoordinates[v1].scalar < isoValue) caseIndex += 2;
                if (gridCoordinates[v2].scalar < isoValue) caseIndex += 4;
                if (gridCoordinates[v3].scalar < isoValue) caseIndex += 8;

                if (gridCoordinates[v4].scalar < isoValue) caseIndex += 16;
                if (gridCoordinates[v5].scalar < isoValue) caseIndex += 32;
                if (gridCoordinates[v6].scalar < isoValue) caseIndex += 64;
                if (gridCoordinates[v7].scalar < isoValue) caseIndex += 128;

                // use case index to fetch an hexadecimal number
                //  the Nth bit of hexadecimal number determines whether the vertex on the Nth edge in current cube will be used to render triangle
                // if the number on Nth bit is 1, then calculate the intersect vertex on this edge, otherwise not
                vertex defaultV;
                vector <vertex> triangleVertexListForThisCube(12, defaultV);
                if (edgeTable[caseIndex] & 1) triangleVertexListForThisCube[0] = calculateIntersection(gridCoordinates[v0], gridCoordinates[v1], isoValue);
                if (edgeTable[caseIndex] & 2) triangleVertexListForThisCube[1] = calculateIntersection(gridCoordinates[v1], gridCoordinates[v2], isoValue);
                if (edgeTable[caseIndex] & 4) triangleVertexListForThisCube[2] = calculateIntersection(gridCoordinates[v2], gridCoordinates[v3], isoValue);
                if (edgeTable[caseIndex] & 8) triangleVertexListForThisCube[3] = calculateIntersection(gridCoordinates[v3], gridCoordinates[v0], isoValue);

                if (edgeTable[caseIndex] & 16) triangleVertexListForThisCube[4] = calculateIntersection(gridCoordinates[v4], gridCoordinates[v5], isoValue);
                if (edgeTable[caseIndex] & 32) triangleVertexListForThisCube[5] = calculateIntersection(gridCoordinates[v5], gridCoordinates[v6], isoValue);
                if (edgeTable[caseIndex] & 64) triangleVertexListForThisCube[6] = calculateIntersection(gridCoordinates[v6], gridCoordinates[v7], isoValue);
                if (edgeTable[caseIndex] & 128) triangleVertexListForThisCube[7] = calculateIntersection(gridCoordinates[v7], gridCoordinates[v4], isoValue);

                if (edgeTable[caseIndex] & 256) triangleVertexListForThisCube[8] = calculateIntersection(gridCoordinates[v0], gridCoordinates[v4], isoValue);
                if (edgeTable[caseIndex] & 512) triangleVertexListForThisCube[9] = calculateIntersection(gridCoordinates[v1], gridCoordinates[v5], isoValue);
                if (edgeTable[caseIndex] & 1024) triangleVertexListForThisCube[10] = calculateIntersection(gridCoordinates[v2], gridCoordinates[v6], isoValue);
                if (edgeTable[caseIndex] & 2048) triangleVertexListForThisCube[11] = calculateIntersection(gridCoordinates[v3], gridCoordinates[v7], isoValue);



                // push the triangle vertices into triVertexUniqueList, preparing for outputing the obj file
                if (edgeTable[caseIndex] & 1) triVertexUniqueList.push_back(triangleVertexListForThisCube[0]);
                if (edgeTable[caseIndex] & 2) triVertexUniqueList.push_back(triangleVertexListForThisCube[1]);
                if (edgeTable[caseIndex] & 4) triVertexUniqueList.push_back(triangleVertexListForThisCube[2]);
                if (edgeTable[caseIndex] & 8) triVertexUniqueList.push_back(triangleVertexListForThisCube[3]);

                if (edgeTable[caseIndex] & 16) triVertexUniqueList.push_back(triangleVertexListForThisCube[4]);
                if (edgeTable[caseIndex] & 32) triVertexUniqueList.push_back(triangleVertexListForThisCube[5]);
                if (edgeTable[caseIndex] & 64) triVertexUniqueList.push_back(triangleVertexListForThisCube[6]);
                if (edgeTable[caseIndex] & 128) triVertexUniqueList.push_back(triangleVertexListForThisCube[7]);

                if (edgeTable[caseIndex] & 256) triVertexUniqueList.push_back(triangleVertexListForThisCube[8]);
                if (edgeTable[caseIndex] & 512) triVertexUniqueList.push_back(triangleVertexListForThisCube[9]);
                if (edgeTable[caseIndex] & 1024) triVertexUniqueList.push_back(triangleVertexListForThisCube[10]);
                if (edgeTable[caseIndex] & 2048) triVertexUniqueList.push_back(triangleVertexListForThisCube[11]);
                /*************************************************************************************************/
       
                /* The following for loop partially references Paul Bourke (1994)
                 url: https://paulbourke.net/geometry/polygonise/ */
                // Get the triangle vertex rendering order from triTable
                // And push them into triangleVertexList
                for (int i = 0; triTable[caseIndex][i] != -1; i++) {

                    int index = triTable[caseIndex][i];

                    vertex triV;
                    triV.x = triangleVertexListForThisCube[index].x;
                    triV.y = triangleVertexListForThisCube[index].y;
                    triV.z = triangleVertexListForThisCube[index].z;
                    vector <float> normal;
                    normal.push_back(triangleVertexListForThisCube[index].normal[0]);
                    normal.push_back(triangleVertexListForThisCube[index].normal[1]);
                    normal.push_back(triangleVertexListForThisCube[index].normal[2]);

                    triV.normal = normal;
                    triangleVertexList.push_back(triV);
                  
                }

    

            }



}

// calculate the gradient as the normal vector for each grid coordinate
void McFunctions::calculateNormalForCubeVertex(vector<vertex>& gridCoordinates, int y, int x, int z, int height, int width, int length)
{
    // the normal of each cube vertex is the gradient 

    int v, vUp, vDown, vLeft, vRight, vFront, vBehind;
     v = y * width * length + x * length + z;

     if (y + 1 < height) vUp = (y + 1) * width * length + x * length + z;
     else vUp = v;
     if (y - 1 >= 0) vDown = (y - 1) * width * length + x * length + z;
     else vDown = v;

     if (z + 1 < length) vLeft = y * width * length + x * length + (z + 1);
     else vLeft = v;
     if (z - 1 >= 0) vRight = y * width * length + x * length + (z - 1);
     else vRight = v;

     if (x + 1 < width) vFront = y * width * length + (x + 1) * length + z;
     else vFront = v;
     if (x - 1 >= 0) vBehind = y * width * length + (x - 1) * length + z;
     else vBehind = v;

    vector<float> gradient;
    float gx = 0.0, gy= 0.0, gz = 0.0;

    gx = (gridCoordinates[vFront].scalar - gridCoordinates[vBehind].scalar) / (2.0 / (width - 1));
    gy = (gridCoordinates[vUp].scalar - gridCoordinates[vDown].scalar) / (2.0 / (height - 1));
    
    // Because we use -z axis, the z coordinate of normal also needs to be multiply -1.0
    gz = -1.0 * (gridCoordinates[vLeft].scalar - gridCoordinates[vRight].scalar) / (2.0 / (length - 1));

    // normalize the gradient
    float vectorLength = sqrt(gx * gx + gy * gy + gz * gz);
    gx = gx / vectorLength;
    gy = gy / vectorLength;
    gz = gz / vectorLength;

    gradient.push_back(gx);
    gradient.push_back(gy);
    gradient.push_back(gz);

    gridCoordinates[v].normal = gradient;
}

// eliminate the repeated vertices shared by two cubes on the same edge
void McFunctions::postProcessingForTriVertices(vector<vertex>& triVertexUniqueList) {

    unordered_set <vertex> vertexSet; 
    for (int i = 0; i < triVertexUniqueList.size(); i++) {
        vertexSet.insert(triVertexUniqueList[i]);
    }
    triVertexUniqueList.clear();
    for (const auto & v: vertexSet) {
        triVertexUniqueList.push_back(v);
    }

}

// Original Marching Cubes function, only need to call this in OpenGL program for producing marching cubes surface, except generating obj file
// All input vector (list) parameters should be empty except SDF data list
// once you call this function, then you can take grid coordinates list and triangle vertex list for OpenGL rendering
void McFunctions::marchingCubes(vector<vertex>& gridCoordinates, vector<vertex>& triangleVertexList, const vector<float>& SDFValues, vector<vertex>& triVertexUniqueList, float isoValue, int height, int width, int length) {

    // Firstly generate the grids coordinates vertices
    generateGridCoordinates(gridCoordinates, SDFValues, height, width, length);
   
    // Then traverse each cube cell to get the triangle surface vertex list. And triVertexUnqiueList is for generating obj file
    traveseCubes(gridCoordinates, triangleVertexList, triVertexUniqueList, isoValue, height, width, length);
    
    postProcessingForTriVertices(triVertexUniqueList);
}

// call this function if you need an output obj file of marching cubes surface 
void McFunctions::generateOBJfile(const vector<vertex>& triangleVertexList, const vector<vertex>& triVertexUniqueList, string filename) {
    
    // Firstly we need to generate the face list, storing the index of triangle vertices
    // We have triangle vertex list (OpenGL rendering order with repeated vertices)
    // We also have triangle unqiue vertex list, compare two list to get the face list

    vector<int> faceList = {};
  
    // create face list based on OpenGL rendering order
    // Using unorder_map (hash map) is more efficient
    unordered_map<vertex, int> vertexMap;

    // put each vertex into hash map
    for (int i = 0; i < triVertexUniqueList.size(); i++) {
        vertexMap[triVertexUniqueList[i]] = i;
    }
    // check the map to get index
    for (int i = 0; i < triangleVertexList.size(); i++) {

        auto iterator = vertexMap.find(triangleVertexList[i]);
        int vertexIndex = iterator->second;
        faceList.push_back(vertexIndex);
    }



    // write the vertex position, normal and face into obj file
    std::ofstream objFile(filename, std::ios::trunc);

    // write vertex positions
    for (vertex v1 : triVertexUniqueList) {
        objFile << "v " << v1.x << " " << v1.y << " " << v1.z << "\n";
    }

    // write normals
    for (vertex v1: triVertexUniqueList) {
        if (v1.normal.size() == 3) {
            objFile << "vn " << v1.normal[0] << " "
                << v1.normal[1] << " "
                << v1.normal[2] << "\n";
        }
        else {
            objFile << "vn 0 0 0\n";
        }
    }

    //write faces
    for (int i = 0; i < faceList.size(); i += 3) {
        int v1 = faceList[i] + 1;    // the index of OBJ start from 1, not 0 
        int v2 = faceList[i + 1] + 1;
        int v3 = faceList[i + 2] + 1;

        objFile << "f " << v1 << "//" << v1 << " "
            << v2 << "//" << v2 << " "
            << v3 << "//" << v3 << "\n";
    }

     objFile.close();
     cout << "obj file saved successfully!" << endl;

}


