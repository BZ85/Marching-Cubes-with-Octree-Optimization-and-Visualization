/*
 Marching Cubes Research Project, Computer Science, USC
 OpenGL rendering part
 Author Name: Xinjie Zhu
 Advised by Professor Oded Stein
 Basic OpenGL environment code referencs the start code from CSCI 420 Computer Graphics,Computer Science, USC
*/

#include "openGLHeader.h"
#include "glutHeader.h"
#include "openGLMatrix.h"
#include "imageIO.h"
#include "pipelineProgram.h"
#include "vbo.h"
#include "vao.h"

#include "mcFunctions.h"
#include <fstream>
#include <iostream>
#include <cstring>

#include"octree.h"
#include"octreeNode.h"
#include <chrono>

#if defined(WIN32) || defined(_WIN32)
  #ifdef _DEBUG
    #pragma comment(lib, "glew32d.lib")
  #else
    #pragma comment(lib, "glew32.lib")
  #endif
#endif

#if defined(WIN32) || defined(_WIN32)
  char shaderBasePath[1024] = SHADER_BASE_PATH;
#else
  char shaderBasePath[1024] = "../openGLHelper";
#endif

using namespace std;

int mousePos[2]; // x,y screen coordinates of the current mouse position

int leftMouseButton = 0; // 1 if pressed, 0 if not 
int middleMouseButton = 0; // 1 if pressed, 0 if not
int rightMouseButton = 0; // 1 if pressed, 0 if not

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROL_STATE;
CONTROL_STATE controlState = ROTATE;
//CONTROL_STATE controlState = TRANSLATE;

// Transformations of the terrain.
float terrainRotate[3] = { 0.0f, 0.0f, 0.0f }; 
// terrainRotate[0] gives the rotation around x-axis (in degrees)
// terrainRotate[1] gives the rotation around y-axis (in degrees)
// terrainRotate[2] gives the rotation around z-axis (in degrees)
float terrainTranslate[3] = { 0.0f, 0.0f, 0.0f };
float terrainScale[3] = { 1.0f, 1.0f, 1.0f };

// Width and height of the OpenGL window, in pixels.
int windowWidth = 1280;
int windowHeight = 720;
char windowTitle[512] = "Marching Cubes Research Project";



bool ifRecording = false;

OpenGLMatrix matrix;
PipelineProgram * pipelineProgram = nullptr;

VAO * vaoPoint = nullptr;
VAO* vaoTri = nullptr;


VAO* vaoOctreePoint = nullptr;
VAO* vaoOctreeTri = nullptr;

VAO* vaoOctreeWireFrame = nullptr;

//The size of grids of marching cubes, determing the degree of refinement of pruduced surface
//Following values are default values of height, width and length

int height = 15;
int width = 15;
int length = 15;

float isoValue = 0.0f; // iso surface value should be 0 because we use SDF for iput
int numGridVertices = height * width * length;
int numTriangleVertices = 0;

int numOctreeGridVertices = 0;
int numOctreeTriangleVertices = 0;
int numOctreeWireFrames = 0;
 VBO * vboVerticesPos;
 VBO * vboVerticesCol;

 VBO* vboTriangleVerticesPos;
 VBO* vboTriangleVerticesCol;

 VBO* vboOctreeVerticesPos;
 VBO* vboOctreeVerticesCol;

 VBO* vboOctreeTriangleVerticesPos;
 VBO* vboOctreeTriangleVerticesCol;

 VBO* vboOctreeWireFramePos;
 VBO* vboOctreeWireFrameCol;
 // rendering mode: 1 is original marching cubes, 2 is octree marching cubes
 int renderingMode = 2;

 bool showSurfaces = true;
 bool showGridsOrWireFrame = true;


 // Default SDF and OBJ file name
 string SDFfilename = "./SDF/15/torus.raw";
 string objFilename = "./OBJ/15/surface.obj";

 ifstream inputFile("input.txt", std::ios::binary);

// Write a screenshot to the specified filename.
void saveScreenshot(const char * filename)
{
  unsigned char * screenshotData = new unsigned char[windowWidth * windowHeight * 3];
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE, screenshotData);

  ImageIO screenshotImg(windowWidth, windowHeight, 3, screenshotData);

  if (screenshotImg.save(filename, ImageIO::FORMAT_JPEG) == ImageIO::OK)
    cout << "File " << filename << " saved successfully." << endl;
  else cout << "Failed to save file " << filename << '.' << endl;

  delete [] screenshotData;
}


void MultiplyMatrices(int m, int p, int n, const float* A, const float* B, float* C)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            float entry = 0.0;
            for (int k = 0; k < p; k++)
                entry += A[k * m + i] * B[j * p + k];
            C[m * j + i] = entry;
        }
    }
}


// number of images saved to disk so far
int sprite = 0;
void idleFunc()
{
  // save the screenshots to disk (to make the animation).
  
    if (ifRecording) {
        
        if (sprite % 200 == 0) {

            char s[20] = "xxx.jpg";

            s[0] = 48 + ((sprite / 200) % 1000) / 100;
            s[1] = 48 + ((sprite / 200) % 100) / 10;
            s[2] = 48 + (sprite / 200) % 10;

            saveScreenshot(s);

        }
        sprite++;
    }

  // Notify GLUT that it should call displayFunc.
  glutPostRedisplay();
}

void reshapeFunc(int w, int h)
{
  glViewport(0, 0, w, h);

  // When the window has been resized, we need to re-set our projection matrix.
  matrix.SetMatrixMode(OpenGLMatrix::Projection);
  matrix.LoadIdentity();
  // You need to be careful about setting the zNear and zFar. 
  // Anything closer than zNear, or further than zFar, will be culled.
  const float zNear = 0.1f;
  const float zFar = 10000.0f;
  const float humanFieldOfView = 60.0f;
  matrix.Perspective(humanFieldOfView, 1.0f * w / h, zNear, zFar);
  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
}

void mouseMotionDragFunc(int x, int y)
{
  // Mouse has moved, and one of the mouse buttons is pressed (dragging).

  // the change in mouse position since the last invocation of this function
  int mousePosDelta[2] = { x - mousePos[0], y - mousePos[1] };

  switch (controlState)
  {
    // translate the terrain
    case TRANSLATE:
      if (leftMouseButton)
      {
        // control x,y translation via the left mouse button
        terrainTranslate[0] += mousePosDelta[0] * 0.01f;
        terrainTranslate[1] -= mousePosDelta[1] * 0.01f;

        //cout << terrainTranslate[0] << endl;
      }
      if (rightMouseButton)
      {
        // control z translation via the right mouse button
        terrainTranslate[2] += mousePosDelta[1] * 0.01f;
      
      }
      break;

    // rotate the terrain
    case ROTATE:
      if (leftMouseButton)
      {
        // control x,y rotation via the left mouse button
        terrainRotate[0] += mousePosDelta[1];
        terrainRotate[1] += mousePosDelta[0];
      }
      if (rightMouseButton)
      {
        // control z rotation via the right mouse button
        // updated: control y rotation via right mouse button, moving horizontally
        terrainRotate[2] += mousePosDelta[0];
       
      }
      break;

    // scale the terrain
    case SCALE:
      if (leftMouseButton)
      {
        // control x,y scaling via the left mouse button
        terrainScale[0] *= 1.0f + mousePosDelta[0] * 0.01f;
        terrainScale[1] *= 1.0f - mousePosDelta[1] * 0.01f;
      }
      if (rightMouseButton)
      {
        // control z scaling via the right mouse button
        terrainScale[2] *= 1.0f - mousePosDelta[1] * 0.01f;
      }
      break;
  }

  // store the new mouse position
  mousePos[0] = x;
  mousePos[1] = y;
}

void mouseMotionFunc(int x, int y)
{
  // Mouse has moved.
  // Store the new mouse position.
  mousePos[0] = x;
  mousePos[1] = y;
}

void mouseButtonFunc(int button, int state, int x, int y)
{
  // A mouse button has has been pressed or depressed.

  // Keep track of the mouse button state, in leftMouseButton, middleMouseButton, rightMouseButton variables.
  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      leftMouseButton = (state == GLUT_DOWN);
    break;

    case GLUT_MIDDLE_BUTTON:
      middleMouseButton = (state == GLUT_DOWN);
    break;

    case GLUT_RIGHT_BUTTON:
      rightMouseButton = (state == GLUT_DOWN);
    break;
  }

  // Keep track of whether CTRL and SHIFT keys are pressed.
  switch (glutGetModifiers())
  {
    case GLUT_ACTIVE_CTRL:
      controlState = TRANSLATE;
    break;

    case GLUT_ACTIVE_SHIFT:
      controlState = SCALE;
    break;

    // If CTRL and SHIFT are not pressed, we are in rotate mode.
    default:
      controlState = ROTATE;
    break;
  }

  // Store the new mouse position.
  mousePos[0] = x;
  mousePos[1] = y;
}

void keyboardFunc(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 27: // ESC key
      exit(0); // exit the program
    break;

    case ' ':
      cout << "You pressed the spacebar." << endl;

      // start or end recording
      if (ifRecording) ifRecording = false;
      else ifRecording = true;
    break;

    case '1':
        renderingMode = 1;
        break;

    case '2':
        renderingMode = 2;
        break;

    case '3':
        showSurfaces = !showSurfaces;
        break;

    case '4':
        showGridsOrWireFrame = !showGridsOrWireFrame;
        break;

    case 'x':
      // Take a screenshot.
      saveScreenshot("screenshot.jpg");

  



    break;

  }
}

void displayFunc()
{
  // This function performs the actual rendering.

  // First, clear the screen.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set up the camera position, focus point, and the up vector.
  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
  matrix.LoadIdentity();
  matrix.LookAt(1.5, 1.0, -1.5,
                0.0, 0.0, 0.0,
                0.0, 1.0, 0.0);

 // Compute the view matrix for the light direction before the modeling tranformation of tge object
 // Because the light direction should not be effected by the rotating, scaling and translating of the object
  float modelViewMatrix[16];
  float viewMatrix[16];
  
  matrix.GetMatrix(viewMatrix);

  float lightDirection[4] = { 0, 1.0, 0, 0 };
  float viewLightDirection[3]{};
  float result[4]{};

  MultiplyMatrices(4, 4, 1, viewMatrix, lightDirection, result);
  viewLightDirection[0] = result[0];
  viewLightDirection[1] = result[1];
  viewLightDirection[2] = result[2];


  pipelineProgram->SetUniformVariable3fv("viewLightDirection", viewLightDirection);
  

  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
  matrix.Translate(terrainTranslate[0], terrainTranslate[1], terrainTranslate[2]);
  matrix.Rotate(terrainRotate[0], 1, 0, 0);
  matrix.Rotate(terrainRotate[1], 0, 1, 0);

 
  matrix.Rotate(terrainRotate[2], 0, 1, 0);
  

  matrix.Scale(terrainScale[0], terrainScale[1], terrainScale[2]);

  matrix.GetMatrix(modelViewMatrix);


  float normalMatrix[16];
  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
  matrix.GetNormalMatrix(normalMatrix);

  //matrix.LoadIdentity();
  float projectionMatrix[16];
  matrix.SetMatrixMode(OpenGLMatrix::Projection);

  //matrix.Perspective(60.0f, 1.0f , 0.01f, 1000.0f);
  matrix.GetMatrix(projectionMatrix);

  // Upload the modelview and projection matrices to the GPU. Note that these are "uniform" variables.
  pipelineProgram->SetUniformVariableMatrix4fv("modelViewMatrix", GL_FALSE, modelViewMatrix);
  pipelineProgram->SetUniformVariableMatrix4fv("projectionMatrix", GL_FALSE, projectionMatrix);
  pipelineProgram->SetUniformVariableMatrix4fv("normalMatrix", GL_FALSE, normalMatrix);

  float La[4] = { 1.0, 1.0, 1.0, 1.0 };
  float Ld[4] = { 1.0, 1.0, 1.0, 1.0 };
  float Ls[4] = { 1.0, 1.0, 1.0, 1.0 };
  float ka[4] = { 0.1, 0.1, 0.1, 1.0 };
  float kd[4] = { 1.0, 1.0, 1.0, 1.0 };
  float ks[4] = { 0.2, 0.2, 0.2, 1.0 };

  float alpha = 1.0;

  pipelineProgram->SetUniformVariable4fv("La", La);
  pipelineProgram->SetUniformVariable4fv("Ld", Ld);
  pipelineProgram->SetUniformVariable4fv("Ls", Ls);
  pipelineProgram->SetUniformVariable4fv("ka", ka);
  pipelineProgram->SetUniformVariable4fv("kd", kd);
  pipelineProgram->SetUniformVariable4fv("ks", ks);

  pipelineProgram->SetUniformVariablef("alpha", alpha);

   
  // Execute the rendering.
  // Bind the VAO that we want to render. Remember, one object = one VAO. 



  if (renderingMode == 1) { // original marchcing cubes
      
      if (showGridsOrWireFrame) {
          vaoPoint->Bind();
          glDrawArrays(GL_POINTS, 0, numGridVertices);
      }

      if (showSurfaces) {
          vaoTri->Bind();
          glDrawArrays(GL_TRIANGLES, 0, numTriangleVertices);
      }
  }

  else if (renderingMode == 2) { // octree marching cubes

      if (showSurfaces) {
          vaoOctreeTri->Bind();
          glDrawArrays(GL_TRIANGLES, 0, numOctreeTriangleVertices);
      }

      if (showGridsOrWireFrame) {
          vaoOctreeWireFrame->Bind();
          glDrawArrays(GL_LINES, 0, numOctreeWireFrames);
      }
      else {
          vaoOctreePoint->Bind();
          glDrawArrays(GL_POINTS, 0, numOctreeGridVertices);
      }
  }
  // Swap the double-buffers.
  glutSwapBuffers();
}



void initScene(int argc, char *argv[])
{
    /* 
    // The command line parameters has been replaced by input file temporarily
    if (argc >= 2) SDFfilename = argv[1];

    if (argc == 3) objFilename = argv[2];

    */
    if (inputFile) {
        inputFile >> height >> width >> length;
        inputFile >> SDFfilename;
        inputFile >> objFilename;
    }
    else {
        cerr << "Failed to open input file for reading. Using default parameters." << std::endl;
    }

    numGridVertices = height * width * length;

    cout << "height:"<< height << " " << "width:" <<width << " " << "length:" << length << endl;
    cout << "SDF file path: " << SDFfilename << endl;
    cout << "output OBJ file path: " << objFilename << endl;

    inputFile.close();

  // Set the background color.
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Black color. White color.

  // Enable z-buffering (i.e., hidden surface removal using the z-buffer algorithm).
  glEnable(GL_DEPTH_TEST);

  // Enable glPolygonOffset
  glEnable(GL_POLYGON_OFFSET_FILL);

  // Create a pipeline program. 

  pipelineProgram = new PipelineProgram(); // Load and set up the pipeline program, including its shaders.
  // Load and set up the pipeline program, including its shaders.
  if (pipelineProgram->BuildShadersFromFiles(shaderBasePath, "vertexShader.glsl", "fragmentShader.glsl") != 0)
  {
    cout << "Failed to build the pipeline program." << endl;
    throw 1;
  } 
  cout << "Successfully built the pipeline program." << endl;
    
  // Bind the pipeline program that we just created. 
  pipelineProgram->Bind();


  vector <vertex> gridCoordinates;
  vector <vertex> triangleVertexList;
  vector <float> SDFValues;
  vector <vertex> triVertexUniqueList;

  // Following are two parts of SDF code, part A and part B
  // part A is to create a new SDF file into the input file name (.raw)
  // part B is reading the existing SDF file from the input file name (.raw)
  // we should only choose one of two parts to uncomment and run

  /************** SDF code Part A ************************/
  /*
  //This part of code is used to generate a new SDF file, and then immediately produce marching cubes algorithm based on newly-created SDF values
  
   // generate new SDF file and directly draw the surfaces
  for (int i = 0; i < height; i++) // y axis
      for (int k = 0; k < width; k++) // x axis
          for (int j = 0; j < length; j++) { // z axis, we use -z axis
           
              // here you can edit the signed distance funtion to create new SDF file
              // Torus
              float dk = k - 32;
              float di = i - 32;
              float dj = j - 32;
              float qk = sqrt(dk * dk + dj * dj) - 28;
           //   float scalar = sqrt(qk * qk + di * di) - 3;
     
             // float scalar = (i - 15) * (i - 15) + (k - 15) * (k - 15) + (j - 15) * (j - 15) - 210;
             // float scalar = (i - 30) * (i - 30) + (k - 30) * (k - 30) - 500;
             // float scalar = i - 30;
              float scalar = i * i + j * j + k * k - 600;
              SDFValues.push_back(scalar);
          }


  ofstream fileOut(SDFfilename, std::ios::binary | std::ios::trunc);
  if (!fileOut) {
      cerr << "Failed to open SDF file for writing." << std::endl;
     exit(EXIT_FAILURE);
  }

  for (float value : SDFValues) {
      fileOut.write(reinterpret_cast<const char*>(&value), sizeof(float)); 
  }

  fileOut.close();
  */
  /*****************************************************/


  /************** SDF code Part B ************************/
  
   // This part of code is reading the existing SDF file  
     
  // read the SDF file
  ifstream fileIn(SDFfilename, std::ios::binary);
  if (!fileIn) {
      cerr << "Failed to open SDF file for reading." << std::endl;
      exit(EXIT_FAILURE);
  }

  SDFValues.resize(numGridVertices);
  
  fileIn.read(reinterpret_cast<char*>(SDFValues.data()), numGridVertices * sizeof(float));
  fileIn.close();
  

  /*****************************************************/
  
  // Reading or Writing SDF part is over, next is do marching cubes algorithm
  // There are two ways, octree marching cubes and original marching cubes, they will be done in the order
  
  /*****************************************************/
  /* Following is octree marching cubes algorithm */
  vector <vertex> octreeTriangleVertexList;
  vector <vertex> octreeTriVertexUniqueList;

  unordered_map<int, vertex> octreeGridCoordinatesMap;
  // compute the maximum depth according to grid size
  int maxDepth = ceil(log(length * width * height) / log(8)); 
 
  auto startTime = chrono::high_resolution_clock::now();
 
  // run the octree marching cubes algorithm
  // first build the octree and get the root node, then traverse the octree
  
  Octree::octreeMarchingCubes(octreeGridCoordinatesMap, octreeTriangleVertexList, octreeTriVertexUniqueList, SDFValues, maxDepth, height, length, width);
  
  auto endTime = chrono::high_resolution_clock::now();
  chrono::duration<double> duration = endTime - startTime;
  cout << "Octree Marching Cubes Running time: " << duration.count() << std::endl;
  // generate OBJ file for octree marching cubes
  Octree::generateOBJfile(octreeTriangleVertexList, octreeTriVertexUniqueList, objFilename);

  // do OpenGL rendering work for surface produced by octree marching cubes
  vector<float> octreeGridCoordinatesPositions = {};
  vector<float> octreeGridCoordinatesColors = {};

  for (const auto& gridVertex: octreeGridCoordinatesMap) {

      octreeGridCoordinatesPositions.push_back(gridVertex.second.x);
      octreeGridCoordinatesPositions.push_back(gridVertex.second.y);
      octreeGridCoordinatesPositions.push_back(gridVertex.second.z);

      octreeGridCoordinatesColors.push_back(0.0f);
      octreeGridCoordinatesColors.push_back(0.0f);
      octreeGridCoordinatesColors.push_back(0.0f);
      octreeGridCoordinatesColors.push_back(1.0f);
  }

  vboOctreeVerticesPos = new VBO((int)octreeGridCoordinatesPositions.size() / 3, 3, octreeGridCoordinatesPositions.data(), GL_STATIC_DRAW); // 3 values per position
  vboOctreeVerticesCol = new VBO((int)octreeGridCoordinatesColors.size() / 4, 4, octreeGridCoordinatesColors.data(), GL_STATIC_DRAW); // 4 values per color


  numOctreeGridVertices = octreeGridCoordinatesMap.size();
  vector<float> octreeTriangleVerticesPositions = {};
  vector<float> octreeTriangleVerticesColors = {};

  for (int i = 0; i < octreeTriangleVertexList.size(); i++) {

      octreeTriangleVerticesPositions.push_back(octreeTriangleVertexList[i].x);
      octreeTriangleVerticesPositions.push_back(octreeTriangleVertexList[i].y);
      octreeTriangleVerticesPositions.push_back(octreeTriangleVertexList[i].z);

      octreeTriangleVerticesColors.push_back(octreeTriangleVertexList[i].normal[0]);
      octreeTriangleVerticesColors.push_back(octreeTriangleVertexList[i].normal[1]);
      octreeTriangleVerticesColors.push_back(octreeTriangleVertexList[i].normal[2]);
      octreeTriangleVerticesColors.push_back(1.0f);
  }
  numOctreeTriangleVertices = octreeTriangleVertexList.size();

  vboOctreeTriangleVerticesPos = new VBO((int)octreeTriangleVerticesPositions.size() / 3, 3, octreeTriangleVerticesPositions.data(), GL_STATIC_DRAW); // 3 values per position
  vboOctreeTriangleVerticesCol = new VBO((int)octreeTriangleVerticesColors.size() / 4, 4, octreeTriangleVerticesColors.data(), GL_STATIC_DRAW); // 4 values per color
  

  // render the wire frame of octree marching cubes
  vector<float> octreeWireFramePositions = {};
  vector<float> octreeWireFrameColors = {};

  for (int i = 0; i < Octree::lineList.size(); i = i + 8) {
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 0]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 0]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 0]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 1]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 1]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 1]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 1]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 1]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 1]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 2]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 2]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 2]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 2]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 2]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 2]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 3]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 3]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 3]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 3]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 3]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 3]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 0]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 0]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 0]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 4]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 4]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 4]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 5]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 5]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 5]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 5]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 5]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 5]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 6]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 6]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 6]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 6]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 6]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 6]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 7]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 7]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 7]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 7]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 7]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 7]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 4]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 4]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 4]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 0]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 0]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 0]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 4]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 4]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 4]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 1]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 1]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 1]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 5]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 5]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 5]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 2]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 2]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 2]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 6]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 6]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 6]].z);

      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 3]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 3]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 3]].z);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 7]].x);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 7]].y);
      octreeWireFramePositions.push_back(octreeGridCoordinatesMap[Octree::lineList[i + 7]].z);
  
      for (int j = 0; j < 12; j++) {
          octreeWireFrameColors.push_back(0.0);
          octreeWireFrameColors.push_back(0.0);
          octreeWireFrameColors.push_back(0.0);
          octreeWireFrameColors.push_back(1.0);
          octreeWireFrameColors.push_back(0.0);
          octreeWireFrameColors.push_back(0.0);
          octreeWireFrameColors.push_back(0.0);
          octreeWireFrameColors.push_back(1.0);
      }
  
  }

  vboOctreeWireFramePos = new VBO((int)octreeWireFramePositions.size() / 3, 3, octreeWireFramePositions.data(), GL_STATIC_DRAW); // 3 values per position
  vboOctreeWireFrameCol = new VBO((int)octreeWireFrameColors.size() / 4, 4, octreeWireFrameColors.data(), GL_STATIC_DRAW); // 4 values per color
  
  numOctreeWireFrames = octreeWireFramePositions.size() / 3;
  /***************************************************/


  /***************************************************/
  /* Following is original marching cubes algorithm */

   startTime = chrono::high_resolution_clock::now();

  // Do the original marching cubes algorithms to produce the vertex list of marching cubes surface
  McFunctions::marchingCubes(gridCoordinates, triangleVertexList, SDFValues, triVertexUniqueList, isoValue, height, width, length);
  
   endTime = chrono::high_resolution_clock::now();
    duration = endTime - startTime;
  cout << "Original Marching Cubes Running time: " << duration.count() << std::endl;
  
  // generate the obj file (has generated in the octree marching cubes, so here just commented)
 // McFunctions::generateOBJfile(triangleVertexList, triVertexUniqueList, objFilename);

  // Marching cubes has produced grid vertices and triangle vertices
  // Then it's time to put the vertices into position and color vectors
  // And then create VBO for OpenGL rendering
  // Following code are all preparing for the OpenGL rendering work
  // and they are not part of marching cubes algorithms itself
  vector<float> gridCoordinatesPositions = {};
  vector<float> gridCoordinatesColors = {};

  for (int i = 0; i < gridCoordinates.size(); i++) {

      gridCoordinatesPositions.push_back(gridCoordinates[i].x);
      gridCoordinatesPositions.push_back(gridCoordinates[i].y);
      gridCoordinatesPositions.push_back(gridCoordinates[i].z);

      gridCoordinatesColors.push_back(0.0f);
      gridCoordinatesColors.push_back(0.0f);
      gridCoordinatesColors.push_back(0.0f);
      gridCoordinatesColors.push_back(1.0f);
  }

  vector<float> triangleVerticesPositions = {};
  vector<float> triangleVerticesColors = {};

  for (int i = 0; i < triangleVertexList.size(); i++) {

      triangleVerticesPositions.push_back(triangleVertexList[i].x);
      triangleVerticesPositions.push_back(triangleVertexList[i].y);
      triangleVerticesPositions.push_back(triangleVertexList[i].z);

      triangleVerticesColors.push_back(triangleVertexList[i].normal[0]);
      triangleVerticesColors.push_back(triangleVertexList[i].normal[1]);
      triangleVerticesColors.push_back(triangleVertexList[i].normal[2]);
      triangleVerticesColors.push_back(1.0f);
  }
     numTriangleVertices = triangleVertexList.size();

     vboVerticesPos = new VBO((int)gridCoordinatesPositions.size()/3, 3, gridCoordinatesPositions.data(), GL_STATIC_DRAW); // 3 values per position
     vboVerticesCol = new VBO((int)gridCoordinatesColors.size()/4, 4, gridCoordinatesColors.data(), GL_STATIC_DRAW); // 4 values per color

     vboTriangleVerticesPos = new VBO((int)triangleVerticesPositions.size() / 3, 3, triangleVerticesPositions.data(), GL_STATIC_DRAW); // 3 values per position
     vboTriangleVerticesCol = new VBO((int)triangleVerticesColors.size() / 4, 4, triangleVerticesColors.data(), GL_STATIC_DRAW); // 4 values per color
     /***************************************************/


  
  vaoPoint = new VAO();

  vaoTri = new VAO();
 
  
  vaoPoint->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVerticesPos, "position");
  vaoPoint->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVerticesCol, "normal");


  vaoTri->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboTriangleVerticesPos, "position");
  vaoTri->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboTriangleVerticesCol, "normal");


  vaoOctreePoint = new VAO();
  vaoOctreePoint->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboOctreeVerticesPos, "position");
  vaoOctreePoint->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboOctreeVerticesCol, "normal");
  
  vaoOctreeTri = new VAO();
  vaoOctreeTri->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboOctreeTriangleVerticesPos, "position");
  vaoOctreeTri->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboOctreeTriangleVerticesCol, "normal");

  vaoOctreeWireFrame = new VAO();
  vaoOctreeWireFrame->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboOctreeWireFramePos, "position");
  vaoOctreeWireFrame->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboOctreeWireFrameCol, "normal");
  
  // Check for any OpenGL errors.
  std::cout << "GL error status is: " << glGetError() << std::endl;
}

int main(int argc, char *argv[])
{
 /*
 // The command line parameters has been replaced by input file temporarily
  
  if (argc != 1 && argc != 2 && argc!= 3)
  {
    cout << "The arguments are incorrect." << endl;
    exit(EXIT_FAILURE);
  }

  if (argc == 1) {
      cout << "Default SDF and OBJ files are enabled." << endl;
  }

  if (argc == 2) {
      cout << "Default output OBJ file name is enabled." << endl;

  }
  */
  cout << "Initializing GLUT..." << endl;
  glutInit(&argc,argv);

  cout << "Initializing OpenGL..." << endl;

  #ifdef __APPLE__
    glutInitDisplayMode(GLUT_3_2_CORE_PROFILE | GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
  #else
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
  #endif

  glutInitWindowSize(windowWidth, windowHeight);
  glutInitWindowPosition(0, 0);  
  glutCreateWindow(windowTitle);

  cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
  cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << endl;
  cout << "Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;

  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(windowWidth - 1, windowHeight - 1);
  #endif

  // Tells GLUT to use a particular display function to redraw.
  glutDisplayFunc(displayFunc);
  // Perform animation inside idleFunc.
  glutIdleFunc(idleFunc);
  // callback for mouse drags
  glutMotionFunc(mouseMotionDragFunc);
  // callback for idle mouse movement
  glutPassiveMotionFunc(mouseMotionFunc);
  // callback for mouse button changes
  glutMouseFunc(mouseButtonFunc);
  // callback for resizing the window
  glutReshapeFunc(reshapeFunc);
  // callback for pressing the keys on the keyboard
  glutKeyboardFunc(keyboardFunc);

  // init glew
  #ifdef __APPLE__
    // nothing is needed on Apple
  #else
    // Windows, Linux
    GLint result = glewInit();
    if (result != GLEW_OK)
    {
      cout << "error: " << glewGetErrorString(result) << endl;
      exit(EXIT_FAILURE);
    }
  #endif

  // Perform the initialization.
  initScene(argc, argv);

  // Sink forever into the GLUT loop.
  glutMainLoop();
}



