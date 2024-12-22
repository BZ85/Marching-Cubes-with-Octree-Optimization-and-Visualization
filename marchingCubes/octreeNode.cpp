#include "octreeNode.h"

using namespace std;

// OctreeNode constructor
OctreeNode::OctreeNode(int minX, int minY, int minZ, int maxX, int maxY, int maxZ, bool leaf) {
	

	maxIndexX = maxX;
	maxIndexY = maxY;
	maxIndexZ = maxZ;
	minIndexX = minX;
	minIndexY = minY;
	minIndexZ = minZ;
	isLeaf = leaf;

}



