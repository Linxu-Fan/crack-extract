#ifndef EXTRACTCRACK_H

#define EXTRACTCRACK_H

#include "mpm-fracture/utils.h"
#include "mpm-fracture/damageGradient.h"
#include "mpm-fracture/grid.h"
#include "mpm-fracture/mpm_utils.h"
#include "mpm-fracture/particles.h"
#include "mpm-fracture/weights.h"
#include "voro++.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#include <algorithm>
#include <set>
#include <iterator>












  

// read obj file
struct meshObjFormat readObj(std::string path);


// void writeOffFile(std::vector<Eigen::Vector3d> vertices, std::vector<std::vector<int>> faces, std::string name);




// Struct of particles
struct Point {
    int index = 0; // index of each point
    Eigen::Vector3d pos = { 0, 0, 0 }; // each point's position
    int numVertices = 0; // number of vertices
    std::vector<Eigen::Vector3d> verticsCoor; // vertices' coordinate
    int numFaces = 0; // number of faces
    std::vector<std::vector<int>> verticesFace; // vertices of each face
    std::vector<Eigen::Vector3d> surfaceNormal; // vertices' coordinate
    std::vector<int> neighbour; // neighbour points that share common faces with this point
    std::vector<int> neighbourCalculated; // neighbour points that have already find shared face

    std::vector<int> neighbourSameSide; // neighbour points that are on the same side with this point
    std::vector<int> neighbourOtherSide; // neighbour points that are on the other side with this point

    Point(int iindex, Eigen::Vector3d ipos, int inumVertices, std::vector<Eigen::Vector3d> iverticsCoor, int inumFaces, std::vector<std::vector<int>> iverticesFace, std::vector<Eigen::Vector3d> isurfaceNormal, std::vector<int> ineighbour)
        : index(iindex)
        , pos(ipos)
        , numVertices(inumVertices)
        , verticsCoor(iverticsCoor)
        , numFaces(inumFaces)
        , verticesFace(iverticesFace)
        , surfaceNormal(isurfaceNormal)
        , neighbour(ineighbour)
    {
    }
};


// User-defined wall of Voro++
class wallShell : public voro::wall {
public:
    wallShell(std::map<int, int>* igridMap, std::vector<Grid>* inodesVec, struct parametersSim iparam, int iw_id = -99)
        : gridMap(igridMap)
        , nodesVec(inodesVec)
        , param(iparam)
        , w_id(iw_id) {};

    bool point_inside(double x, double y, double z)
    {
        bool inside = false;

        Eigen::Vector3d pos = { x, y, z };
        Eigen::Vector3i ppIndex = { 0, 0, 0 };

        int numDamage = (*nodesVec).size();
        for (int k = 0; k < numDamage; k++) {
            Eigen::Vector3d nodePos = (*nodesVec)[k].posIndex.cast<double>() * param.dx;
            double diff = (nodePos - pos).norm();

            if (diff <= 0.00001) {
                ppIndex = (*nodesVec)[k].posIndex;
                goto stop0;
            }
        }

        //cout << pos[0] << " " << pos[1] << " " << pos[2] << endl;

    stop0:
        for (int axis = 0; axis < 3; axis++) // three axis directions
        {
            for (int i = -1; i < 2; i++) {
                for (int j = -1; j < 2; j++) {
                    Eigen::Vector3i normal = { 0, 0, 0 };
                    if (axis == 0) {
                        normal = { 0, i, j };
                    } else if (axis == 1) {
                        normal = { i, 0, j };
                    } else {
                        normal = { i, j, 0 };
                    }

                    Eigen::Vector3i nodePosIndex = ppIndex + normal;
                    int ID = calculateID(nodePosIndex[0], nodePosIndex[1], nodePosIndex[2], param.length, param.dx);
                    if ((*gridMap).find(ID) != (*gridMap).end()) // if this point is a neighbour of a fully damaged particle(the distance is smaller than dx)
                    {
                        int eid = (*gridMap)[ID];
                        if ((*nodesVec)[eid].Di == 2) // cannot be shell particles
                        {
                            inside = true;
                            goto stop1;
                        }
                    }
                }
            }
        }

    stop1:
        return inside;
    }

    template <class vc_class>
    inline bool cut_cell_base(vc_class& c, double x, double y, double z)
    {
        Eigen::Vector3d pos = { x / param.dx, y / param.dx, z / param.dx };
        Eigen::Vector3i ppIndex = { 0, 0, 0 };
        ppIndex[0] = round(pos[0]);
        ppIndex[1] = round(pos[1]);
        ppIndex[2] = round(pos[2]);

        //cout << pos[0] << " " << pos[1] << " " << pos[2] << " " << endl;

        for (int axis = 0; axis < 3; axis++) // three axis directions
        {
            for (int pn = 0; pn < 2; pn++) {
                int direction = 2 * pn - 1;
                int distance = 1;
                bool reachBoundary = false;
                do {
                    Eigen::Vector3i normal = { 0, 0, 0 };
                    normal[axis] = direction * distance;

                    Eigen::Vector3i nodePosIndex = ppIndex + normal;

                    int ID = calculateID(nodePosIndex[0], nodePosIndex[1], nodePosIndex[2], param.length, param.dx);
                    if ((*gridMap).find(ID) != (*gridMap).end()) {
                        int eid = (*gridMap)[ID];
                        if ((*nodesVec)[eid].Di == 3) // if this node is a boundary shell
                        {
                            Eigen::Vector3d posNodeOther = normal.cast<double>() * 2 * param.dx;
                            bool tempCut = c.nplane(posNodeOther[0], posNodeOther[1], posNodeOther[2], w_id);
                            reachBoundary = true;
                        }
                    }

                    distance += 1;
                } while (reachBoundary == false);
            }
        }

        return true;
    }

    bool cut_cell(voro::voronoicell& c, double x,
        double y, double z)
    {
        return cut_cell_base(c, x, y, z);
    }
    bool cut_cell(voro::voronoicell_neighbor& c, double x,
        double y, double z)
    {
        return cut_cell_base(c, x, y, z);
    }

private:
    const int w_id;
    std::map<int, int>* gridMap;
    std::vector<Grid>* nodesVec;
    struct parametersSim param;
};

// split a line from a text file
std::vector<std::string> split(const std::string&, const std::string&);

// Set the damage phase of a grid node std::vector into a specific value
void setNodeValue(std::vector<Grid>*, int);

// Find the bounding box boundary nodes and set its damage phase into a specific value
void findBoundaryNodes(std::vector<Particle>*, std::vector<Grid>*, std::map<int, int>*, struct parametersSim, int);

// Read particles' positions and damage phases
void readParticles(std::string, std::vector<Particle>*, bool, struct parametersSim);

// Calculate the damage value of any point and return the value
double ifFullyDamaged(Eigen::Vector3d, parametersSim, std::map<int, int>*, std::vector<Grid>*);

// Store the index of each vertex and its position
int findIndexVertex(Eigen::Vector3d, std::vector<Eigen::Vector3d>*);

// Read voronoi tesselation
void readPoints(std::vector<Point>*);

// Read all structured nodes and calculate the damage gradient
void readParticlesAndCalGradient(std::string, std::vector<Particle>*, parametersSim, std::map<int, int>*, std::vector<Grid>*);

// Find paths between two nodes
bool findPath(Eigen::Vector3d, Eigen::Vector3d, parametersSim, std::vector<Grid>*, std::map<int, int>*, std::vector<int>*);

// Find if a pair of nodes belong to critical nodes. The function return true if one node is a critical node
bool ifCriticalNode(Eigen::Vector3d, parametersSim, std::vector<int>*);

// Find the nearest boundary node of a critical node
Eigen::Vector3i findNearestBoundaryNode(int, std::vector<Point>*, std::vector<int>*, std::vector<Eigen::Vector3i>*, std::map<int, int>*, parametersSim, std::vector<int>*);

// Judge if a pair of points are on different sides of a crack
bool ifTwoSides(int, int, std::vector<Point>*, std::vector<int>*, std::vector<Eigen::Vector3i>*, std::map<int, int>*, parametersSim, std::vector<int>*);



//// Extract the crack surface
//crackSurface extractCrackSurface(std::string, std::vector<Eigen::Vector3d>*, struct parametersSim);






// the following functions clean the mesh and sample points
// the following functions clean the mesh and sample points



// clean faces from voro output
std::tuple< std::vector<Eigen::Vector3d>, std::vector<std::vector<int>>, std::vector<std::vector<int>> > cleanMeshAfterVoro(std::vector<Eigen::Vector3d>* vertices, std::vector<std::vector<int>>* faces);



// insert points between two points
std::vector<Eigen::Vector3d> insertPoints(Eigen::Vector3d startPos, Eigen::Vector3d endPos, double dx);



// mesh struct
struct meshFace
{
	// the index of each face
	int index = -99;

	// each vertex in the face
	std::vector<int> vertices;

	// updated vertices in each edge
	std::vector<int> verticesUpdated;

	// updated vertices in each edge
	std::vector<Eigen::Vector2d> verticesUpdatedCoor2D;

	// sample points index inside of the face
	std::vector<int> sampledPoints;

	// sample points coordinates in 2D
	std::vector<Eigen::Vector2d> sampledPointsCoor2D;

	// inserted points after each vertex
	std::vector<std::vector<int>> insertPointsAfterVertex;

	// key is edge name, value is the first vertex and second vertex
	// the element is string. If the start point and end point are pos1 and pos2 respectively, the edge is {string(pos1) + "*" + string(pos2)} 
	std::map<std::string, std::pair<int, int>> edge;

	// key is vertex index, value is the position where the vertex is in the vertices vector
	std::map<int, int> vertexPosition;

	// key is edge name, value is the edge length
	std::map<std::string, double> edgeLength;

	// neighbouring faces of each edge
	std::map<std::string, std::vector<int>> neighbours;

	// inserted pints between two vertices
	// key is edge name, value is the index of each poin
	std::map<std::string, std::vector<int>> edgeInsertedPoints;

	// indicate if an edge has already been interpolated or not
	std::map<std::string, bool> edgeInsertedOrNot;


	// https://answers.unity.com/questions/1522620/converting-a-3d-polygon-into-a-2d-polygon.html
	// face normal
	Eigen::Vector3d normal = { 0, 0, 0 };

	// first unit vector in the face plane
	Eigen::Vector3d u = { 0, 0, 0 };

	// second unit vector in the face plane
	Eigen::Vector3d v = { 0, 0, 0 };


	// triangulation results
	// store triangles
	std::vector<Eigen::Vector3i> trianglesInThisFace;
};



// refined mesh struct
struct meshFaceRefined
{
	// branching vertices
	std::set<int> branchingEdgeVertices; // vertices that are in the branching region

	// refined vertices
	std::vector<Eigen::Vector3d> verticesRefined;

    // refined faces
	std::vector<Eigen::Vector3i> facesRefined;

    // refined faces
	std::vector<std::vector<int>> facesRefinedFail;

	// neighbouring vertices
	std::vector<std::set<int>> neighbouringVertex;

	// the average normal of each vertex
	std::vector<Eigen::Vector3d> VertexNormal;

};






// refine mesh
struct meshFaceRefined meshRefine(std::vector<Eigen::Vector3d>* vertices, std::vector<std::vector<int>>* faces, std::vector<std::vector<int>>* verticesInfor, double dx); // dx is the refined resolution



// Extract the crack surface
// if a crack surface is found, the crack surface with partial cut, the crack surface with full cut,  each fragment volume in .obj format
std::tuple<bool, meshObjFormat, meshObjFormat, std::vector<meshObjFormat> >  extractCrackSurface(std::vector<Particle>* particlesRaw, struct parametersSim param);







////////////////////////////////
// function for multi-velocity field
////////////////////////////////
crackSurfaceInfor  extractCrackSurface_MultiVel(std::string fileName, std::vector<Eigen::Vector3d>* mpmParticlePositions, struct parametersSim param);







////////////////////////////////
// function for getting partial cracks
////////////////////////////////
meshFaceRefined extractCrackSurface_PartialCracks(std::string fileName, std::vector<Eigen::Vector3d>* mpmParticlePositions, struct parametersSim param);






////////////////////////////////
// cut objects with ftetwild
////////////////////////////////
void cutObject(std::string nameObject, std::vector<meshObjFormat>* fragments, std::vector<meshObjFormat>* fragmentCollision ); // 1) the object's name, ie. bunny_F1, 2) fragments


////////////////////////////////
// triangulate a polygon mesh with libigl
////////////////////////////////
void triangulate_polygonMesh(meshObjFormat* polygonMesh);


////////////////////////////////
// intersection of two objects with mcut
////////////////////////////////
std::vector<std::pair<meshObjFormat, std::pair<std::set<int> , std::map<int, int> > > > intersection_MCUT(trimesh::TriMesh* collisionMesh, meshObjFormat* fragmentVolumeth);



////////////////////////////////
// upsample a triangle mesh from cleaned ftetwild output
////////////////////////////////
meshObjFormat upsampleMesh(parametersSim param, meshObjFormat *inputMesh);



////////////////////////////////
// cut objects with mcut
////////////////////////////////
void cutObject_MCUT(parametersSim param, std::string objectName, meshObjFormat* collisionMeshParent, std::vector<meshObjFormat>* fragments, std::vector<meshObjFormat>* finalFragments);











#endif
