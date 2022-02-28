#pragma once


#include <vector>
#include <Eigen/Core>

struct Object;
struct GenericMesh;
struct Params;
struct meshObjFormat;
struct chevronMarksGridAndMap;
/* Run the fracture simulation using the previously added contacts as boundary conditions
* rbTimeCode is used to distinguish output files from different calls to this method
* returns negative values for errors zero otherwise
* "fractureSurfacePoints" is the points with which to construct the final fracture surface volume 
*/
extern "C" int runFractureSim(Object* pBrb, int timeStepIndex, GenericMesh* finalCrackSurfaceMesh /*output*/ , const Params & simParams);


extern "C" int runFractureSim_MeshCutting(Object* pBrb, int timeStepIndex, std::vector<meshObjFormat>* fragmentsCollision, std::vector<meshObjFormat>* fragmentsCollision_NoTri, std::vector<meshObjFormat>* fragmentsRendering, std::vector<std::vector<Eigen::Vector3d>>* perturbVerticesMag,const Params & simParams);


extern "C" int runFractureSim_OpenVdbFullyCut(Object* pBrb, int timeStepIndex, std::vector<GenericMesh>* finalCrackSurfaceMesh /*output*/,  const Params& simParams);


extern "C" int fixLionMesh();