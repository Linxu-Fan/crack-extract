#pragma once

#include "mpm-fracture/params.h"
#include "mpm-fracture/utils.h"
#include "mpm-fracture/weights.h"

#include <igl/AABB.h>
//#include <igl/polygons_to_triangles.h>
#include <igl/facet_components.h>
#include <igl/remove_unreferenced.h>
#include <igl/gaussian_curvature.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/centroid.h>

// Used by the function Object::getCurvatureAndElementID because
// trimesh library seems broken when locating the closest point.
struct iglTriMesh
{
    // vertices of recentered clean input mesh
    Eigen::MatrixXf V;
    Eigen::MatrixXi F;

    // BVH of recentered mesh
    igl::AABB<Eigen::MatrixXf, 3> aabbTree;

    // per vertex in V
    Eigen::VectorXf meanCurvature;
    Eigen::Vector3f centerOfMass; // of original (i.e. non-recentered) input mesh.

    double volume;

    iglTriMesh(const std::string &path)
    {
        printf("Read mesh: %s\n", path.c_str());

        Eigen::MatrixXf v;
        Eigen::MatrixXi f;

        if (!igl::readOBJ(path, v, f))
        {
            printf("Error: failed to open file: %s\n", path.c_str());
            std::abort();
        }

        V=v;
        F=f;

        preProcess();
    }

    iglTriMesh(const Eigen::MatrixXf &v, const Eigen::MatrixXi &f)
    {
        V=v;
        F=f;
        preProcess();
    }

    void write(const std::string &path)
    {
        printf("Write mesh: %s\n", path.c_str());

        if (!igl::writeOBJ(path, V, F))
        {
            printf("Error: failed to open file: %s\n", path.c_str());
            std::abort();
        }
    }

private:
    void preProcess()
    {

        //computeCurvature();
        
        igl::centroid(V, F, centerOfMass, volume);

        // recenter
        for (int i = 0; i < V.rows(); ++i)
        {
            Eigen::Vector3f x = V.row(i);
           // V.row(i) = x - centerOfMass;
        }

        printf("create AABB tree\n");
        aabbTree.init(V, F);
    }

    void computeCurvature()
    {
        printf("Compute mean curvature of input mesh\n");

        // Alternative discrete mean curvature
        Eigen::MatrixXf HN;
        Eigen::SparseMatrix<float> L, M, Minv;
        igl::cotmatrix(V, F, L);
        igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_FULL, M);
        igl::invert_diag(M, Minv);
        // Laplace-Beltrami of position
        HN = -Minv * (L * V);
        // Extract magnitude as mean curvature
        Eigen::VectorXf H = HN.rowwise().norm();

        // Compute curvature directions via quadric fitting
        Eigen::MatrixXf PD1, PD2;
        Eigen::VectorXf PV1, PV2;
        igl::principal_curvature(V, F, PD1, PD2, PV1, PV2, 8);

        // mean curvature
        meanCurvature = 0.5 * (PV1 + PV2);
    }
};

struct ColliderData
{

    btTriangleIndexVertexArray *iface = nullptr;
    // for collision detection etc. (holds data used by bullet in "iface")
    trimesh::TriMesh *trimeshLowRes = nullptr;
    // bool usedByActiveRigidBody = false;
    btCollisionShape *shape = nullptr;

    void deleteAll()
    {

        if (iface)
        {
            delete iface;
        }

        if (trimeshLowRes)
        {
            delete trimeshLowRes;
        }
    }

    void set(const ColliderData &other)
    {
        shape = other.shape;
        iface = other.iface;
        trimeshLowRes = other.trimeshLowRes;
    }

    void setNull()
    {
        iface = nullptr;
        trimeshLowRes = nullptr;
        shape = nullptr;
    }
};

// defined in main.cpp
extern std::vector<ColliderData> storedData;

struct Object
{
    // parameters for mesh cutting
    /////////////////////////////
    double fragmentVolume = 0;

    /////////////////////////////

    iglTriMesh *iglRenderingMesh = nullptr; // for closest point queries etc.

    float trimeshHighResVolume = 0;
    Eigen::Vector3f centerOfMass;
    trimesh::TriMesh *renderingMesh = nullptr; // for rendering etc.
    trimesh::KDtree *trimeshHighResKDTree = nullptr; // for queries
    float trimeshHighResMaxEdgeLen = 0;



    trimesh::TriMesh *cuttingMesh = nullptr; // for cutting object etc
    meshObjFormat *cuttingMesh_NoTri = nullptr;


    std::vector<Eigen::Vector3d> perturbVerticesMag; // the chevron marks are lost because we cut with the cuttingMesh_NoTri rather than the rendering mesh


    float levelSetVolumeOrig = 0;
    float levelSetVolume = 0;
    openvdb::FloatGrid::Ptr levelSetGridPtr = nullptr;

    float mass = 0;
    const BreakableObjectParams *m_params;
    int fragmentCounter = 0;
    int updateCounter = 0;

    ColliderData collUpdate;
    ColliderData collInUse;
    Eigen::Vector3f updateCOM;
    float updateVolume = 0;
    std::map<int, Eigen::Vector3f> balancedTractions; // MPM boundary conditions

    // MapKey=faceIndex
    // MapValue=traction vector
    std::map<int, Eigen::Vector3f> contactTractions;
    std::map<int, Eigen::Vector3f> elemCtr; // store centroid of surface elements
    std::map<int, float> elemArea;
    float contactDuration = 0;
    btRigidBody *rigidBodyPtr = nullptr;

    void initVolumeMeshData();

    static void storeData(ColliderData collMesh);

    /* notify the BulletWrapper that this objects has a new collision shape
		 * that should be updated, the wrapper must first remove this objects RB
		 * from the dynamics world, then call doCollisionUpdate() and finally
		 * add this object's RB back into the dynamics world
		 */
    bool pendingCollisionUpdate();

    void doCollisionUpdate();

    float getEovNuSq();

    float getCurvatureAndElementID(const Eigen::Vector3f &point, int &closestFaceIndex);

    int getClosestFace(const trimesh::TriMesh *m, const trimesh::KDtree *kd, const trimesh::point &p, float maxdist2);

    trimesh::point getClosestPointOnTriangle(const trimesh::TriMesh *m, int i, const trimesh::point &p);

    // Helper for dist2mesh:
    // Find closest point to p on segment from v0 to v1
    trimesh::point closestPointOnSegment(const trimesh::point &v0, const trimesh::point &v1, const trimesh::point &p);

    void clearContacts();

    /* Add a contact point which will contribute a boundary condition for the fracture simulation
		 * p is the point of impact in this rigid body's local coordinate system (alternatively a triangle-ID can be specified)
		 * d is the direction of the impact (also in local coords), should point towards the inside for collisions
		 * the simulation will run for as many time-steps as required by the contact with the longest duration
		 */

    void addContact(unsigned int tri, Eigen::Vector3f &d, float impulse, float duration);

    /* Check if the current set of contact points is sufficient to run a fracture simulation
		 * returns 0 if the longest contact duration is too short for the time-step size
		 * otherwise returns the total magnitude of deformational forces
		 * the caller should then decide whether the force magnitude warrants running the simulation
     * This function must be called after "calculateBalancedContactTractions"
		 */
    float getTotalContactForce(float fractureSimTimeStepSize);

    // balance tractions --> goal is to have sum(forces)=0 and sum(torques)=0 over the surface
    // reads from member variable contactTractions
    // This function must be called before "getTotalContactForce"
    void calculateBalancedContactTractions();

    // volume of triangle mesh + its volume center of mass
    static btVector3 getMeshVolumeAndCentreOfMass(const trimesh::TriMesh *mesh, float &meshVolume);

    static void decimateMesh(
        trimesh::TriMesh *&mesh,
        const int targetTris);

    static float getCenteredMeshFromVBD(
        trimesh::TriMesh *&pMesh,
        Eigen::Vector3f &centreOfMass,
        const openvdb::FloatGrid::Ptr& gridPtr,
        const float vdbVolToMeshAdaptivity);

    static void setCollMeshFromRenderingMesh(
    ColliderData &coll,
    const trimesh::TriMesh *renderingMesh,
    std::string writeMeshFile,
    unsigned int targetTris,
    const float vdbVolToMeshAdaptivity);

    static void setCollMeshFromRenderingMesh(
        ColliderData &coll,
        const trimesh::TriMesh *renderingMesh,
        openvdb::FloatGrid::Ptr& gridPtr,
        std::string writeMeshFile);

    btRigidBody *createConvexRigidBody(
        openvdb::FloatGrid::Ptr& grid,
        Object &parent,
        std::string writeMeshFile,
        double gridVolume,
        int renderMeshRemeshTarget);

    btRigidBody *createFragmentRB(btRigidBody *parent, btCollisionShape *shape, const btVector3 &btCOM, btScalar mass);

#if 0
	double FractureRB::updateMass(){
		double volume = openvdb::tools::levelSetVolume(*fractSim->getLevelSet().getObjectGrid());
		double mass = volume*mat->getDensity();
		btVector3 inertia;
		rb->getCollisionShape()->calculateLocalInertia(mass, inertia);
		rb->setMassProps(mass, inertia);
		return volume;
	}
#endif

#ifdef USE_LEVEL_SET
int splitFragments(
    const Params& simParams,
    const GenericMesh& crackSurfaceMesh,
    const std::string& outDir,
    std::vector<Object*>& largeFragments,
    std::vector<btRigidBody*>& smallFragments);

#endif

#ifdef USE_MCUT


    // set collision mesh
    static void setCollMesh_toColliderData(
        ColliderData &coll,
        trimesh::TriMesh *collisionMesh_Tri,
        std::string writeMeshFile);

    // split fragments
    int splitFragments_MeshCutting(
        const Params &simParams,
        const std::string &outDir,
        std::vector<meshObjFormat> *fragmentsCollision,
        std::vector<meshObjFormat> *fragmentsCollision_NoTri,
        std::vector<meshObjFormat> *fragmentsRendering,
        std::vector<std::vector<Eigen::Vector3d>> perturbVerticesMagVec,
        std::vector<Object *> &largeFragments,
        std::vector<btRigidBody *> &smallFragments);

#endif

#ifdef USE_LEVEL_SET_FULLYCUT
int splitFragments_OpenVdbFullyCut(
    const Params& simParams,
    const std::vector<GenericMesh>& crackSurfaceMesh,
    const std::string& outDir,
    std::vector<Object*>& largeFragments,
    std::vector<btRigidBody*>& smallFragments);


#endif

    const float vdbVolToMeshAdaptivityForRenderMesh = 0.2;
    const float vdbVolToMeshAdaptivityForCollisionMesh = 0.8;
};




static void meshObjToTriMesh(const meshObjFormat& mesh, trimesh::TriMesh *triMeshOutput)
{
    for (int i = 0; i < mesh.vertices.size(); ++i)
    {
        trimesh::point3 p(mesh.vertices[i][0], mesh.vertices[i][1], mesh.vertices[i][2]);
        (*triMeshOutput).vertices.push_back(p);
    }
    for (int i = 0; i < mesh.faces.size(); ++i)
    {
        trimesh::TriMesh::Face f(mesh.faces[i][0], mesh.faces[i][1], mesh.faces[i][2]);
        (*triMeshOutput).faces.push_back(f);
    }
}




static void TriMeshToMeshObj(trimesh::TriMesh *triMesh, meshObjFormat *&meshOut)
{
    meshOut = new meshObjFormat;
    for (int i = 0; i < triMesh->vertices.size(); ++i)
    {
        Eigen::Vector3d vert = {triMesh->vertices[i][0], triMesh->vertices[i][1], triMesh->vertices[i][2]};
        (*meshOut).vertices.push_back(vert);
    }
    for (int i = 0; i < triMesh->faces.size(); ++i)
    {
        std::vector<int> face;
        face.push_back(triMesh->faces[i][0]);
        face.push_back(triMesh->faces[i][1]);
        face.push_back(triMesh->faces[i][2]);
        (*meshOut).faces.push_back(face);
    }
}





static trimesh::TriMesh addChevronMarksToTriMesh(trimesh::TriMesh *mesh, chevronMarksGridAndMap &chevronMarksGridoutput, Eigen::Vector3f com)
{

    trimesh::TriMesh meshWithChevron;


    std::map<int, int> gridMap = chevronMarksGridoutput.gridMap;
	std::vector<Grid> nodesVec =  chevronMarksGridoutput.nodesVec;
    parametersSim param = chevronMarksGridoutput.param;
    std::vector<Eigen::Vector3d> faceNormal(mesh->faces.size());
    std::vector< std::vector<int>> verticesFace(mesh->vertices.size());
    for (int i = 0; i < mesh->faces.size(); ++i)
    {
        int indexV1 = mesh->faces[i][0];
        Eigen::Vector3d v1 = {mesh->vertices[indexV1][0], mesh->vertices[indexV1][1], mesh->vertices[indexV1][2]};
        int indexV2 = mesh->faces[i][1];
        Eigen::Vector3d v2 = {mesh->vertices[indexV2][0], mesh->vertices[indexV2][1], mesh->vertices[indexV2][2]};
        int indexV3 = mesh->faces[i][2];
        Eigen::Vector3d v3 = {mesh->vertices[indexV3][0], mesh->vertices[indexV3][1], mesh->vertices[indexV3][2]};
        Eigen::Vector3d normal = (v1 - v2).cross(v2 - v3).normalized();
        faceNormal[i] = normal;

        verticesFace[indexV1].push_back(i);
        verticesFace[indexV2].push_back(i);
        verticesFace[indexV3].push_back(i);
    }


    for (int i = 0; i < mesh->vertices.size(); ++i)
    {
        Eigen::Vector3d vertex = {mesh->vertices[i][0], mesh->vertices[i][1], mesh->vertices[i][2]};
        vertex += com.cast<double>();
        vertex -= chevronMarksGridoutput.minCoordinate;

        double pertMag = 0;
        Eigen::Vector3d base = vertex / param.smoothRadius - Eigen::Vector3d::Constant(0.5);
        Eigen::Vector3i ppIndex = base.cast<int>();
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    int ID = calculateID(ppIndex[0] + i, ppIndex[1] + j, ppIndex[2] + k, param.length, param.smoothRadius);
                    if (gridMap.find(ID) != gridMap.end())
                    {
                        int eid = gridMap[ID];
                        struct weightAndDreri  WD = calWeight(param.dx, vertex);
                        double weight = WD.weight(0, i) * WD.weight(1, j) * WD.weight(2, k);
                        pertMag += weight * nodesVec[eid].Di;

                        
                    }

                };

            };
        };

        Eigen::Vector3d pertNormal = {0,0,0};
        for(int m = 0; m < verticesFace[i].size(); m++)
        {
            pertNormal += faceNormal[verticesFace[i][m]];
        }
		vertex += chevronMarksGridoutput.param.pertubationMagitude * pertMag * pertNormal.normalized();


        
        vertex -= com.cast<double>();
        vertex += chevronMarksGridoutput.minCoordinate;

        trimesh::point3 p(vertex[0], vertex[1], vertex[2]);
        meshWithChevron.vertices.push_back(p);
    }


    for (int i = 0; i < mesh->faces.size(); ++i)
    {
        trimesh::TriMesh::Face f(mesh->faces[i][0], mesh->faces[i][1], mesh->faces[i][2]);
        meshWithChevron.faces.push_back(f);
    }

    return meshWithChevron;


}





static void getCenteredMesh_Meshcutting(trimesh::TriMesh *pMesh, btVector3 ComFragment)
{
    // recenter the object's mesh around the center of mass
    for (int i = 0; i < pMesh->vertices.size(); ++i)
    {
        pMesh->vertices[i][0] -= ComFragment[0];
        pMesh->vertices[i][1] -= ComFragment[1];
        pMesh->vertices[i][2] -= ComFragment[2];
    }
}


static void getCenteredMesh_Obj_Meshcutting(meshObjFormat *pMesh, btVector3 ComFragment)
{
    // recenter the object's mesh around the center of mass
    for (int i = 0; i < pMesh->vertices.size(); ++i)
    {
        pMesh->vertices[i][0] -= ComFragment[0];
        pMesh->vertices[i][1] -= ComFragment[1];
        pMesh->vertices[i][2] -= ComFragment[2];
    }
}
