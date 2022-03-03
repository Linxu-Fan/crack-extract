#pragma once

// std

// this include must come first in our program to ensure that "#define _USE_MATH_DEFINES" is specified
// before "#include <cmath>". Otherwise, math constants will not be defined
#ifdef _WIN32
// https://docs.microsoft.com/en-us/cpp/c-runtime-library/math-constants?redirectedfrom=MSDN&view=vs-2019
#define _USE_MATH_DEFINES
#include <cmath>
#endif


#include <omp.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <streambuf>
#include <string>
#include <vector>
#include <stdio.h> 

// trimesh
#include "KDtree.h"
#include "TriMesh.h"
#include "lineqn.h"

// pystring
#include "pystring.h"

// eigen
#include <Eigen/Core>
#include <Eigen/Dense>

// bullet
#ifdef _WIN32
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <btBulletDynamicsCommon.h>
#else
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <btBulletDynamicsCommon.h>
#include <BulletDynamics/ConstraintSolver/btSequentialImpulseConstraintSolverMt.h>
#include <BulletCollision/CollisionDispatch/btCollisionDispatcherMt.h>
#include <BulletDynamics/Dynamics/btDiscreteDynamicsWorldMt.h>
//btDiscreteDynamicsWorldMt
#endif

// openvdb
#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetFracture.h>
#include <openvdb/tools/LevelSetMeasure.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/tools/VolumeToMesh.h>


// libigl
#include "igl/upsample.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/polygons_to_triangles.h>
        


// triangle.h
#define REAL double
#define VOID void
#include "triangle/triangle.h"


#include "mpm-fracture/mpm_utils.h"



#define MULTI_BULLET_THREADED 1
#define MULTI_MPM_THREADED 1


#define TIMEASNAME // the previous method for assigning name is too long. Use time as the name

// Using mcut to generate fragments by enabling this marco, 
//#define USE_MCUT

// This corresponds to the method we used before. Partial cracks are enabled.
#define USE_LEVEL_SET

// This is the new method. No partial crack.
//#define USE_LEVEL_SET_FULLYCUT

// asertion macros
#define ASSERT_(Expr, Msg, ...)                                      \
    if (!(Expr)) {                                                   \
        std::printf("Assert failed:\t\n");                           \
        std::printf("Expected:\t%s\n", #Expr);                       \
        std::printf("Source:\t\t%s, line %d\n", __FILE__, __LINE__); \
        std::printf(Msg, ##__VA_ARGS__);                             \
        std::abort();                                                \
    }

#define ASSERT(Expr) ASSERT_(Expr, " \b")

// Uncomment if you wish do dump (large) openvdb grids to file for debugging
// #define DUMP_DEBUG_MESH_DATA_TO_FILE 1

///
#define PROFILING_BUILD

#if defined(PROFILING_BUILD)
#include <chrono>
#include <stack>
#include <memory>

#define TIMESTACK_PUSH(name) \
if(g_timestack.empty()){ \
    g_timestack.push(std::unique_ptr<mini_timer>(new mini_timer(name, &g_time_profiles)));\
}else{\
    g_timestack.push(std::unique_ptr<mini_timer>(new mini_timer(g_timestack.top()->get_name() + "::" + name, &g_time_profiles)));\
}

#define TIMESTACK_POP() \
    g_timestack.pop()

#define TIMESTACK_RESET()                 \
    while (!g_timestack.empty())          \
    {                                     \
        g_timestack.top()->set_invalid(); \
        g_timestack.pop();                \
    }



#define SCOPED_TIMER(name) \
  raii_timer ttt_(name);

#else
#define SCOPED_TIMER(name)
#define TIMESTACK_PUSH(name)
#define TIMESTACK_POP()
#define TIMESTACK_RESET()
#endif

#if defined(PROFILING_BUILD)

class mini_timer
{
    std::chrono::time_point<std::chrono::steady_clock> m_start;
    std::chrono::time_point<std::chrono::steady_clock> m_pause_start;
    std::chrono::time_point<std::chrono::steady_clock> m_pause_end;

    std::string m_name;
    bool m_valid = true;
    bool m_was_paused = false;
    std::map<std::string, std::vector<unsigned long long>>* tmap;

public:
    mini_timer(const std::string &name, std::map<std::string, std::vector<unsigned long long>>* tmap_) : 
    m_start(std::chrono::steady_clock::now()), 
    m_name(name),
    tmap(tmap_)
    {
    }

    ~mini_timer()
    {
        if (m_valid)
        {
            const std::chrono::time_point<std::chrono::steady_clock> now = std::chrono::steady_clock::now();
            const std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_start);
            unsigned long long elapsed_ = elapsed.count();
            unsigned long long elapsed_pause_time_ = 0;

            if (m_was_paused)
            {
                std::chrono::milliseconds elapsed_pause_time = std::chrono::duration_cast<std::chrono::milliseconds>(m_pause_end - m_pause_start);
                elapsed_pause_time_ = elapsed_pause_time.count();
                elapsed_ -= elapsed_pause_time_; // remove noise due to nested tasks
            }

            //printf("[PROFILE]: %s (%llums | %llums)\n", m_name.c_str(), elapsed_, elapsed_pause_time_);

            (*tmap)[m_name].push_back(elapsed_);
        }
    }

    void pause()
    {
        m_was_paused = true;
        m_pause_start = std::chrono::steady_clock::now();
    }

    void resume()
    {
        m_pause_end = std::chrono::steady_clock::now();
    }

    void set_invalid()
    {
        m_valid = false;
    }
    std::string get_name()
    {
        return m_name;
    }
    void set_name(const std::string name)
    {
        m_name = name;
    }
};

extern std::map<std::string, std::vector<unsigned long long>> g_time_profiles;
extern std::stack<std::unique_ptr<mini_timer>> g_timestack;

struct raii_timer{
    raii_timer(const std::string& name)
    {
        TIMESTACK_PUSH(name);
    }

    ~raii_timer()
    {
        TIMESTACK_POP();
    }
};

#endif

extern std::vector<std::pair<std::string,int>> g_mpmFractureSimTimestepCounts;
extern std::map<std::string,double> g_rigidBodyVolumes;

// return two vectors which define the crack surface
struct crackSurface {
    std::vector<Eigen::Vector3d> vertices; // vertices of the crack surface
    std::vector<std::vector<int>> faces; // faces of the crack surface

    crackSurface(std::vector<Eigen::Vector3d> ivertices, std::vector<std::vector<int>> ifaces)
        : vertices(ivertices)
        , faces(ifaces)
    {
    }

    crackSurface() { }
};


// name of rigid body object
static std::string &getName(const btCollisionObject *obj)
{
    return *((std::string *)obj->getUserPointer());
}


struct meshObjFormat
{
	// input mesh vertices and faces
	std::vector<Eigen::Vector3d> vertices;
	std::vector<std::vector<int>> faces;
	std::vector<int> faceFromVoroCell; // indicate the voronoi cell which the face belongs to
	std::vector<int> faceFromtheOtherVoroCell; // indicate the voronoi cell which the face belongs to(the other side)
};



// refined mesh struct
struct crackSurfaceInfor
{
    int numOfFragments;

	// refined vertices
	std::vector<Eigen::Vector3d> vertices;

    // refined faces
	std::vector<std::vector<int>> faces;

	// refined vertices
	std::vector<Eigen::Vector3d> facesNormal;

};



struct GenericMesh {
    GenericMesh()
    {
    }

    // same as "crackSurfaceMesh.faces" but only contains triangles
    // NOTE: do not populate this vector yourself, use "triangulateGenericMesh()".
    // Indices point into "crackSurfaceMesh.vertices"
    std::vector<int> triangulatedFaces;

    crackSurface m;
};

static std::string readTextFile(const std::string& fpath)
{
    std::ifstream t(fpath.c_str());
    std::string str((std::istreambuf_iterator<char>(t)),
        std::istreambuf_iterator<char>());
    return str;
}

static btVector3 toBullet(const trimesh::point& in)
{
    return btVector3(in[0], in[1], in[2]);
}

static btVector3 toBullet(const Eigen::Vector3f& in)
{
    return btVector3(in[0], in[1], in[2]);
}

static trimesh::point toTrimesh(const Eigen::Vector3f& in)
{
    return trimesh::point((in[0]), (in[1]), (in[2]));
}

static Eigen::Vector3f toEigen(const btVector3& in)
{
    return Eigen::Vector3f(in[0], in[1], in[2]);
}

static Eigen::Vector3f toEigen(const trimesh::point& in)
{
    return Eigen::Vector3f(in[0], in[1], in[2]);
}

static btVector3 fromEigen(const Eigen::Vector3f& in)
{
    return btVector3(in[0], in[1], in[2]);
}

static btVector3 fromTrimesh(const trimesh::point& in)
{
    return btVector3(in[0], in[1], in[2]);
}

static trimesh::point toTrimesh(const btVector3& in)
{
    return trimesh::point(in[0], in[1], in[2]);
}

extern void toTrimesh(trimesh::TriMesh& tm, const GenericMesh& gm);

extern void writeVDBMesh(const char* filename, openvdb::tools::VolumeToMesh& mesh);
extern void saveOpenVDBGrids(const openvdb::GridPtrVec& grids, const std::string& outDir, const std::string& name);
extern void saveOpenVDBGrid(const openvdb::FloatGrid::Ptr& grid, const std::string& outDir);

extern openvdb::FloatGrid::Ptr meshToVDBLevelSetGrid(
    const trimesh::TriMesh* mesh,
    const float voxelSize,
    const std::string& name,
    float& volume);

// converts "crackSurfaceMesh" to a level set. This is produces to a thin
// volume around the provided surface, whose thickness is the diagonal length of
// voxel. The diagonal length is determined by the voxel size.
extern openvdb::FloatGrid::Ptr crackSurfaceMeshToLevelSetGrid(
    const GenericMesh& crackSurfaceMesh,
    const float voxelSize,
    const std::string& name);


// convert a mesh with chevron marks to openvdb grid
extern openvdb::FloatGrid::Ptr chevronMeshToLevelSetGrid(
    trimesh::TriMesh& chevronMesh,
    const float voxelSize,
    Eigen::Vector3f com,
    const std::string name);

extern void cleanGenericMesh(GenericMesh& gm);
// updates the member-variable called "triangles" in "m"
extern void triangulateGenericMesh(GenericMesh& m);









///////////////////////////////////////////////////////
// mcut dependency




#include "mcut/mcut.h"

#include <map>
#include <stdio.h>
#include <stdlib.h>
// libigl dependencies
#include <Eigen/Core>
#include <igl/barycentric_coordinates.h>
#include <igl/barycentric_interpolation.h>
#include <igl/remove_unreferenced.h>
#include <igl/facet_components.h>


struct InputMesh
{
    // variables for reading .obj file data with libigl
    std::vector<std::vector<double>> V, TC, N;
    std::vector<std::vector<int>> F, FTC, FN;
    std::vector<std::tuple<std::string, unsigned, unsigned>> FM;

    // variables for mesh data in a format suited for MCUT
    std::string fpath;                      // path to mesh file
    std::vector<uint32_t> faceSizesArray;   // vertices per face
    std::vector<uint32_t> faceIndicesArray; // face indices
    std::vector<double> vertexCoordsArray;  // vertex coords
};

///////////////////////////////////////////////////////





static void getLargestConnectedComponent(trimesh::TriMesh *&pMesh)
{
    Eigen::MatrixXd ftetwildMesh_V(pMesh->vertices.size(), 3);  
    Eigen::MatrixXi ftetwildMesh_F(pMesh->faces.size(), 3);;  
    for(int i = 0; i < pMesh->vertices.size(); i++)
    {
        Eigen::RowVector3d vert = {pMesh->vertices[i][0], pMesh->vertices[i][1], pMesh->vertices[i][2]};
        ftetwildMesh_V.row(i) = vert;
    }
    for(int i = 0; i < pMesh->faces.size(); i++)
    {
        Eigen::RowVector3i face = {pMesh->faces[i][0], pMesh->faces[i][1],pMesh->faces[i][2] };
        ftetwildMesh_F.row(i) = face;
    }



    Eigen::MatrixXd V;  
    Eigen::MatrixXi F;  


    { // remove any unreference vertices (if any) from the surface mesh of the tet-mesh
        Eigen::VectorXi _1;
        Eigen::VectorXi _2;
        igl::remove_unreferenced(ftetwildMesh_V, ftetwildMesh_F, V, F, _1, _2);
    }

        {
            // We depend on "fTetWild", which sometimes produces multiple
            // connected components within a single surface mesh. This is
            // problematic when computing geometric quantities on the
            // remeshed input mesh like curviture. Also, MCUT does not take input
            // meshes that contain multiple connected components.
            // Thus, our solution is to keep only the largest connected component.
            Eigen::VectorXi F_cc;
            igl::facet_components(F, F_cc);

            std::map<int, int> CC_to_face_count;
            for (int i = 0; i < F_cc.rows(); ++i)
            {
                int cc_id = F_cc(i);
                if (CC_to_face_count.count(cc_id) == 0)
                {
                    CC_to_face_count[cc_id] = 0; // init
                }

                CC_to_face_count[cc_id] += 1;
            }

            printf("[PROG]: Found %d connected components in input mesh\n", (int)CC_to_face_count.size());

            if (CC_to_face_count.size() > 1)
            {
                int largest_CC = 0;
                int faceCountMax = 0;

                for (auto &e : CC_to_face_count)
                {
                    if (e.second > faceCountMax)
                    {
                        faceCountMax = e.second;
                        largest_CC = e.first;
                    }

                    //printf("[PROG]: Connected component %d has %d faces \n", (int)e.first, (int)e.second);
                }

                {
                    // create a new mesh which contains the faces of the largest connected
                    // component. this new mesh will be the formal "input" mesh that we use.
                    Eigen::MatrixXi newF(CC_to_face_count.at(largest_CC), 3);
                    int face_index = 0;
                    for (int i = 0; i < F_cc.rows(); ++i)
                    {
                        int cc_id = F_cc(i);
                        if (cc_id == largest_CC)
                        {
                            newF.row(face_index) = F.row(i);
                            face_index++;
                        }
                    }

                    { // keep only the vertices that are used by the largest connected component.
                        Eigen::MatrixXd cleanV;
                        Eigen::MatrixXi cleanF;
                        Eigen::VectorXi _1;
                        Eigen::VectorXi _2;
                        igl::remove_unreferenced(V, newF, cleanV, cleanF, _1, _2);
                        V = cleanV;
                        F = cleanF;
                    }
                }
            }

            printf("[PROG]: Final input mesh vertices = %d, faces = %d\n", (int)V.rows(), (int)F.rows());
        }


    pMesh->clear();
    for(int i = 0; i < V.rows(); i++)
    {
        trimesh::vec3 vert;
        vert[0] = V(i, 0);
        vert[1] = V(i, 1);
        vert[2] = V(i, 2);
        pMesh->vertices.push_back(vert);
    }
    for(int i = 0; i < F.rows(); i++)
    {
        trimesh::TriMesh::Face face;
        face[0] = F(i,0);
        face[1] = F(i,1);     
        face[2] = F(i,2);
        pMesh->faces.push_back(face);
    }


}





static meshObjFormat cleanFtetWild(meshObjFormat *ftetwildMesh)
{
    meshObjFormat outMesh;
    Eigen::MatrixXd ftetwildMesh_V(ftetwildMesh->vertices.size(), 3);  
    Eigen::MatrixXi ftetwildMesh_F(ftetwildMesh->faces.size(), 3);;  
    for(int i = 0; i < ftetwildMesh->vertices.size(); i++)
    {
        ftetwildMesh_V.row(i) = ftetwildMesh->vertices[i];
    }
    for(int i = 0; i < ftetwildMesh->faces.size(); i++)
    {
        Eigen::Vector3i face = {ftetwildMesh->faces[i][0], ftetwildMesh->faces[i][1],ftetwildMesh->faces[i][2] };
        ftetwildMesh_F.row(i) = face;
    }



    Eigen::MatrixXd V;  
    Eigen::MatrixXi F;  


    { // remove any unreference vertices (if any) from the surface mesh of the tet-mesh
        Eigen::VectorXi _1;
        Eigen::VectorXi _2;
        igl::remove_unreferenced(ftetwildMesh_V, ftetwildMesh_F, V, F, _1, _2);
    }

        {
            // We depend on "fTetWild", which sometimes produces multiple
            // connected components within a single surface mesh. This is
            // problematic when computing geometric quantities on the
            // remeshed input mesh like curviture. Also, MCUT does not take input
            // meshes that contain multiple connected components.
            // Thus, our solution is to keep only the largest connected component.
            Eigen::VectorXi F_cc;
            igl::facet_components(F, F_cc);

            std::map<int, int> CC_to_face_count;
            for (int i = 0; i < F_cc.rows(); ++i)
            {
                int cc_id = F_cc(i);
                if (CC_to_face_count.count(cc_id) == 0)
                {
                    CC_to_face_count[cc_id] = 0; // init
                }

                CC_to_face_count[cc_id] += 1;
            }

            printf("[PROG]: Found %d connected components in input mesh\n", (int)CC_to_face_count.size());

            if (CC_to_face_count.size() > 1)
            {
                int largest_CC = 0;
                int faceCountMax = 0;

                for (auto &e : CC_to_face_count)
                {
                    if (e.second > faceCountMax)
                    {
                        faceCountMax = e.second;
                        largest_CC = e.first;
                    }

                    //printf("[PROG]: Connected component %d has %d faces \n", (int)e.first, (int)e.second);
                }

                {
                    // create a new mesh which contains the faces of the largest connected
                    // component. this new mesh will be the formal "input" mesh that we use.
                    Eigen::MatrixXi newF(CC_to_face_count.at(largest_CC), 3);
                    int face_index = 0;
                    for (int i = 0; i < F_cc.rows(); ++i)
                    {
                        int cc_id = F_cc(i);
                        if (cc_id == largest_CC)
                        {
                            newF.row(face_index) = F.row(i);
                            face_index++;
                        }
                    }

                    { // keep only the vertices that are used by the largest connected component.
                        Eigen::MatrixXd cleanV;
                        Eigen::MatrixXi cleanF;
                        Eigen::VectorXi _1;
                        Eigen::VectorXi _2;
                        igl::remove_unreferenced(V, newF, cleanV, cleanF, _1, _2);
                        V = cleanV;
                        F = cleanF;
                    }
                }
            }

            printf("[PROG]: Final input mesh vertices = %d, faces = %d\n", (int)V.rows(), (int)F.rows());
        }


    for(int i = 0; i < V.rows(); i++)
    {
        Eigen::Vector3d vert = V.row(i);
        outMesh.vertices.push_back(vert);
    }
    for(int i = 0; i < F.rows(); i++)
    {
        std::vector<int> face;
        face.push_back(F(i,0));
        face.push_back(F(i,1));        
        face.push_back(F(i,2));
        outMesh.faces.push_back(face);
    }

    return outMesh;

}


// Struct of grid
struct Grid {
    // property of velocity-field 0
    double m = 0; // each node's mass
    Eigen::Vector3d mom = { 0, 0, 0 }; // each node's momentum
    Eigen::Vector3d velocity = { 0, 0, 0 }; // each node's velocity
    Eigen::Vector3d force = { 0, 0, 0 }; // each node's force

    // general grid node property
    Eigen::Vector3i posIndex = { 0, 0, 0 };
    Eigen::Vector3d deltaDi = { 0, 0, 0 }; // gradient of damage field
    double Di = 0; // value of damage field
    double sw = 0; // sum of particle-grid weight

    // particle index in the support radius of this node. The order of the vector is important
	std::vector<int> supportParticles; // store the position of the particle in vector "particles"; 
	std::vector<double> supportParticlesWeight; // store the weight of particle to the grid node



    // set of crack surface points withing the grid cell
	std::vector<int> crackPoints;
	int nearestPoint = -1000; // (nearestPoint < 0) means it is far away from the crack surface
	Eigen::Vector3d crackSurfaceNormal = { 0,0,0 }; // the vector pointing from the nearest point on the crack surface to the grid node


    // parameters of contact algorithm
	double mass_0 = 0;
	Eigen::Vector3d mom_0 = { 0 , 0 , 0 }; 
	Eigen::Vector3d velocity_0 = { 0 , 0 , 0 };
	Eigen::Vector3d force_0 = { 0 , 0 , 0 };

	double mass_1 = 0;
	Eigen::Vector3d mom_1 = { 0 , 0 , 0 };
	Eigen::Vector3d velocity_1 = { 0 , 0 , 0 };
	Eigen::Vector3d force_1 = { 0 , 0 , 0 }; 



    Grid(double im)
        : m(im)
    {
    }
};

struct chevronMarksGridAndMap
{
    std::vector<Eigen::Vector3d> vertices; // vertices of the crack surface
    std::vector<std::vector<int>> faces; // faces of the crack surface
	std::map<int, int> gridMap;
	std::vector<Grid> nodesVec;
    Eigen::Vector3d minCoordinate = {0,0,0 };
    parametersSim param;
};

////////////////////////////////////////////////////////


// Given the vertices and faces information, write the off file.
static void writeOffFile(std::vector<Eigen::Vector3d> vertices, std::vector<std::vector<int>> faces, std::string name)
{
    std::string path = name + ".off";
    std::cout << "write: " << path << std::endl;
	std::ofstream outfile2(path, std::ios::trunc);
	outfile2 << "OFF" << std::endl;
	outfile2 << vertices.size() << " " << faces.size() << " 0" << std::endl;
	for (int vert = 0; vert < vertices.size(); ++vert)
	{
		outfile2 << std::scientific << std::setprecision(8) << vertices[vert][0] << " " << vertices[vert][1] << " " << vertices[vert][2] << " " << std::endl;
	}

	for (int face = 0; face < faces.size(); ++face)
	{
		outfile2 << faces[face].size() << " ";
		for (int vert = 0; vert < faces[face].size(); ++vert) {
			outfile2 << faces[face][vert] << " ";
		}
		outfile2 << std::endl;
	}


	outfile2.close();

}


// Given the vertices and faces information, write the off file.
static void writeObjFile(std::vector<Eigen::Vector3d> vertices, std::vector<std::vector<int>> faces, std::string name)
{
    std::string path = name + ".obj";
    std::cout << "write: " << path << std::endl;
	std::ofstream outfile2(path, std::ios::trunc);
	for (int vert = 0; vert < vertices.size(); ++vert)
	{
		outfile2 << std::scientific << std::setprecision(8)<<"v " << vertices[vert][0] << " " << vertices[vert][1] << " " << vertices[vert][2] << " " << std::endl;
	}

	for (int face = 0; face < faces.size(); ++face)
	{
		outfile2  << "f ";
		for (int vert = 0; vert < faces[face].size(); ++vert) {
			outfile2 << faces[face][vert] + 1 << " ";
		}
		outfile2 << std::endl;
	}


	outfile2.close();

}
