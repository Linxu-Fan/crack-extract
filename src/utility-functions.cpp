#include "mpm-fracture/utils.h"
#include <mapbox/earcut.hpp>

void toTrimesh(trimesh::TriMesh& tm, const GenericMesh& gm)
{
    ASSERT(gm.m.vertices.size() >= 3);
    for (int i = 0; i < gm.m.vertices.size(); ++i) {
        const Eigen::Vector3d& v = gm.m.vertices[i];
        tm.vertices.push_back(trimesh::vec3(v[0], v[1], v[2]));
    }

    ASSERT(gm.triangulatedFaces.size() >= 3);

    for (int i = 0; i < gm.triangulatedFaces.size() / 3; ++i) {
        trimesh::TriMesh::Face f(gm.triangulatedFaces[i * 3 + 0], gm.triangulatedFaces[i * 3 + 1], gm.triangulatedFaces[i * 3 + 2]);
        tm.faces.push_back(f);
    }
}



// convert a mesh with chevron marks to openvdb grid
openvdb::FloatGrid::Ptr chevronMeshToLevelSetGrid(
    trimesh::TriMesh& chevronMesh,
    const float voxelSize,
    Eigen::Vector3f com,
    const std::string name)
{
    // copy vertices and triangles

    std::vector<openvdb::Vec3f> myMeshPoints;
    for (int i = 0; i < chevronMesh.vertices.size(); ++i) 
    {
        myMeshPoints.push_back(openvdb::Vec3f(chevronMesh.vertices[i][0] + com[0], chevronMesh.vertices[i][1] + com[1], chevronMesh.vertices[i][2] + com[2]));
    }

    std::vector<openvdb::Vec3I> myMeshTris;
    for (int i = 0; i < chevronMesh.faces.size(); ++i) 
    {
        myMeshTris.push_back(openvdb::Vec3I(chevronMesh.faces[i][0], chevronMesh.faces[i][1], chevronMesh.faces[i][2]));
    }

    // convert crack surface mesh to level set

    openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(voxelSize);


    // create a float grid containing the level representation of the sphere mesh
    openvdb::FloatGrid::Ptr myCrackMeshLevelSetGrid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
        *transform,
        myMeshPoints,
        myMeshTris);

    myCrackMeshLevelSetGrid->setGridClass(openvdb::GRID_LEVEL_SET);
    myCrackMeshLevelSetGrid->setName(name);

    if (myCrackMeshLevelSetGrid->empty()) {
        fprintf(stderr, "fatal error: grid empty (%s)\n", name.c_str());
        std::exit(1);
    }


    return myCrackMeshLevelSetGrid;
}



openvdb::FloatGrid::Ptr crackSurfaceMeshToLevelSetGrid(
    const GenericMesh& crackSurfaceMesh,
    const float voxelSize,
    const std::string& name)
{
    // the size of the narrowband of the level set defining an unsigned distance field
    // around the crack surface.
    //float bandwidth = voxelSize * crackMeshLevelSetBandwidthScalingFactor;
    //printf("narrowband bandwidth=%f\n", bandwidth);

    // copy vertices and triangles

    std::vector<openvdb::Vec3f> myMeshPoints;
    ASSERT(crackSurfaceMesh.m.vertices.empty() == false);
    for (int i = 0; i < crackSurfaceMesh.m.vertices.size(); ++i) {
        const Eigen::Vector3d p = crackSurfaceMesh.m.vertices[i];
        myMeshPoints.push_back(openvdb::Vec3f(p[0], p[1], p[2]));
    }

    std::vector<openvdb::Vec3I> myMeshTris;
    ASSERT(crackSurfaceMesh.triangulatedFaces.empty() == false);

    for (int i = 0; i < crackSurfaceMesh.triangulatedFaces.size() / 3; ++i) {
        const int* f = crackSurfaceMesh.triangulatedFaces.data() + (i * 3);
        myMeshTris.push_back(openvdb::Vec3I(f[0], f[1], f[2]));
    }

    // convert crack surface mesh to level set

    openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(voxelSize);

#ifdef USE_LEVEL_SET_FULLYCUT

    // create a float grid containing the level representation of the sphere mesh
    openvdb::FloatGrid::Ptr myCrackMeshLevelSetGrid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
        *transform,
        myMeshPoints,
        myMeshTris);

#else 

    //note that we first create an [unsigned] level set here
    openvdb::FloatGrid::Ptr myCrackMeshLevelSetGrid = openvdb::tools::meshToUnsignedDistanceField<openvdb::FloatGrid>(
        *transform,
        myMeshPoints,
        myMeshTris,
        std::vector<openvdb::Vec4I>(),
        3);


#endif








    if (myCrackMeshLevelSetGrid->empty()) {
        fprintf(stderr, "fatal error: grid empty (%s)\n", name.c_str());
        std::exit(1);
    }

    std::string fullName = name + "_CrackSurfaceVolGrid";
    myCrackMeshLevelSetGrid->setName(fullName);





#ifdef USE_LEVEL_SET_FULLYCUT

    std::cout<<""<<std::endl;

#else 

    // now, we create a [signed] level set which leads to our crack surface volume.

    // Visit and update all of the grid's active values, which correspond to
    // voxels on the narrow band.
    for (openvdb::FloatGrid::ValueOnIter iter = myCrackMeshLevelSetGrid->beginValueOn(); iter; ++iter) {
        float dist = iter.getValue();

        float value = dist - std::sqrt(3 * std::pow(voxelSize, 2));
        iter.setValue(value);
    }


#endif








    myCrackMeshLevelSetGrid->setGridClass(openvdb::GRID_LEVEL_SET);





    return myCrackMeshLevelSetGrid;
}

// only outputs .obj file
void writeVDBMesh(const char* filename, openvdb::tools::VolumeToMesh& mesh)
{
    std::cout << "write: " << filename << std::endl;
    std::ofstream file;
    file.open(filename);

    openvdb::tools::PointList* verts = &mesh.pointList();
    openvdb::tools::PolygonPoolList* polys = &mesh.polygonPoolList();

    for (size_t i = 0; i < mesh.pointListSize(); i++) {
        openvdb::Vec3s& v = (*verts)[i];
        file << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    for (size_t i = 0; i < mesh.polygonPoolListSize(); i++) {

        for (size_t ndx = 0; ndx < (*polys)[i].numTriangles(); ndx++) {
            openvdb::Vec3I* p = &((*polys)[i].triangle(ndx));
            file << "f " << p->x() << " " << p->y() << " " << p->z() << std::endl;
        }

        for (size_t ndx = 0; ndx < (*polys)[i].numQuads(); ndx++) {
            openvdb::Vec4I* p = &((*polys)[i].quad(ndx));
            file << "f " << p->z() + 1 << " " << p->y() + 1 << " " << p->x() + 1 << std::endl;
            file << "f " << p->w() + 1 << " " << p->z() + 1 << " " << p->x() + 1 << std::endl;
        }
    }

    file.close();
}

void saveOpenVDBGrid(const openvdb::FloatGrid::Ptr& grid, const std::string& outDir)
{
#if defined(DUMP_DEBUG_MESH_DATA_TO_FILE)
    std::string name = grid->getName();
    std::string fname = pystring::os::path::join({ outDir, "_" + name + ".vdb" }); // dash needed incase "name" is empty
    printf("%% save grid %s\n", fname.c_str());
    openvdb::io::File outfile(fname);
    outfile.write({ grid });
    outfile.close();
#else
    printf("%% skip logging vdb grid to file\n");
#endif
}

void saveOpenVDBGrids(const openvdb::GridPtrVec& grids, const std::string& outDir, const std::string& name)
{
#if defined(DUMP_DEBUG_MESH_DATA_TO_FILE)
    std::string fname = pystring::os::path::join({ outDir, "_" + name + ".vdb" }); // dash needed incase "name" is empty
    printf("%% save grids %s\n", fname.c_str());
    openvdb::io::File outfile(fname);
    outfile.write(grids);
    outfile.close();
#else
    printf("%% skip logging vdb grid to file\n");
#endif
}

// watertight mesh
openvdb::FloatGrid::Ptr meshToVDBLevelSetGrid(
    const trimesh::TriMesh* mesh,
    const float voxelSize,
    const std::string& name,
    float& volume)
{
    std::vector<openvdb::Vec3s> vdbMeshPoints;
    std::vector<openvdb::Vec3I> vdbMeshTriangles;

    vdbMeshPoints.resize(mesh->vertices.size());
    for (int i = 0; i < mesh->vertices.size(); ++i) {
        const trimesh::vec3& v = mesh->vertices[i];
        vdbMeshPoints[i] = openvdb::Vec3s(v[0], v[1], v[2]);
    }

    vdbMeshTriangles.resize(mesh->faces.size());
    for (int i = 0; i < mesh->faces.size(); ++i) {
        const trimesh::TriMesh::Face& f = mesh->faces[i];
        vdbMeshTriangles[i] = openvdb::Vec3I(f[0], f[1], f[2]);
    }

    // Associate a scaling transform with the grid that sets the voxel size
    // to "voxelSize" units in world space .
    openvdb::math::Transform::Ptr linearTransf = openvdb::math::Transform::createLinearTransform(voxelSize);

    // create a float grid containing the level representation of the sphere mesh
    openvdb::FloatGrid::Ptr levelSetGridPtr = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
        *linearTransf,
        vdbMeshPoints,
        vdbMeshTriangles);

    levelSetGridPtr->setGridClass(openvdb::GRID_LEVEL_SET);
    levelSetGridPtr->setName(name);

    ASSERT(levelSetGridPtr->hasUniformVoxels());

    // If we need to know the the voxel size is, we need to inspect the transform
    // of this grid. If the transform has uniform scaling (use grid.hasUniformVoxels())
    // then the voxel size is simply by grid.voxelSize()[0].

    ASSERT(levelSetGridPtr->voxelSize()[0] == voxelSize)

    volume = openvdb::tools::levelSetVolume(*levelSetGridPtr);

    return levelSetGridPtr;
}

// remove duplicate vertices and update connectivity
void cleanGenericMesh(GenericMesh& gm)
{
    if (gm.m.vertices.empty() || gm.m.faces.empty()) {
        printf("mesh not cleaned. skip.\n");
        return;
    }

    // check if zero-based indexing

    // for each face
    int minidx = 9999;
    for (int i = 0; i < gm.m.faces.size(); ++i) {
        std::vector<int>& face = gm.m.faces[i];
        for (int j = 0; j < face.size(); ++j) {
            minidx = std::min(face[j], minidx);
        }
    }

    ASSERT(minidx == 0 || minidx == 1);
    bool usingAbsIndexing = (minidx == 1);

    if (usingAbsIndexing) {
        // change to zero-based indexing
        for (int i = 0; i < gm.m.faces.size(); ++i) { // for each face
            std::vector<int>& face = gm.m.faces[i];

            for (int j = 0; j < face.size(); ++j) { // for each vertex in face
                face[j] -= 1;
            }
        }
    }

    // eliminate duplicate vertices in each face
    for (int i = 0; i < gm.m.faces.size(); ++i) { // for each face
        std::vector<int>& face = gm.m.faces[i];
        std::vector<int> faceWithoutDuplicates; //
        std::vector<int> duplicatesInFace;
        int N = face.size();
        for (int j = 0; j < N; ++j) { // for each vertex in face
            int idxCurrent = face[j];
            const Eigen::Vector3d& vertex = gm.m.vertices[idxCurrent];
            bool currentIsDuplicate = false;
            /// compare "vertex" with all others in "face"

            for (int k = 1; k < N; ++k) { // for each "other" vertex
                int idxOther = face[(j + k) % N]; // wrap index
                const Eigen::Vector3d& otherVertex = gm.m.vertices[idxOther];
                bool isDuplicate = std::fabs(vertex[0] - otherVertex[0]) < 1e-8 && std::fabs(vertex[1] - otherVertex[1]) < 1e-8 && std::fabs(vertex[2] - otherVertex[2]) < 1e-8;
                if (isDuplicate) {
                    currentIsDuplicate = true; // at least once
                    duplicatesInFace.push_back(idxOther);
                }
            }

            // if the current is a duplicate, we want to know if it is the first in the sequence (i.e. it comes first, so that we keep it)
            bool isRegisteredDuplicate = std::find(duplicatesInFace.begin(), duplicatesInFace.end(), idxCurrent) != duplicatesInFace.end();
            bool keep = ((currentIsDuplicate && isRegisteredDuplicate) == false);
            if (keep) {
                faceWithoutDuplicates.push_back(idxCurrent);
            }
        }
        face = faceWithoutDuplicates;
    }
}

std::vector<std::vector<std::array<double, 2>>> to2Dpoly(std::vector<Eigen::Vector3d>& poly3d)
{
    const Eigen::Vector3d& v0 = poly3d[0];

    int num_verts = (int)poly3d.size();
    Eigen::Vector3d normal(0.0, 0.0, 0.0); // polygon normal
    // longest computation vector in-plane (using available vertices)
    Eigen::Vector3d longestVectorInPlane(0);

    // calculate normal of polygon
    for (int i = 1; i < num_verts; ++i) {
        const Eigen::Vector3d diff = poly3d[i] - v0;
        normal = normal + diff.cross(poly3d[(i + 1) % num_verts] - v0);

        if (diff.norm() > longestVectorInPlane.norm()) {
            longestVectorInPlane = diff;
        }
    }

    // https://answers.unity.com/questions/1522620/converting-a-3d-polygon-into-a-2d-polygon.html

    // first unit vector on the plane
    Eigen::Vector3d u = longestVectorInPlane.normalized();

    // second unit vector just calculate
    Eigen::Vector3d v = (u.cross(normal)).normalized();

    // the 2d polygon
    std::vector<std::vector<std::array<double, 2>>> poly2d;
    poly2d.resize(1); // just one polygon (stuck with type required by earcut lib)
    poly2d.back().resize(poly3d.size()); // same number of vertices

    for (int i = 0; i < (int)poly3d.size(); ++i) {
        const Eigen::Vector3d& point = poly3d[i];
        Eigen::Vector2d p(point.dot(u), point.dot(v)); // 2d point
        poly2d[0][i] = { p[0], p[1] };
    }

    return poly2d;
}

#include "mpm-fracture/polygon_triangulate.hpp"

void triangulateGenericMesh(GenericMesh& m)
{
    if (m.m.vertices.empty() || m.m.faces.empty()) {
        std::fprintf(stderr, "mesh is empty. skipping triangulation\n");
        return;
    }

    std::vector<int>& triangles = m.triangulatedFaces;
    triangles.clear();
    triangles.reserve(m.m.faces.size() / 3);

    //
    //  Triangulate
    // https://github.com/mapbox/earcut.hpp/blob/master/include/mapbox/earcut.hpp

    for (int i = 0; i < (int)m.m.faces.size(); ++i) { // for each face
        // std::cout << "face=" << i << std::endl;
        const std::vector<int>& face = m.m.faces[i];
        const int numVertsInFace = (int)face.size();

        if (numVertsInFace > 3) {

            std::vector<Eigen::Vector3d> poly3d; // vertices of current face
            std::map<int, uint32_t> imap; // index map
            //std::cout << "numVertsInFace=" << numVertsInFace << std::endl;
            //std::set<uint32_t> indicesSetTest;
            for (int j = 0; j < numVertsInFace; ++j) {
                uint32_t idx = m.m.faces[i][j]; // -1;
                //indicesSetTest.insert(m.m.faces[i][j]);
                const Eigen::Vector3d p = m.m.vertices[idx];
                imap[j] = idx; // save index in "meshVertices" (need for remapping later, below)
                //std::cout << p[0] << " " << p[1] << " " << p[2] << " " << std::endl;
                poly3d.push_back(p);
            }
            // ASSERT(indicesSetTest.size() == face.size());
            // std::cout << std::endl;
            // convert our 3d polygon to 2d
            std::vector<std::vector<std::array<double, 2>>> poly2D = to2Dpoly(poly3d);

            // Run tessellation using earcut lib
            // Three subsequent indices form a triangle. Output triangles are clockwise.
            //
//    There are N-3 triangles in the triangulation.
//
//    For the first N-2 triangles, the first edge listed is always an
//    internal diagonal.

            std::vector<uint32_t> indices = mapbox::earcut<uint32_t>(poly2D);

            std::set<uint32_t> indicesSet;
            for (int j = 0; j < (int)indices.size(); ++j) {
                ASSERT(imap.find(indices[j]) != imap.cend());
                indices[j] = imap[indices[j]]; // from local to global
                indicesSet.insert(indices[j]);
            }

            ASSERT(indicesSet.size() == imap.size());

            // NOTE: inserting indices in reverse order because "mapbox" winding order is clockwise
            triangles.insert(triangles.end(), indices.crbegin(), indices.crend());
        } else {

            int idx0 = m.m.faces[i][0];
            int idx1 = m.m.faces[i][1];
            int idx2 = m.m.faces[i][2];

            triangles.push_back(idx0);
            triangles.push_back(idx1);
            triangles.push_back(idx2);
        }
    }
}
