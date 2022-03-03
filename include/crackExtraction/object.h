#pragma once

#include "crackExtraction/utils.h"
#include "crackExtraction/weights.h"

#include <igl/AABB.h>
#include <Eigen/Sparse>
#include <igl/facet_components.h>
#include <igl/remove_unreferenced.h>
#include <igl/gaussian_curvature.h>



#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>    // level set CSG operations
#include <openvdb/tools/LevelSetUtil.h> // segmentSDF(...)

#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/decimate.h>
#include <igl/edge_flaps.h>
#include <igl/parallel_for.h>
#include <igl/qslim.h>
#include <igl/shortest_edge_and_midpoint.h>


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


void toLibiglMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, trimesh::TriMesh &tm)
{
    V.resize(tm.vertices.size(), 3);
    for (int i = 0; i < tm.vertices.size(); ++i)
    {
        const trimesh::point &p = tm.vertices[i];
        V(i, 0) = p[0];
        V(i, 1) = p[1];
        V(i, 2) = p[2];
    }

    F.resize(tm.faces.size(), 3);
    for (int i = 0; i < tm.faces.size(); ++i)
    {
        const trimesh::TriMesh::Face &f = tm.faces[i];
        F(i, 0) = f[0];
        F(i, 1) = f[1];
        F(i, 2) = f[2];
    }
}


void toTrimesh(trimesh::TriMesh &tm, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    tm.vertices.clear();
    for (int i = 0; i < V.rows(); ++i)
    {
        trimesh::point3 p(V(i, 0), V(i, 1), V(i, 2));
        tm.vertices.push_back(p);
    }

    tm.faces.clear();
    for (int i = 0; i < F.rows(); ++i)
    {
        trimesh::TriMesh::Face f(F(i, 0), F(i, 1), F(i, 2));
        tm.faces.push_back(f);
    }
}

void reduceTriangles(trimesh::TriMesh *&mesh, int target)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    toLibiglMesh(V, F, *mesh);

    Eigen::MatrixXd U;
    Eigen::MatrixXi G;
    Eigen::VectorXi J;
    Eigen::VectorXi I;
    printf("remesh target: %d\nrun igl::decimate\n", target);
    //igl::qslim(V, F, (std::size_t)target, U, G, J, I);
    igl::decimate(V, F, (std::size_t)target, U, G, J, I);

    toTrimesh(*mesh, U, G);
}

