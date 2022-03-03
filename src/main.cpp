
#include "crackExtraction/object.h"
#include "crackExtraction/extractCrack.h"
#include "crackExtraction/utils.h" 



#include <igl/decimate.h>
#include <thread>





// find the surface mesh from a vdb grid
void getSurfaceMeshFromVdbGrid(openvdb::FloatGrid::Ptr bareMeshVdbGridPtr,int decimateTarget, meshObjFormat& fragmentVolume)
{
    openvdb::tools::VolumeToMesh volumeToMeshHandle;
    volumeToMeshHandle(*bareMeshVdbGridPtr);

    trimesh::TriMesh* pMesh = new trimesh::TriMesh;

    openvdb::tools::PointList *verts = &volumeToMeshHandle.pointList();
    openvdb::tools::PolygonPoolList *polys = &volumeToMeshHandle.polygonPoolList();


    for (size_t i = 0; i < volumeToMeshHandle.pointListSize(); i++)
    {
        openvdb::Vec3s &v = (*verts)[i];
        pMesh->vertices.push_back(trimesh::vec3(v[0], v[1], v[2]));

    }

    for (size_t i = 0; i < volumeToMeshHandle.polygonPoolListSize(); i++)
    {

        for (size_t ndx = 0; ndx < (*polys)[i].numTriangles(); ndx++)
        {
            openvdb::Vec3I *p = &((*polys)[i].triangle(ndx));
        }

        for (size_t ndx = 0; ndx < (*polys)[i].numQuads(); ndx++)
        {
            openvdb::Vec4I *p = &((*polys)[i].quad(ndx));

            trimesh::TriMesh::Face f0;
            f0[0] = p->z();
            f0[1] = p->y();
            f0[2] = p->x();
            trimesh::TriMesh::Face f1;
            f1[0] = p->w();
            f1[1] = p->z();
            f1[2] = p->x();

            pMesh->faces.push_back(f0);
            pMesh->faces.push_back(f1);

        }
    }


    // decimate mesh using libigl
    Eigen::MatrixXd V(pMesh->vertices.size(), 3);
    Eigen::MatrixXi F(pMesh->faces.size(), 3);
    for(int i = 0; i < pMesh->vertices.size(); i++)
    {
        V(i,0) = pMesh->vertices[i][0];
        V(i,1) = pMesh->vertices[i][1];
        V(i,2) = pMesh->vertices[i][2];
    }
    for(int i = 0; i < pMesh->faces.size(); i++)
    {
        F(i,0) = pMesh->faces[i][0];
        F(i,1) = pMesh->faces[i][1];
        F(i,2) = pMesh->faces[i][2];
    }



    Eigen::MatrixXd U;
    Eigen::MatrixXi G;
    Eigen::VectorXi J;
    Eigen::VectorXi I;
    igl::decimate(V, F, decimateTarget, U, G, J, I);
    for (int i = 0; i < U.rows(); ++i)
    {
        Eigen::Vector3d p = {U(i, 0), U(i, 1), U(i, 2)};
        fragmentVolume.vertices.push_back(p);
    }
    for (int i = 0; i < G.rows(); ++i)
    {
        std::vector<int> f;
        f.push_back(G(i, 0));
        f.push_back(G(i, 1));
        f.push_back(G(i, 2));
        fragmentVolume.faces.push_back(f);
    }




}


// the input particles may locate in the negative domain. This function projects them to the positive domain
void preprocessing(std::string crackFilePath, std::string cutObjectFilePath, parametersSim* parameters, std::vector<Particle>* particleVec, meshObjFormat* objectMesh)
{
    // project damaged particles to the positive domain
    // read damaged particles
    std::ifstream inf;
    inf.open(crackFilePath);
    std::string sline;
    std::string s0, s1, s2, s3;
    while (getline(inf, sline)) {

        std::istringstream sin(sline);
        sin >> s0 >> s1 >> s2 >> s3;
        Eigen::Vector3d ipos = { atof(s0.c_str()), atof(s1.c_str()), atof(s2.c_str()) };
        Eigen::Vector3d ivel = { 0, 0, 0 };
        double iDp = atof(s3.c_str());
        (*particleVec).push_back(Particle(ipos, ivel, 0, 0, iDp));
		
    }
 

    
    double xmin = (*particleVec)[0].pos[0], xmax = (*particleVec)[0].pos[0], ymin = (*particleVec)[0].pos[1], ymax = (*particleVec)[0].pos[1], zmin = (*particleVec)[0].pos[2], zmax = (*particleVec)[0].pos[2];
    for (int km = 0; km < (*particleVec).size(); km++) 
    {
        Eigen::Vector3d ipos = (*particleVec)[km].pos;
        xmin = std::min(xmin, ipos[0]);
        xmax = std::max(xmax, ipos[0]);
        ymin = std::min(ymin, ipos[1]);
        ymax = std::max(ymax, ipos[1]);
        zmin = std::min(zmin, ipos[2]);
        zmax = std::max(zmax, ipos[2]);
    }
    Eigen::Vector3d minCoordinate = { xmin - 10 * (*parameters).dx, ymin - 10 * (*parameters).dx, zmin - 10 * (*parameters).dx };
    for (int km = 0; km < (*particleVec).size(); km++) {
        (*particleVec)[km].pos = (*particleVec)[km].pos - minCoordinate;
    }
    Eigen::Vector3d length = {xmax - xmin + 20* (*parameters).dx,  ymax - ymin + 20* (*parameters).dx, zmax - zmin + 20 * (*parameters).dx};
    (*parameters).length = length;
    (*parameters).minCoordinate = minCoordinate;



    Eigen::MatrixXd  V;
    Eigen::MatrixXi  F;
    bool success = igl::readOBJ(cutObjectFilePath,V,F);
    if(!success)
    {
        std::cout<<"Fail to load object mesh"<<std::endl;
    }
    else
    {

        for (int i = 0; i < V.rows(); ++i)
        {
            Eigen::Vector3d pos = {V(i,0), V(i,1),V(i,2)};
            (*objectMesh).vertices.push_back(pos - minCoordinate);
        }

        for (int i = 0; i < F.rows(); ++i)
        {
            std::vector<int> face;
            for (int j = 0; j < F.cols(); ++j)
            {
                face.push_back(F(i,j));
            }
            (*objectMesh).faces.push_back(face);
        }
    }
   
}


// the input particles may locate in the negative domain. This function projects them to the positive domain
void postprocessing(parametersSim* parameters, std::tuple<bool, meshObjFormat,  meshObjFormat, std::vector<meshObjFormat> >* crackSurface, meshObjFormat* objectMesh)
{
    // project damaged particles to the original domain
    for(int m = 0; m < std::get<1>(*crackSurface).vertices.size(); m++)
    {
        std::get<1>(*crackSurface).vertices[m] = std::get<1>(*crackSurface).vertices[m] + (*parameters).minCoordinate;
    }


    for(int m = 0; m < std::get<2>(*crackSurface).vertices.size(); m++)
    {
        std::get<2>(*crackSurface).vertices[m] = std::get<2>(*crackSurface).vertices[m] + (*parameters).minCoordinate;
    }


    for(int m = 0; m < std::get<3>(*crackSurface).size(); m++)
    {
        for(int n = 0; n < std::get<3>(*crackSurface)[m].vertices.size(); n++)
        {
            std::get<3>(*crackSurface)[m].vertices[n] = std::get<3>(*crackSurface)[m].vertices[n] + (*parameters).minCoordinate;
        }
    }


    // project cutting objects to the original domain
    for(int m = 0; m < (*objectMesh).vertices.size(); m++)
    {
        (*objectMesh).vertices[m] = (*objectMesh).vertices[m] + (*parameters).minCoordinate;
    }

}




int main(int argc, char *argv[])
{
    std::vector<std::string> inputPara;
    std::string inputFile = argv[1];
    std::ifstream inf;
    inf.open(inputFile);
    std::string sline;
    while (getline(inf, sline)) 
    {
        inputPara.push_back(sline);		
    }
    inf.close();
 


    // the resolution of crack faces
    // parameters
    parametersSim parameters;
    parameters.dx = std::stod(inputPara[4]); 
    parameters.vdbVoxelSize = std::stod(inputPara[5]);
    // std::string crackFilePath = "/home/floyd/Linxu/crackExtraction/crackExtraction/build/output/particles.txt";
    // std::string cutObjectFilePath = "/home/floyd/Linxu/crackExtraction/crackExtraction/build/output/object.obj";
    std::string crackFilePath = inputPara[0] + "/" + inputPara[1];
    std::string cutObjectFilePath = inputPara[0] + "/" + inputPara[2];
    std::vector<Particle> particleVec;
    meshObjFormat objectMesh;


    // extract the crack surface
    preprocessing(crackFilePath, cutObjectFilePath, &parameters, &particleVec, &objectMesh); 
    
    std::tuple<bool, meshObjFormat,  meshObjFormat, std::vector<meshObjFormat> > result =  extractCrackSurface(&particleVec, parameters);
    
    
    if(std::get<0>(result) == false)
    {
        std::cout<<"No crack found!"<<std::endl;
    }
    else
    {
        postprocessing(&parameters, &result, &objectMesh);

        // output the crack surface and fragments
        meshObjFormat crackSurfacePartialCut = std::get<1>(result);
        writeObjFile(crackSurfacePartialCut.vertices, crackSurfacePartialCut.faces, inputPara[0] + "/" + "partialCutSurface");
        meshObjFormat crackSurfaceFullCut = std::get<2>(result);
        writeObjFile(crackSurfaceFullCut.vertices, crackSurfaceFullCut.faces, inputPara[0] + "/" + "fullCutSurface");
        std::vector<meshObjFormat> fragments = std::get<3>(result);

        
        // define the cutting method
        // case 0: complete cut with MCUT
        // case 1: complete cut with openVDB
        // case 2: partial cut with openVDB
        std::string cuttingMethod = inputPara[3];
        if(cuttingMethod == "MCUT")
        {
            std::vector<meshObjFormat> fragmentsFinal;
            cutObject_MCUT(parameters,"tmpCutPbject", &objectMesh, &fragments, &fragmentsFinal);

            // for each fragment level
            for (unsigned int i = 0; i < fragmentsFinal.size(); ++i) 
            {
                writeObjFile(fragmentsFinal[i].vertices, fragmentsFinal[i].faces, inputPara[0] + "/" + "fullCut_MCUT_Fragment_"+std::to_string(i));
            }
        }
        else if(cuttingMethod == "OPENVDB_FULL")
        {

        
            // define openvdb linear transformation
            openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(parameters.vdbVoxelSize);

            
            // convert crack surface mesh to vdb grid
            GenericMesh crackSurfaceGeneric;
            crackSurfaceGeneric.m.vertices = crackSurfaceFullCut.vertices;
            crackSurfaceGeneric.m.faces = crackSurfaceFullCut.faces;
            triangulateGenericMesh(crackSurfaceGeneric);


            std::vector<openvdb::Vec3f> myMeshPoints_crack;
            for (int i = 0; i < crackSurfaceGeneric.m.vertices.size(); ++i) 
            {
                const Eigen::Vector3d p = crackSurfaceGeneric.m.vertices[i];
                myMeshPoints_crack.push_back(openvdb::Vec3f(p[0], p[1], p[2]));
            }
            std::vector<openvdb::Vec3I> myMeshTris_crack;
            for (int i = 0; i < crackSurfaceGeneric.triangulatedFaces.size() / 3; ++i) 
            {
                myMeshTris_crack.push_back(openvdb::Vec3I(crackSurfaceGeneric.triangulatedFaces[3 * i + 0], crackSurfaceGeneric.triangulatedFaces[3 * i + 1], crackSurfaceGeneric.triangulatedFaces[3 * i + 2]));
            }



            openvdb::FloatGrid::Ptr crackLevelSetGrid = openvdb::tools::meshToUnsignedDistanceField<openvdb::FloatGrid>(
                *transform,
                myMeshPoints_crack,
                myMeshTris_crack,
                std::vector<openvdb::Vec4I>(),
                3);

            for (openvdb::FloatGrid::ValueOnIter iter = crackLevelSetGrid->beginValueOn(); iter; ++iter) {
                float dist = iter.getValue();
                float value = dist - std::sqrt(3 * std::pow(parameters.vdbVoxelSize, 2));
                iter.setValue(value);
            }
            crackLevelSetGrid->setGridClass(openvdb::GRID_LEVEL_SET);




            // convert object mesh to vdb grid
            GenericMesh objectMeshGeneric;
            objectMeshGeneric.m.vertices = objectMesh.vertices;
            objectMeshGeneric.m.faces = objectMesh.faces;
            triangulateGenericMesh(objectMeshGeneric);


            std::vector<openvdb::Vec3f> myMeshPoints_object;
            for (int i = 0; i < objectMeshGeneric.m.vertices.size(); ++i) 
            {
                const Eigen::Vector3d p = objectMeshGeneric.m.vertices[i];
                myMeshPoints_object.push_back(openvdb::Vec3f(p[0], p[1], p[2]));
            }
            std::vector<openvdb::Vec3I> myMeshTris_object;
            for (int i = 0; i < objectMeshGeneric.triangulatedFaces.size() / 3; ++i) 
            {
                myMeshTris_object.push_back(openvdb::Vec3I(objectMeshGeneric.triangulatedFaces[3 * i + 0], objectMeshGeneric.triangulatedFaces[3 * i + 1], objectMeshGeneric.triangulatedFaces[3 * i + 2]));
            }


            // create a float grid containing the level representation of the sphere mesh
            openvdb::FloatGrid::Ptr objectLevelSetGrid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
                *transform,
                myMeshPoints_object,
                myMeshTris_object);
            objectLevelSetGrid->setGridClass(openvdb::GRID_LEVEL_SET);



            // do the boolean operation
            openvdb::FloatGrid::Ptr copyOfCrackGrid = crackLevelSetGrid->deepCopy();
            openvdb::FloatGrid::Ptr copyOfObjGrid = objectLevelSetGrid->deepCopy();
            // Compute the difference (A / B) of the two level sets.
            openvdb::tools::csgDifference(*copyOfObjGrid, *copyOfCrackGrid);
            openvdb::FloatGrid::Ptr csgSubtractedObjGrid = copyOfObjGrid; // cutted piece
            // list of fragment pieces after cutting (as level set grids)
            std::vector<openvdb::FloatGrid::Ptr> fragmentLevelSetGridPtrList;
            openvdb::tools::segmentSDF(*csgSubtractedObjGrid, fragmentLevelSetGridPtrList);




            // for each fragment level
            for (unsigned int i = 0; i < fragmentLevelSetGridPtrList.size(); ++i) 
            {
                meshObjFormat fullCutFragment;
                getSurfaceMeshFromVdbGrid(fragmentLevelSetGridPtrList[i],1000000, fullCutFragment);
                writeObjFile(fullCutFragment.vertices, fullCutFragment.faces, inputPara[0] + "/" + "fullCutFragment_"+std::to_string(i));
            }

            
           
        }
        else if(cuttingMethod == "OPENVDB_PARTIAL")
        {

            // define openvdb linear transformation
            openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(parameters.vdbVoxelSize);

            
            // convert crack surface mesh to vdb grid
            GenericMesh crackSurfaceGeneric;
            crackSurfaceGeneric.m.vertices = crackSurfacePartialCut.vertices;
            crackSurfaceGeneric.m.faces = crackSurfacePartialCut.faces;
            triangulateGenericMesh(crackSurfaceGeneric);


            std::vector<openvdb::Vec3f> myMeshPoints_crack;
            for (int i = 0; i < crackSurfaceGeneric.m.vertices.size(); ++i) 
            {
                const Eigen::Vector3d p = crackSurfaceGeneric.m.vertices[i];
                myMeshPoints_crack.push_back(openvdb::Vec3f(p[0], p[1], p[2]));
            }
            std::vector<openvdb::Vec3I> myMeshTris_crack;
            for (int i = 0; i < crackSurfaceGeneric.triangulatedFaces.size() / 3; ++i) 
            {
                myMeshTris_crack.push_back(openvdb::Vec3I(crackSurfaceGeneric.triangulatedFaces[3 * i + 0], crackSurfaceGeneric.triangulatedFaces[3 * i + 1], crackSurfaceGeneric.triangulatedFaces[3 * i + 2]));
            }



            openvdb::FloatGrid::Ptr crackLevelSetGrid = openvdb::tools::meshToUnsignedDistanceField<openvdb::FloatGrid>(
                *transform,
                myMeshPoints_crack,
                myMeshTris_crack,
                std::vector<openvdb::Vec4I>(),
                3);

            for (openvdb::FloatGrid::ValueOnIter iter = crackLevelSetGrid->beginValueOn(); iter; ++iter) {
                float dist = iter.getValue();
                float value = dist - std::sqrt(3 * std::pow(parameters.vdbVoxelSize, 2));
                iter.setValue(value);
            }
            crackLevelSetGrid->setGridClass(openvdb::GRID_LEVEL_SET);




            // convert object mesh to vdb grid
            GenericMesh objectMeshGeneric;
            objectMeshGeneric.m.vertices = objectMesh.vertices;
            objectMeshGeneric.m.faces = objectMesh.faces;
            triangulateGenericMesh(objectMeshGeneric);


            std::vector<openvdb::Vec3f> myMeshPoints_object;
            for (int i = 0; i < objectMeshGeneric.m.vertices.size(); ++i) 
            {
                const Eigen::Vector3d p = objectMeshGeneric.m.vertices[i];
                myMeshPoints_object.push_back(openvdb::Vec3f(p[0], p[1], p[2]));
            }
            std::vector<openvdb::Vec3I> myMeshTris_object;
            for (int i = 0; i < objectMeshGeneric.triangulatedFaces.size() / 3; ++i) 
            {
                myMeshTris_object.push_back(openvdb::Vec3I(objectMeshGeneric.triangulatedFaces[3 * i + 0], objectMeshGeneric.triangulatedFaces[3 * i + 1], objectMeshGeneric.triangulatedFaces[3 * i + 2]));
            }


            // create a float grid containing the level representation of the sphere mesh
            openvdb::FloatGrid::Ptr objectLevelSetGrid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
                *transform,
                myMeshPoints_object,
                myMeshTris_object);
            objectLevelSetGrid->setGridClass(openvdb::GRID_LEVEL_SET);



            // do the boolean operation
            openvdb::FloatGrid::Ptr copyOfCrackGrid = crackLevelSetGrid->deepCopy();
            openvdb::FloatGrid::Ptr copyOfObjGrid = objectLevelSetGrid->deepCopy();
            // Compute the difference (A / B) of the two level sets.
            openvdb::tools::csgDifference(*copyOfObjGrid, *copyOfCrackGrid);
            openvdb::FloatGrid::Ptr csgSubtractedObjGrid = copyOfObjGrid; // cutted piece
            // list of fragment pieces after cutting (as level set grids)
            std::vector<openvdb::FloatGrid::Ptr> fragmentLevelSetGridPtrList;
            openvdb::tools::segmentSDF(*csgSubtractedObjGrid, fragmentLevelSetGridPtrList);




            // for each fragment level
            for (unsigned int i = 0; i < fragmentLevelSetGridPtrList.size(); ++i) 
            {
                meshObjFormat fullCutFragment;
                getSurfaceMeshFromVdbGrid(fragmentLevelSetGridPtrList[i],1000000, fullCutFragment);
                writeObjFile(fullCutFragment.vertices, fullCutFragment.faces, inputPara[0] + "/" + "partialCutFragment_"+std::to_string(i));
            }
        
        }
        else
        {
            std::cout<<"The cutting method is incorrect!"<<std::endl;
        }
   

        
        

    }

    
    
    return 0;
}
