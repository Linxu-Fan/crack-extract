// #include "mpm-fracture/extractCrack.h"
// #include "mpm-fracture/mpm-fracture.h"
// #include "mpm-fracture/object.h"
// #include "mpm-fracture/params.h"
// #include "mpm-fracture/particles.h"
// #include "mpm-fracture/utils.h"

// #include <igl/decimate.h>
// #include <openvdb/points/PointCount.h>
// #include <openvdb/points/PointDataGrid.h>
// #include <openvdb/points/PointScatter.h>





// void generateMaterialParticles(std::vector<Eigen::Vector3d>& mpmParticlePositions, Object* pBrb)
// {

//     SCOPED_TIMER(__FUNCTION__);

//     openvdb::FloatGrid::Ptr objectGridCpy = pBrb->levelSetGridPtr->deepCopy();


//     // Ensure all tiles have been voxelized
//     // Propagate the outside/inside sign information from the narrow band
//     // throughout the grid.
//     openvdb::tools::signedFloodFill(objectGridCpy->tree());

//     for (auto iter = objectGridCpy->beginValueOff(); iter; ++iter) {
//         if (iter.getValue() < 0.0) {
//             iter.setActiveState(true);
//         }
//     }

//     // NOTE: we assume that the grid class is already set based on the naming of "meshToLevelSet"
//     objectGridCpy->setGridClass(openvdb::GRID_LEVEL_SET);


//     // Name the grid .
//     objectGridCpy->setName(objectGridCpy->getName() + std::string("COPY"));

//     openvdb::Index64 count = pBrb->m_params->materialPointCount; // TODO: this must be a parameter (possibly scaled by volume)

//     printf("materialPointCount = %d\n", (int)pBrb->m_params->materialPointCount);

//     //openvdb::FloatGrid::Ptr sampledBunnyMeshLevelSetGrid = bunnyMeshLevelSetGrid->deepCopy();
//     //sampledBunnyMeshLevelSetGrid->setName("sampledBunnyMeshLevelSetGrid");
//     printf("sdfInteriorMask\n");
//     openvdb::BoolGrid::Ptr interiorMaskBoolGridPtr = openvdb::tools::sdfInteriorMask(*objectGridCpy);

//     printf("topologyUnion\n");
//     objectGridCpy->topologyUnion(*interiorMaskBoolGridPtr);

//     printf("uniformPointScatter\n");
//     openvdb::points::PointDataGrid::Ptr points = openvdb::points::uniformPointScatter(*objectGridCpy, count, 0);
//     points->setName(objectGridCpy->getName() + "points");

//     // Verify the point count.
//     //openvdb::Index count = openvdb::points::pointCount(points->tree());
//     std::cout << "activeVoxelCount=" << objectGridCpy->tree().activeVoxelCount() << std::endl;
//     // std::cout << "leafCountPointTree=" << leafCountPointTree << std::endl;
//     //std::cout << "PointCount=" << count << std::endl;
// #if 0
//     // Create a VDB file object.
//     openvdb::io::File outfile(objectGridCpy->getName() + ".vdb");
//     // Add the grid pointer to a container.
//     openvdb::GridPtrVec grids;
//     grids.push_back(objectGridCpy);
//     grids.push_back(points);
//     // Write out the contents of the container.
//     outfile.write(grids);
// #endif
//     // save the points to a text file
//     // --------------------------

//     std::ofstream offFile("points.off");
//     offFile << "OFF\n";
//     offFile << count << " 0 0\n";

//     // Iterate over all the leaf nodes in the grid.
//     for (auto leafIter = points->tree().cbeginLeaf(); leafIter; ++leafIter) {
//         // Verify the leaf origin.
//         // std::cout << "Leaf" << leafIter->origin() << std::endl;
//         // Extract the position attribute from the leaf by name (P is position).
//         const openvdb::points::AttributeArray& positionArray = leafIter->constAttributeArray("P");
//         // Extract the radius attribute from the leaf by name (pscale is radius).
//         // const openvdb::points::AttributeArray& radiusArray =
//         //    leafIter->constAttributeArray("pscale");
//         // Create read-only handles for position and radius.
//         openvdb::points::AttributeHandle<openvdb::Vec3f> positionHandle(positionArray);
//         //openvdb::points::AttributeHandle<float> radiusHandle(radiusArray);
//         // Iterate over the point indices in the leaf.
//         for (auto indexIter = leafIter->beginIndexOn(); indexIter; ++indexIter) {
//             // Extract the voxel-space position of the point.
//             openvdb::Vec3f voxelPosition = positionHandle.get(*indexIter);
//             // Extract the world-space position of the voxel.
//             openvdb::Vec3d xyz = indexIter.getCoord().asVec3d();
//             // Compute the world-space position of the point.
//             openvdb::Vec3f worldPosition = points->transform().indexToWorld(voxelPosition + xyz);

//             offFile << worldPosition.x() << " " << worldPosition.y() << " " << worldPosition.z() << std::endl;
//             mpmParticlePositions.push_back(Eigen::Vector3d(worldPosition.x(), worldPosition.y(), worldPosition.z()));
//         }
//     }
// }

// // do simulation

// // generate crack surface

// // do mesh cutting and save mesh fragments


// // find vertices that are not connected by two faces
// std::set<int> findDupVertices( std::vector<Eigen::Vector3d> vertices , std::vector<Eigen::Vector3i> faces)
// {   
//     std::set<int> dupVertices;


//     std::map<std::string, int> edgeNum;
//     std::map<std::string, std::pair<int, int>> edgeVertices;
//     for(int f = 0; f < faces.size(); f++)
//     {
//         for(int fv = 0; fv < 3; fv++)
//         {
//             int vert1, vert2;
//             if(fv == 2)
//             {
//                 vert1 = faces[f][2];
//                 vert2 = faces[f][0];
//             }
//             else
//             {
//                 vert1 = faces[f][fv];
//                 vert2 = faces[f][fv + 1];
//             }
//             int vmin = std::min(vert1, vert2);
//             int vmax = std::max(vert1, vert2);
//             std::string edge = std::to_string(vmin) + "#" + std::to_string(vmax);
//             std::pair<int, int> edgeName = std::make_pair(vmin, vmax);

//             if(edgeNum.find(edge) == edgeNum.end())
//             {
//                 edgeNum[edge] = 1;
//                 edgeVertices[edge] = edgeName;
//             }
//             else
//             {
//                 edgeNum[edge] += 1;
//             }

//         }
        
//     }


//     std::map<std::string, int> ::iterator it;
//     for (it = edgeNum.begin(); it != edgeNum.end(); it++)
//     {
//         if(it->second == 1)
//         {
//             std::pair<int, int> twoVertices = edgeVertices[it->first];
//             dupVertices.insert(twoVertices.first);
//             dupVertices.insert(twoVertices.second);
//         }
//     }


//     return dupVertices;

// }




// int runFractureSim(Object* pBrb, int timeStepIndex, GenericMesh* finalCrackSurfaceMesh /*output*/, const Params& simParams)
// {
//     int numFragments = 0;
//     SCOPED_TIMER(__FUNCTION__);



//     // 1) init material particles from object level set
//     std::cout << "init material particles from object level set" << std::endl;
//     std::vector<Eigen::Vector3d> mpmParticlePositions;
//     generateMaterialParticles(mpmParticlePositions, pBrb);



//     // 3) get material parameters
//     float youngsModulus = pBrb->m_params->youngsModulus;
//     float poissonsRatio = pBrb->m_params->poissonsRatio;
//     float density = pBrb->m_params->density;
//     float volume = pBrb->levelSetVolume;
//     float mass = pBrb->mass; // density * volume
//     Eigen::Vector3d centerOfMass = pBrb->centerOfMass.cast<double>();

//     // 4) init MPM simulator

//     std::cout << "init MPM simulator" << std::endl;
//     struct parametersSim param;
// #if defined(MULTI_MPM_THREADED)
//     param.numOfMpmThreads = std::thread::hardware_concurrency();
// #else
//     param.numOfMpmThreads = 1;
// #endif 

//     std::cout<< "param.numOfMpmThreads = " <<param.numOfMpmThreads<<std::endl;
//     param.dt = pBrb->m_params->mpmTimestep;
//     param.dpx = std::pow(volume / mpmParticlePositions.size(), 1.0 / 3.0);
//     param.dx = param.dpx * 2;

//     param.DP = param.dx * param.dx / 4.0;
//     param.dcx = param.dpx * 1.5;
//     param.vol = volume / mpmParticlePositions.size();
//     param.particle_mass = mass / mpmParticlePositions.size();

//     param.E = youngsModulus;
//     param.nu = poissonsRatio;
//     param.mu = param.E / (2 * (1 + param.nu));
//     param.lambda = param.E * param.nu / ((1 + param.nu) * (1 - 2 * param.nu));
//     param.K = 2 / 3 * param.mu + param.lambda;
    

//     // local damage field parameters
//     param.thetaf = pBrb->m_params->referenceStress;
//     param.Gf = pBrb->m_params->modeIFractureEnergy;
//     param.lch = std::sqrt(3) * param.dx;
//     param.HsBar = param.thetaf * param.thetaf / 2 / param.E / param.Gf;
//     param.Hs = param.HsBar * param.lch / (1 - param.HsBar * param.lch);
//     param.damageThreshold = pBrb->m_params->particleDamageUpdateThreshold; // after this threshold, the damage value is updated by a sigmoid function to let it approach 1 but never reach 1
//     param.sigmoidK = 5; // this parameter control the curevature of the sigmoid function. It is recommend that it is bigger than 5

//     // chevron marks parameters
//     param.drx = param.dpx * pBrb->m_params->meshRefineSpacing; // refined mesh resolution
//     param.pathLength = pBrb->m_params->chevronMarkMaxLength; // the maximum path length
//     param.stepSize = param.drx; // stepsize
//     param.divertAngle = 0.3; // cos(angle) where angle is the angle between this step and next step
//     param.influenceRadius = pBrb->m_params->influenceRadius * param.drx; // the radius during which no vertices should be visited
//     param.smoothRadius = pBrb->m_params->smoothRadius * param.drx; // background grid spacing for smooth height value
//     param.pertubationMagitude = pBrb->m_params->pertubationMagitude * param.drx; // pertubation magnitude for each vertex
//     param.ratioOfDeepenMarks = pBrb->m_params->ratioOfDeepenMarks; // pertubation magnitude for each vertex
//     param.deepenMagnitude = pBrb->m_params->deepenMagnitude; // pertubation magnitude for each vertex



//     std::vector<Particle> particleVec; // store particles' positions in a vector and get the bounding box corner
//     std::vector<Eigen::Vector3d> mpmParticlePositions_Proj; // projected particle position
//     // project particles' positions into the area starting from 0
//     double xmin = mpmParticlePositions[0][0], xmax = mpmParticlePositions[0][0], ymin = mpmParticlePositions[0][1], ymax = mpmParticlePositions[0][1], zmin = mpmParticlePositions[0][2], zmax = mpmParticlePositions[0][2];
//     for (int km = 0; km < mpmParticlePositions.size(); km++) 
//     {
//         Eigen::Vector3d ipos = mpmParticlePositions[km];
//         xmin = std::min(xmin, ipos[0]);
//         xmax = std::max(xmax, ipos[0]);
//         ymin = std::min(ymin, ipos[1]);
//         ymax = std::max(ymax, ipos[1]);
//         zmin = std::min(zmin, ipos[2]);
//         zmax = std::max(zmax, ipos[2]);
//     }
//     Eigen::Vector3d minCoordinate = { xmin - 50 * param.dx, ymin - 50 * param.dx, zmin - 50 * param.dx };
//     for (int km = 0; km < mpmParticlePositions.size(); km++) {
//         Eigen::Vector3d ipos = mpmParticlePositions[km] - minCoordinate;
//         Eigen::Vector3d ivel = { 0, 0, 0 };
//         particleVec.push_back(Particle(ipos, ivel, param.particle_mass, 0, 0));
//         mpmParticlePositions_Proj.push_back(ipos);
//     }
//     Eigen::Vector3d length = {xmax - xmin + 100* param.dx,  ymax - ymin + 100* param.dx, zmax - zmin + 100 * param.dx};
//     param.length = length;



//     // 2) init boundary conditions
//     // force value & location (in object's local coordinates)
//     std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> boundaryConditions;
//     for (std::map<int, Eigen::Vector3f>::const_iterator i = pBrb->contactTractions.cbegin();
//          i != pBrb->contactTractions.cend();
//          ++i) {

//         if (i->second.norm() > 0) {
//             int triangleID = i->first;
//             Eigen::Vector3d force = i->second.cast<double>();
//             Eigen::Vector3d location = pBrb->elemCtr.at(triangleID).cast<double>();
//             boundaryConditions.emplace_back(force, location - minCoordinate);
//         }
//     }
//     param.appliedForce = boundaryConditions;


//     if(boundaryConditions.empty())
//     {
//         printf("no boundary conditions, return\n");
//         return -1;
//     }



//     // // 0) write the object in .obj file
//     // std::string objectName = getName(pBrb->rigidBodyPtr).c_str();
//     // trimesh::TriMesh* collisionMeshOut =  pBrb->collInUse.trimeshLowRes;
//     // collisionMeshOut->write("original.obj");



//     // 5) run MPM simulator

//     std::cout << "run MPM simulator" << std::endl;
//     int timeStep = std::min((int)round(pBrb->contactDuration / param.dt), pBrb->m_params->mpmTimeStepsMax); // use the minimum value of calculated timestep and predefined timestep 
//     advanceStep(param, &particleVec, timeStep, pBrb->m_params->fractureSimDamageRatioThreshold, pBrb->m_params->forceScaleFactor);



//     // 6) run crack extraction algorithm
//     std::cout << "run crack extraction algorithm and refine mesh.." << std::endl;
//     std::string fielName = "particles" + std::to_string(timeStep) + ".txt";
//     std::ofstream outfile0(fielName, std::ios::trunc);
//     for (int f = 0; f < particleVec.size(); f++) 
//     {
//         outfile0 << std::scientific << std::setprecision(8) << mpmParticlePositions_Proj[f][0] << "  " << mpmParticlePositions_Proj[f][1] << " " << mpmParticlePositions_Proj[f][2] << " " << particleVec[f].Dp << " " << particleVec[f].Dg << "\n";
//     };
//     outfile0.close();




    
//     struct meshFaceRefined refinedMeshResult;
//     std::cout << "run crack extraction algorithm and refine mesh.." << std::endl;
//     try 
//     {
// 		refinedMeshResult = extractCrackSurface_PartialCracks(fielName, &mpmParticlePositions_Proj, param);
// 	}
// 	catch (...) 
//     {
// 		std::cout << "Failed to run triangle.h library" << std::endl;
// 	}

//     std::vector<Eigen::Vector3d> vertices;
//     std::vector<std::vector<int>> faces;

//     try 
//     {
// 		if (refinedMeshResult.verticesRefined.size() != 0) {
//         std::cout << "calculate order information.." << std::endl;
//         std::vector<Particle> particleVec_Copy = particleVec;
//         for (int i = 0; i < particleVec_Copy.size(); i++) {
//             particleVec_Copy[i].Dp = particleVec[i].Dg;
//         }
//         struct std::vector<verticesInformation> verticesVec = calculateOrderInfor(param, &particleVec_Copy, &refinedMeshResult);

//         std::cout << "finding chevron marks.." << std::endl;
//         crackSurface res = findChevronMarks_partialCracks(param, &verticesVec, &refinedMeshResult);

//         vertices = res.vertices; // vertices of the crack surface
//         for (int i = 0; i < vertices.size(); i++) // project vertices back
//         {
//             vertices[i] += minCoordinate;
//         }
//         faces = res.faces; // faces of the crack surface
//     }
// 	}
// 	catch (...) 
//     {
// 		std::cout << "Failed to compute chevronmarks" << std::endl;
// 	}



//     //writeOffFile(vertices, faces, "cracksurface");


//     if (vertices.size() == 0 || faces.size() == 0) // the crack algorithm may not generate any crack surface[this is one limitation of the algorithm]
//     {
//         return -1;
//     }

//     // 7) make sure crack mesh is clean and return result
//     TIMESTACK_PUSH("cleanCrackMesh");
//     std::cout << " make sure crack mesh is clean and return result" << std::endl;
//     finalCrackSurfaceMesh->m.vertices = vertices; // copy
//     finalCrackSurfaceMesh->m.faces = faces;
//     cleanGenericMesh(*finalCrackSurfaceMesh);
//     triangulateGenericMesh(*finalCrackSurfaceMesh); // triangulate faces for OpenVDB conversion

//     trimesh::TriMesh tm;
//     toTrimesh(tm, *finalCrackSurfaceMesh);

//     //tm.write("TriCrackMeshPreview.obj");

// #if 0
//     float t_max = maxTime, t_step; // time-steps
//     int addedCracks = 0, chk;
//     if (contactDuration < t_max) {
//         t_max = contactDuration
//     };

//     // controls how many timesteps to execute before we stop teh fracture sim
//     float maxPropagationLen = 1; // TODO
//     t_step = pBrb->getFractureTimeStepSize(maxPropagationSteps);
//     int steps = (int)((t_max / t_step) + 0.5); // round

//     printf("\n%% fracture simulation for %s_%d ... %d time steps", ((string*)rb->getUserPointer())->c_str(), rbTimeCode, steps);

//     if (steps == 0) {
//         printf("\n%% contact duration or max time too short (%.3lgs)", t_max);
//         return -1;
//     }

    

//     if (steps > 0) {
//         std::map<int, Eigen::Vector3d> q;130ncedContactTractions(q);

//         // TODO: Set simulation boundary conditions (force etc.)

//         bool done = false;
//         for (int k = 0; k < steps && !done; ++k) { // do fracture steps
 
//             printf("\n%% fracturing (%d/%d) ... ", k + 1, steps);

//             // TODO: do simulation, calculate stress, update damage etc,
            
//             printf("\n%% compute MPM timestep    ...");
//         }

//         // TODO: do mesh cutting and save mesh fragments
//     }
// #endif
//     TIMESTACK_POP();
//     return numFragments;
// }



// int runFractureSim_MeshCutting(Object* pBrb, int timeStepIndex, std::vector<meshObjFormat>* fragmentsCollision, std::vector<meshObjFormat>* fragmentsCollision_NoTri, std::vector<meshObjFormat>* fragmentsRendering, std::vector<std::vector<Eigen::Vector3d>>* perturbVerticesMag, const Params& simParams)
// {

//     int numFragments = 0;
//     SCOPED_TIMER(__FUNCTION__);


//     // 1) init material particles from object level set
//     std::cout << "init material particles from object level set" << std::endl;
//     std::vector<Eigen::Vector3d> mpmParticlePositions;
//     generateMaterialParticles(mpmParticlePositions, pBrb);



//     // 3) get material parameters
//     float youngsModulus = pBrb->m_params->youngsModulus;
//     float poissonsRatio = pBrb->m_params->poissonsRatio;
//     float density = pBrb->m_params->density;
//     float volume = pBrb->levelSetVolume;
//     float mass = pBrb->mass; // density * volume
//     Eigen::Vector3d centerOfMass = pBrb->centerOfMass.cast<double>();

//     // 4) init MPM simulator

//     std::cout << "init MPM simulator" << std::endl;

//     struct parametersSim param;
// #if defined(MULTI_MPM_THREADED)
//     param.numOfMpmThreads = (int)std::thread::hardware_concurrency();
// #else
//     param.numOfMpmThreads = 1;
// #endif
//     param.dt = pBrb->m_params->mpmTimestep;
//     param.dpx = std::pow(volume / mpmParticlePositions.size(), 1.0 / 3.0);
//     param.dx = param.dpx * 2;

//     param.DP = param.dx * param.dx / 4.0;
//     param.dcx = param.dpx * 1.5;
//     param.vol = volume / mpmParticlePositions.size();
//     param.particle_mass = mass / mpmParticlePositions.size();

//     param.E = youngsModulus;
//     param.nu = poissonsRatio;
//     param.mu = param.E / (2 * (1 + param.nu));
//     param.lambda = param.E * param.nu / ((1 + param.nu) * (1 - 2 * param.nu));
//     param.K = 2 / 3 * param.mu + param.lambda;
    

//     // local damage field parameters
//     param.thetaf = pBrb->m_params->referenceStress;
//     param.Gf = pBrb->m_params->modeIFractureEnergy;
//     param.lch = std::sqrt(3) * param.dx;
//     param.HsBar = param.thetaf * param.thetaf / 2 / param.E / param.Gf;
//     param.Hs = param.HsBar * param.lch / (1 - param.HsBar * param.lch);
//     param.damageThreshold = pBrb->m_params->particleDamageUpdateThreshold; // after this threshold, the damage value is updated by a sigmoid function to let it approach 1 but never reach 1
//     param.sigmoidK = 5; // this parameter control the curevature of the sigmoid function. It is recommend that it is bigger than 5

//     // chevron marks parameters
//     param.drx = param.dpx * pBrb->m_params->meshRefineSpacing; // refined mesh resolution
//     param.pathLength = pBrb->m_params->chevronMarkMaxLength; // the maximum path length
//     param.stepSize = param.drx; // stepsize
//     param.divertAngle = 0.3; // cos(angle) where angle is the angle between this step and next step
//     param.influenceRadius = pBrb->m_params->influenceRadius * param.drx; // the radius during which no vertices should be visited
//     param.smoothRadius = pBrb->m_params->smoothRadius * param.drx; // background grid spacing for smooth height value
//     param.pertubationMagitude = pBrb->m_params->pertubationMagitude * param.drx; // pertubation magnitude for each vertex
//     param.ratioOfDeepenMarks = pBrb->m_params->ratioOfDeepenMarks; // pertubation magnitude for each vertex
//     param.deepenMagnitude = pBrb->m_params->deepenMagnitude; // pertubation magnitude for each vertex


//     std::vector<Particle> particleVec; // store particles' positions in a vector and get the bounding box corner
//     std::vector<Eigen::Vector3d> mpmParticlePositions_Proj; // projected particle position
//     // project particles' positions into the area starting from 0
//     double xmin = mpmParticlePositions[0][0], xmax = mpmParticlePositions[0][0], ymin = mpmParticlePositions[0][1], ymax = mpmParticlePositions[0][1], zmin = mpmParticlePositions[0][2], zmax = mpmParticlePositions[0][2];
//     for (int km = 0; km < mpmParticlePositions.size(); km++) 
//     {
//         Eigen::Vector3d ipos = mpmParticlePositions[km];
//         xmin = std::min(xmin, ipos[0]);
//         xmax = std::max(xmax, ipos[0]);
//         ymin = std::min(ymin, ipos[1]);
//         ymax = std::max(ymax, ipos[1]);
//         zmin = std::min(zmin, ipos[2]);
//         zmax = std::max(zmax, ipos[2]);
//     }
//     Eigen::Vector3d minCoordinate = { xmin - 50 * param.dx, ymin - 50 * param.dx, zmin - 50 * param.dx };
//     for (int km = 0; km < mpmParticlePositions.size(); km++) {
//         Eigen::Vector3d ipos = mpmParticlePositions[km] - minCoordinate;
//         Eigen::Vector3d ivel = { 0, 0, 0 };
//         particleVec.push_back(Particle(ipos, ivel, param.particle_mass, 0, 0));
//         mpmParticlePositions_Proj.push_back(ipos);
//     }
//     Eigen::Vector3d length = {xmax - xmin + 100* param.dx,  ymax - ymin + 100* param.dx, zmax - zmin + 100 * param.dx};
//     param.length = length;



//     // 2) init boundary conditions
//     // force value & location (in object's local coordinates)
//     std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> boundaryConditions;
//     for (std::map<int, Eigen::Vector3f>::const_iterator i = pBrb->contactTractions.cbegin();
//          i != pBrb->contactTractions.cend();
//          ++i) {

//         if (i->second.norm() > 0) {
//             int triangleID = i->first;
//             Eigen::Vector3d force = i->second.cast<double>();
//             Eigen::Vector3d location = pBrb->elemCtr.at(triangleID).cast<double>();
//             boundaryConditions.emplace_back(force, location - minCoordinate);
//         }
//     }
//     param.appliedForce = boundaryConditions;

//     if(boundaryConditions.empty())
//     {
//         printf("no boundary conditions, return\n");
//         return -1;
//     }




// #if defined(USE_MCUT)
//     // 0) write the object in .obj file
//     std::string objectName = getName(pBrb->rigidBodyPtr).c_str();
//     meshObjFormat collisionMeshParent;
//     for (int i = 0; i < pBrb->cuttingMesh_NoTri->vertices.size(); i++)
// 	{
//         Eigen::Vector3d vert = { pBrb->cuttingMesh_NoTri->vertices[i][0] - minCoordinate[0], pBrb->cuttingMesh_NoTri->vertices[i][1] - minCoordinate[1], pBrb->cuttingMesh_NoTri->vertices[i][2] - minCoordinate[2] };
//         collisionMeshParent.vertices.push_back(vert);
//     } 
//     collisionMeshParent.faces = pBrb->cuttingMesh_NoTri->faces; 
//     std::vector<Eigen::Vector3d> perturbVerticesMagParent = pBrb->perturbVerticesMag;
// #else
//     // 0) write the object in .obj file
//     std::string objectName = getName(pBrb->rigidBodyPtr).c_str();
//     meshObjFormat collisionMeshParent;
//     for (int i = 0; i < pBrb->cuttingMesh->vertices.size(); ++i)
//     {
//         Eigen::Vector3d vert = {pBrb->cuttingMesh->vertices[i][0] - minCoordinate[0], pBrb->cuttingMesh->vertices[i][1] - minCoordinate[1], pBrb->cuttingMesh->vertices[i][2] - minCoordinate[2]};
//         collisionMeshParent.vertices.push_back(vert);
//     }
//     for (int i = 0; i < pBrb->cuttingMesh->faces.size(); ++i)
//     {
//         std::vector<int> face;
//         face.push_back(pBrb->cuttingMesh->faces[i][0]);
//         face.push_back(pBrb->cuttingMesh->faces[i][1]);
//         face.push_back(pBrb->cuttingMesh->faces[i][2]);
//         collisionMeshParent.faces.push_back(face);
//     }
//     std::vector<Eigen::Vector3d> perturbVerticesMagParent = pBrb->perturbVerticesMag;
// #endif




// std::cout << "mpm  = "<< 4 << std::endl;

//     // 5) run MPM simulator

//     std::cout << "run MPM simulator" << std::endl;
//     int timeStep = std::min((int)round(pBrb->contactDuration / param.dt), pBrb->m_params->mpmTimeStepsMax); // use the minimum value of calculated timestep and predefined timestep 
//     advanceStep(param, &particleVec, timeStep, pBrb->m_params->fractureSimDamageRatioThreshold, pBrb->m_params->forceScaleFactor);



//     // 6) run crack extraction algorithm
//     std::cout << "run crack extraction algorithm and refine mesh.." << std::endl;
//     std::string fielName = "particles" + std::to_string(timeStep) + ".txt";
//     std::ofstream outfile0(fielName, std::ios::trunc);
//     for (int f = 0; f < particleVec.size(); f++) 
//     {
//         ASSERT(f < mpmParticlePositions_Proj.size());
//         outfile0 << std::scientific << std::setprecision(8) << mpmParticlePositions_Proj[f][0] << "  " << mpmParticlePositions_Proj[f][1] << " " << mpmParticlePositions_Proj[f][2] << " " << particleVec[f].Dp << " " << particleVec[f].Dg << "\n";
//     };
//     outfile0.close();
//     std::tuple<bool, meshFaceRefined, std::vector<meshObjFormat> > result = extractCrackSurface(fielName, &mpmParticlePositions_Proj, param);
//     struct meshFaceRefined refinedMeshResult = std::get<1>(result);
//     std::vector<meshObjFormat> fragmentsRaw = std::get<2>(result);

//     std::cout<<"0"<<std::endl;
//     // 7) get final fragments




// #if defined(USE_MCUT)
    

//     if(fragmentsRaw.size() == 1)
//     {
//         // collision mesh / cutting mesh
//         meshObjFormat collisionMeshParent;
//         for (int i = 0; i < pBrb->renderingMesh->vertices.size(); i++)
//         {
//             Eigen::Vector3d vert = { pBrb->renderingMesh->vertices[i][0], pBrb->renderingMesh->vertices[i][1], pBrb->renderingMesh->vertices[i][2]};
//             collisionMeshParent.vertices.push_back(vert);
//         }
//         for (int i = 0; i < pBrb->renderingMesh->faces.size(); i++)
//         {
//             std::vector<int> face;
//             face.push_back(pBrb->renderingMesh->faces[i][0]);
//             face.push_back(pBrb->renderingMesh->faces[i][1]);
//             face.push_back(pBrb->renderingMesh->faces[i][2]);
//             collisionMeshParent.faces.push_back(face);
//         }
//         (*fragmentsCollision).push_back(collisionMeshParent);


//         // cutting mesh without triangulation
//         meshObjFormat collisionMeshOut_NoTri;
//         collisionMeshOut_NoTri.vertices = pBrb->cuttingMesh_NoTri->vertices;
//         collisionMeshOut_NoTri.faces = pBrb->cuttingMesh_NoTri->faces;
//         (*fragmentsCollision_NoTri).push_back(collisionMeshOut_NoTri);


//         // rendering mesh
//         trimesh::TriMesh* renderingMeshOut =  pBrb->renderingMesh;
//         meshObjFormat renderingMeshParent;
//         for (int i = 0; i < renderingMeshOut->vertices.size(); i++)
//         {
//             Eigen::Vector3d vert = { renderingMeshOut->vertices[i][0], renderingMeshOut->vertices[i][1], renderingMeshOut->vertices[i][2]};
//             renderingMeshParent.vertices.push_back(vert);
//         }
//         for (int i = 0; i < renderingMeshOut->faces.size(); i++)
//         {
//             std::vector<int> face;
//             face.push_back(renderingMeshOut->faces[i][0]);
//             face.push_back(renderingMeshOut->faces[i][1]);
//             face.push_back(renderingMeshOut->faces[i][2]);
//             renderingMeshParent.faces.push_back(face);
//         }
//         (*fragmentsRendering).push_back(renderingMeshParent);


        
//     }
//     else
//     {
//         std::vector<std::pair<std::set<int> , std::map<int, int> >> crackSurVertices;
//         cutObject_MCUT(param ,objectName, &collisionMeshParent, &fragmentsRaw,  fragmentsCollision, &crackSurVertices, fragmentsCollision_NoTri);


//         std::cout<<"1"<<std::endl;


//         if(std::get<0>(result) == false) // fail to use triangle.h
//         {
//             std::cout << "Failed to use triangle.h. Use cutted objects directly" << std::endl;
//             for(int mm =0;mm < fragmentsCollision->size(); mm++)
//             {
//                 ASSERT(mm < (*fragmentsCollision).size());

//                 for(int nn =0;nn < (*fragmentsCollision)[mm].vertices.size(); nn++)
//                 {
//                     ASSERT(nn < (*fragmentsCollision)[mm].vertices.size());
//                     (*fragmentsCollision)[mm].vertices[nn] += minCoordinate;
//                 }

//                 for(int nn =0;nn < (*fragmentsCollision_NoTri)[mm].vertices.size(); nn++)
//                 {
//                      ASSERT(nn < (*fragmentsCollision_NoTri)[mm].vertices.size());
//                     (*fragmentsCollision_NoTri)[mm].vertices[nn] += minCoordinate;
//                 }

//                 (*fragmentsRendering).push_back((*fragmentsCollision)[mm]);
//             }

//         }
//         else
//         {

//             std::cout<<"2"<<std::endl;
//             std::vector<Particle> particleVec_Copy = particleVec;
//             #pragma omp parallel for num_threads(param.numOfMpmThreads)
//             for (int i = 0; i < particleVec_Copy.size(); i++) 
//             {
//                 ASSERT(i < particleVec_Copy.size());
//                 ASSERT(i < mpmParticlePositions_Proj.size());
//                 particleVec_Copy[i].pos = mpmParticlePositions_Proj[i];
//                 particleVec_Copy[i].Dp = particleVec[i].Dg;
//             }

//             std::cout<<"3"<<std::endl;
//             struct std::vector<verticesInformation> gridInformation = calculateOrderInfor(param, &particleVec_Copy, &refinedMeshResult);
//             std::map<int, int> gridMapChevron;
//             std::vector<Grid> nodesVecChevron;
//             std::cout<<"4"<<std::endl;
//             try 
//             {
//                 std::cout<<"5"<<std::endl;
//                 chevronMarksGrid(param, &gridInformation, &refinedMeshResult, &gridMapChevron, &nodesVecChevron);
//                 // add chevron marks to each fragments
//                 std::cout<<"6"<<std::endl;

//                 addChevronFragments(param, fragmentsCollision, &crackSurVertices,fragmentsRendering, &perturbVerticesMagParent, perturbVerticesMag, &gridMapChevron, &nodesVecChevron);
//             }
//             catch (...) 
//             {
//                 std::cout << "Failed to compute chevronmarks. Use cutted objects directly" << std::endl;
//                 fragmentsRendering = fragmentsCollision;
//             }

//             std::cout<<"7"<<std::endl;
//             #pragma omp parallel for num_threads(param.numOfMpmThreads)
//             for(int mm =0;mm < (*fragmentsCollision).size(); mm++)
//             {
//                 ASSERT(mm < (*fragmentsCollision).size());
//                 for(int nn =0;nn < (*fragmentsCollision)[mm].vertices.size(); nn++)
//                 {
//                     ASSERT(nn < (*fragmentsCollision)[mm].vertices.size());
//                     (*fragmentsCollision)[mm].vertices[nn] += minCoordinate;
//                 }
//             }



//             std::cout<<"8"<<std::endl;
//             #pragma omp parallel for num_threads(param.numOfMpmThreads)
//             for(int mm =0;mm < (*fragmentsRendering).size(); mm++)
//             {
//                 for(int nn =0;nn <  (*fragmentsRendering)[mm].vertices.size(); nn++)
//                 {
//                     (*fragmentsRendering)[mm].vertices[nn] += minCoordinate;
//                 }
//             }

//             #pragma omp parallel for num_threads(param.numOfMpmThreads)
//             for(int mm =0;mm < (*fragmentsCollision_NoTri).size(); mm++)
//             {
//                 for(int nn =0;nn <  (*fragmentsCollision_NoTri)[mm].vertices.size(); nn++)
//                 {
//                     (*fragmentsCollision_NoTri)[mm].vertices[nn] += minCoordinate;
//                 }
//             }

//         }


//     }

    
// #else
    

// #endif


    


//     for(int k = 0; k <  (*fragmentsCollision).size(); k++)
//     {
//         std::string name = "fragment_Collision" + std::to_string(k);
// 		writeOffFile((*fragmentsCollision)[k].vertices, (*fragmentsCollision)[k].faces, name);
//     }

//     for(int k = 0; k <  (*fragmentsCollision).size(); k++)
//     {
//         std::string name = "fragment_Collision_NoTri" + std::to_string(k);
// 		writeOffFile((*fragmentsCollision_NoTri)[k].vertices, (*fragmentsCollision_NoTri)[k].faces, name);
//     }

//     for(int k = 0; k <  (*fragmentsRendering).size(); k++)
//     {
//         std::string name = "fragment_chevron" + std::to_string(k);
// 		writeOffFile( (*fragmentsRendering)[k].vertices,  (*fragmentsRendering)[k].faces, name);
//     }


//     numFragments = (*fragmentsCollision).size();

//     if ((*fragmentsCollision).size() == 0) // there may be no fragments
//     {
//         return -1;
//     }
//     return numFragments;
// }





// void getSurfaceMeshFromVdb(openvdb::FloatGrid::Ptr bareMeshVdbGridPtr,int decimateTarget, meshObjFormat& fragmentVolume)
// {
//     openvdb::tools::VolumeToMesh volumeToMeshHandle;
//     volumeToMeshHandle(*bareMeshVdbGridPtr);

//     trimesh::TriMesh* pMesh = new trimesh::TriMesh;

//     openvdb::tools::PointList *verts = &volumeToMeshHandle.pointList();
//     openvdb::tools::PolygonPoolList *polys = &volumeToMeshHandle.polygonPoolList();


//     for (size_t i = 0; i < volumeToMeshHandle.pointListSize(); i++)
//     {
//         openvdb::Vec3s &v = (*verts)[i];
//         // NOTE: This is a hack! The y coord is negated because of mixup with blender's XZY coordinate system ****************************************************************************************
//         pMesh->vertices.push_back(trimesh::vec3(v[0], v[1], v[2]));

//     }

//     for (size_t i = 0; i < volumeToMeshHandle.polygonPoolListSize(); i++)
//     {

//         for (size_t ndx = 0; ndx < (*polys)[i].numTriangles(); ndx++)
//         {
//             openvdb::Vec3I *p = &((*polys)[i].triangle(ndx));
//         }

//         for (size_t ndx = 0; ndx < (*polys)[i].numQuads(); ndx++)
//         {
//             openvdb::Vec4I *p = &((*polys)[i].quad(ndx));

//             trimesh::TriMesh::Face f0;
//             f0[0] = p->z();
//             f0[1] = p->y();
//             f0[2] = p->x();
//             trimesh::TriMesh::Face f1;
//             f1[0] = p->w();
//             f1[1] = p->z();
//             f1[2] = p->x();

//             pMesh->faces.push_back(f0);
//             pMesh->faces.push_back(f1);

//         }
//     }


//     // decimate mesh using libigl
//     Eigen::MatrixXd V(pMesh->vertices.size(), 3);
//     Eigen::MatrixXi F(pMesh->faces.size(), 3);
//     for(int i = 0; i < pMesh->vertices.size(); i++)
//     {
//         V(i,0) = pMesh->vertices[i][0];
//         V(i,1) = pMesh->vertices[i][1];
//         V(i,2) = pMesh->vertices[i][2];
//     }
//     for(int i = 0; i < pMesh->faces.size(); i++)
//     {
//         F(i,0) = pMesh->faces[i][0];
//         F(i,1) = pMesh->faces[i][1];
//         F(i,2) = pMesh->faces[i][2];
//     }

//     std::cout<<"pMesh->vertices.size() = " << pMesh->vertices.size()<<std::endl;
//     std::cout<<"pMesh->faces.size() = " << pMesh->faces.size()<<std::endl;


//     Eigen::MatrixXd U;
//     Eigen::MatrixXi G;
//     Eigen::VectorXi J;
//     Eigen::VectorXi I;
//     igl::decimate(V, F, decimateTarget, U, G, J, I);
//     for (int i = 0; i < U.rows(); ++i)
//     {
//         Eigen::Vector3d p = {U(i, 0), U(i, 1), U(i, 2)};
//         fragmentVolume.vertices.push_back(p);
//     }
//     for (int i = 0; i < G.rows(); ++i)
//     {
//         std::vector<int> f;
//         f.push_back(G(i, 0));
//         f.push_back(G(i, 1));
//         f.push_back(G(i, 2));
//         fragmentVolume.faces.push_back(f);
//     }

//     std::cout<<"U.rows() = " << U.rows()<<std::endl;
//     std::cout<<"G.rows() = " << G.rows()<<std::endl;


// }



// int runFractureSim_OpenVdbFullyCut(Object* pBrb, int timeStepIndex, std::vector<GenericMesh>* finalCrackSurfaceMesh /*output*/,  const Params& simParams)
// {
//     int numFragments = 0;
//     SCOPED_TIMER(__FUNCTION__);


//     std::string objectName1 = getName(pBrb->rigidBodyPtr).c_str();
//     if(objectName1 == "lionBall_") // becasue the bullet is set as a breakable object while we don't want to break it
//     {
//         return -1;
//     }

//     // 1) init material particles from object level set
//     std::cout << "init material particles from object level set" << std::endl;
//     std::vector<Eigen::Vector3d> mpmParticlePositions;
//     generateMaterialParticles(mpmParticlePositions, pBrb);



//     // 3) get material parameters
//     float youngsModulus = pBrb->m_params->youngsModulus;
//     float poissonsRatio = pBrb->m_params->poissonsRatio;
//     float density = pBrb->m_params->density;
//     float volume = pBrb->levelSetVolume;
//     float mass = pBrb->mass; // density * volume
//     Eigen::Vector3d centerOfMass = pBrb->centerOfMass.cast<double>();

//     // 4) init MPM simulator

//     std::cout << "init MPM simulator" << std::endl;
//     struct parametersSim param;
// #if defined(MULTI_MPM_THREADED)
//     param.numOfMpmThreads = std::thread::hardware_concurrency();
// #else
//     param.numOfMpmThreads = 1;
// #endif 

//     std::cout<< "param.numOfMpmThreads = " <<param.numOfMpmThreads<<std::endl;
//     param.dt = pBrb->m_params->mpmTimestep;
//     param.dpx = std::pow(volume / mpmParticlePositions.size(), 1.0 / 3.0);
//     param.dx = param.dpx * 2;

//     param.DP = param.dx * param.dx / 4.0;
//     param.dcx = param.dpx * 1.5;
//     param.vol = volume / mpmParticlePositions.size();
//     param.particle_mass = mass / mpmParticlePositions.size();

//     param.E = youngsModulus;
//     param.nu = poissonsRatio;
//     param.mu = param.E / (2 * (1 + param.nu));
//     param.lambda = param.E * param.nu / ((1 + param.nu) * (1 - 2 * param.nu));
//     param.K = 2 / 3 * param.mu + param.lambda;
    

//     // local damage field parameters
//     param.thetaf = pBrb->m_params->referenceStress;
//     param.Gf = pBrb->m_params->modeIFractureEnergy;
//     param.lch = std::sqrt(3) * param.dx;
//     param.HsBar = param.thetaf * param.thetaf / 2 / param.E / param.Gf;
//     param.Hs = param.HsBar * param.lch / (1 - param.HsBar * param.lch);
//     param.damageThreshold = pBrb->m_params->particleDamageUpdateThreshold; // after this threshold, the damage value is updated by a sigmoid function to let it approach 1 but never reach 1
//     param.sigmoidK = 5; // this parameter control the curevature of the sigmoid function. It is recommend that it is bigger than 5

//     // chevron marks parameters
//     param.drx = param.dpx * pBrb->m_params->meshRefineSpacing; // refined mesh resolution
//     param.pathLength = pBrb->m_params->chevronMarkMaxLength; // the maximum path length
//     param.stepSize = param.drx; // stepsize
//     param.divertAngle = 0.3; // cos(angle) where angle is the angle between this step and next step
//     param.influenceRadius = pBrb->m_params->influenceRadius * param.drx; // the radius during which no vertices should be visited
//     param.smoothRadius = pBrb->m_params->smoothRadius * param.drx; // background grid spacing for smooth height value
//     param.pertubationMagitude = pBrb->m_params->pertubationMagitude * param.drx; // pertubation magnitude for each vertex
//     param.ratioOfDeepenMarks = pBrb->m_params->ratioOfDeepenMarks; // pertubation magnitude for each vertex
//     param.deepenMagnitude = pBrb->m_params->deepenMagnitude; // pertubation magnitude for each vertex



//     std::vector<Particle> particleVec; // store particles' positions in a vector and get the bounding box corner
//     std::vector<Eigen::Vector3d> mpmParticlePositions_Proj; // projected particle position
//     // project particles' positions into the area starting from 0
//     double xmin = mpmParticlePositions[0][0], xmax = mpmParticlePositions[0][0], ymin = mpmParticlePositions[0][1], ymax = mpmParticlePositions[0][1], zmin = mpmParticlePositions[0][2], zmax = mpmParticlePositions[0][2];
//     for (int km = 0; km < mpmParticlePositions.size(); km++) 
//     {
//         Eigen::Vector3d ipos = mpmParticlePositions[km];
//         xmin = std::min(xmin, ipos[0]);
//         xmax = std::max(xmax, ipos[0]);
//         ymin = std::min(ymin, ipos[1]);
//         ymax = std::max(ymax, ipos[1]);
//         zmin = std::min(zmin, ipos[2]);
//         zmax = std::max(zmax, ipos[2]);
//     }
//     Eigen::Vector3d minCoordinate = { xmin - 50 * param.dx, ymin - 50 * param.dx, zmin - 50 * param.dx };
//     for (int km = 0; km < mpmParticlePositions.size(); km++) {
//         Eigen::Vector3d ipos = mpmParticlePositions[km] - minCoordinate;
//         Eigen::Vector3d ivel = { 0, 0, 0 };
//         particleVec.push_back(Particle(ipos, ivel, param.particle_mass, 0, 0));
//         mpmParticlePositions_Proj.push_back(ipos);
//     }
//     Eigen::Vector3d length = {xmax - xmin + 100* param.dx,  ymax - ymin + 100* param.dx, zmax - zmin + 100 * param.dx};
//     param.length = length;
//     param.minCoordinate = minCoordinate;



//     // 2) init boundary conditions
//     // force value & location (in object's local coordinates)
//     std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> boundaryConditions;
//     for (std::map<int, Eigen::Vector3f>::const_iterator i = pBrb->contactTractions.cbegin();
//          i != pBrb->contactTractions.cend();
//          ++i) {

//         if (i->second.norm() > 0) {
//             int triangleID = i->first;
//             Eigen::Vector3d force = i->second.cast<double>();
//             Eigen::Vector3d location = pBrb->elemCtr.at(triangleID).cast<double>();
//             boundaryConditions.emplace_back(force, location - minCoordinate);
//         }
//     }
//     param.appliedForce = boundaryConditions;


//     if(boundaryConditions.empty())
//     {
//         printf("no boundary conditions, return\n");
//         return -1;
//     }


//     // 0) write the object in .obj file
//     std::string objectName = getName(pBrb->rigidBodyPtr).c_str();
//     trimesh::TriMesh* collisionMeshOut =  pBrb->collInUse.trimeshLowRes;
//     collisionMeshOut->write("original.obj");


//     // 5) run MPM simulator

//     std::cout << "run MPM simulator" << std::endl;
//     int timeStep = std::min((int)round(pBrb->contactDuration / param.dt), pBrb->m_params->mpmTimeStepsMax); // use the minimum value of calculated timestep and predefined timestep 
//     advanceStep(param, &particleVec, timeStep, pBrb->m_params->fractureSimDamageRatioThreshold, pBrb->m_params->forceScaleFactor);



//     // 6) run crack extraction algorithm
//     std::cout << "run crack extraction algorithm and refine mesh.." << std::endl;
//     std::string fielName = "particles" + std::to_string(timeStep) + ".txt";
//     std::ofstream outfile0(fielName, std::ios::trunc);
//     for (int f = 0; f < particleVec.size(); f++) 
//     {
//         outfile0 << std::scientific << std::setprecision(8) << mpmParticlePositions_Proj[f][0] << "  " << mpmParticlePositions_Proj[f][1] << " " << mpmParticlePositions_Proj[f][2] << " " << particleVec[f].Dp << " " << particleVec[f].Dg << "\n";
//     };
//     outfile0.close();
//     std::tuple<bool, meshFaceRefined, std::vector<meshObjFormat> > result = extractCrackSurface(fielName, &mpmParticlePositions_Proj, param);



//     // 7) get final fragments

//     if(std::get<0>(result) == false) // fail to use triangle.h or crack surface not found
//     {
//         return -1;
//     }
//     else
//     {
//         SCOPED_TIMER("addFractureDetail");

//         struct meshFaceRefined refinedMeshResult = std::get<1>(result);
//         chevronMarksGridAndMap chevron;
//         try 
//         {
//             // find chevron marks
//             struct meshFaceRefined refinedMeshResult = std::get<1>(result);
            
//             if (refinedMeshResult.verticesRefined.size() != 0) 
//             {
//                 std::cout << "calculate order information.." << std::endl;
//                 std::vector<Particle> particleVec_Copy = particleVec;
//                 for (int i = 0; i < particleVec_Copy.size(); i++) {
//                     particleVec_Copy[i].Dp = particleVec[i].Dg;
//                 }
//                 struct std::vector<verticesInformation> verticesVec = calculateOrderInfor(param, &particleVec_Copy, &refinedMeshResult);
//                 std::cout << "finish finding chevron " << std::endl;
//                 chevron = findChevronMarks(param, &verticesVec, &refinedMeshResult);
                
//             }


//         }
//         catch (...) 
//         {

//             std::cout<<"Fail to find chevron marks" <<std::endl;

//         }


//         std::cout << "start adding chevron " << std::endl;
//         // for each fragment, add chevron marks
//         std::vector<meshObjFormat> allFragments = std::get<2>(result);
//         std::string parentObjectName = getName(pBrb->rigidBodyPtr).c_str();
//         std::cout<<"allFragments.size() = " <<  allFragments.size() <<std::endl;
//         for(int hi = 0; hi < allFragments.size(); hi++)
//         {

//             triangulate_polygonMesh(&allFragments[hi]);
//             GenericMesh bareMesh;
//             bareMesh.m.vertices = allFragments[hi].vertices;
//             bareMesh.m.faces = allFragments[hi].faces;
//             for(int k = 0; k < bareMesh.m.faces.size(); k++)
//             {
//                 bareMesh.triangulatedFaces.push_back(bareMesh.m.faces[k][0]);
//                 bareMesh.triangulatedFaces.push_back(bareMesh.m.faces[k][1]);
//                 bareMesh.triangulatedFaces.push_back(bareMesh.m.faces[k][2]);
//             }
//             std::string bareMeshName = "bareMesh";
//             std::cout<<"Build bare mesh VDB grid.."<<std::endl;
//             openvdb::FloatGrid::Ptr bareMeshVdbGridPtr = crackSurfaceMeshToLevelSetGrid(bareMesh, param.drx, bareMeshName);
//             std::cout<<"Decimate bare mesh.."<<std::endl;
//             meshObjFormat fragmentVolume;
//             getSurfaceMeshFromVdb(bareMeshVdbGridPtr,100000, fragmentVolume);



//             std::vector<Eigen::Vector3d> faceNormal(fragmentVolume.faces.size());
//             std::vector< std::vector<int>> verticesFace(fragmentVolume.vertices.size());
//             for (int i = 0; i < fragmentVolume.faces.size(); ++i)
//             {
//                 int indexV1 = fragmentVolume.faces[i][0];
//                 Eigen::Vector3d v1 = fragmentVolume.vertices[indexV1];
//                 int indexV2 = fragmentVolume.faces[i][1];
//                 Eigen::Vector3d v2 = fragmentVolume.vertices[indexV2];
//                 int indexV3 = fragmentVolume.faces[i][2];
//                 Eigen::Vector3d v3 = fragmentVolume.vertices[indexV3];
//                 Eigen::Vector3d normal = (v1 - v2).cross(v2 - v3).normalized();
//                 faceNormal[i] = normal;

//                 verticesFace[indexV1].push_back(i);
//                 verticesFace[indexV2].push_back(i);
//                 verticesFace[indexV3].push_back(i);
//             }

// #pragma omp parallel for num_threads(param.numOfMpmThreads)
//             for(int j = 0; j < fragmentVolume.vertices.size(); j++)
//             {
//                 double value = calOrderValuePoint(fragmentVolume.vertices[j], param, param.smoothRadius, &chevron.gridMap, &chevron.nodesVec);
//                 Eigen::Vector3d pertNormal = {0,0,0};
//                 for(int m = 0; m < verticesFace[j].size(); m++)
//                 {
//                     pertNormal += faceNormal[verticesFace[j][m]];
//                 }
//                 fragmentVolume.vertices[j] += param.pertubationMagitude * value * pertNormal.normalized();

//                 fragmentVolume.vertices[j] += minCoordinate;
//             }


//             GenericMesh meshChevron;
//             meshChevron.m.vertices = fragmentVolume.vertices;
//             meshChevron.m.faces = fragmentVolume.faces;
//             for(int k = 0; k < meshChevron.m.faces.size(); k++)
//             {
//                 meshChevron.triangulatedFaces.push_back(meshChevron.m.faces[k][0]);
//                 meshChevron.triangulatedFaces.push_back(meshChevron.m.faces[k][1]);
//                 meshChevron.triangulatedFaces.push_back(meshChevron.m.faces[k][2]);
//             }
//             finalCrackSurfaceMesh->push_back(meshChevron);
//         }




//     }






//     // 7) make sure crack mesh is clean and return result
//     TIMESTACK_PUSH("cleanCrackMesh");

// #if 0
//     float t_max = maxTime, t_step; // time-steps
//     int addedCracks = 0, chk;
//     if (contactDuration < t_max) {
//         t_max = contactDuration
//     };

//     // controls how many timesteps to execute before we stop teh fracture sim
//     float maxPropagationLen = 1; // TODO
//     t_step = pBrb->getFractureTimeStepSize(maxPropagationSteps);
//     int steps = (int)((t_max / t_step) + 0.5); // round

//     printf("\n%% fracture simulation for %s_%d ... %d time steps", ((string*)rb->getUserPointer())->c_str(), rbTimeCode, steps);

//     if (steps == 0) {
//         printf("\n%% contact duration or max time too short (%.3lgs)", t_max);
//         return -1;
//     }

    

//     if (steps > 0) {
//         std::map<int, Eigen::Vector3d> q;130ncedContactTractions(q);

//         // TODO: Set simulation boundary conditions (force etc.)

//         bool done = false;
//         for (int k = 0; k < steps && !done; ++k) { // do fracture steps
 
//             printf("\n%% fracturing (%d/%d) ... ", k + 1, steps);

//             // TODO: do simulation, calculate stress, update damage etc,
            
//             printf("\n%% compute MPM timestep    ...");
//         }

//         // TODO: do mesh cutting and save mesh fragments
//     }
// #endif
//     TIMESTACK_POP();
//     return 1;



// }

