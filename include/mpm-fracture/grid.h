#ifndef GRID_H

#define GRID_H

#include <Eigen/Core>
#include <math.h>
#include <vector>


// // Struct of grid
// struct Grid {
//     // property of velocity-field 0
//     double m = 0; // each node's mass
//     Eigen::Vector3d mom = { 0, 0, 0 }; // each node's momentum
//     Eigen::Vector3d velocity = { 0, 0, 0 }; // each node's velocity
//     Eigen::Vector3d force = { 0, 0, 0 }; // each node's force

//     // general grid node property
//     Eigen::Vector3i posIndex = { 0, 0, 0 };
//     Eigen::Vector3d deltaDi = { 0, 0, 0 }; // gradient of damage field
//     double Di = 0; // value of damage field
//     double sw = 0; // sum of particle-grid weight

//     // particle index in the support radius of this node. The order of the vector is important
// 	std::vector<int> supportParticles; // store the position of the particle in vector "particles"; 
// 	std::vector<double> supportParticlesWeight; // store the weight of particle to the grid node 



//     Grid(double im)
//         : m(im)
//     {
//     }
// };

#endif
