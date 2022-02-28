#ifndef DAMAGEGRADIENT_H

#define DAMAGEGRADIENT_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <set>
#include <time.h>


#include "mpm-fracture/mpm_utils.h"
#include "mpm-fracture/utils.h"
#include "mpm-fracture/particles.h"

// calculate the damage gradient of all particles and grid nodes.
void calDamageGradient(std::vector<Particle>*, parametersSim, double, std::map<int, int>*, std::vector<Grid>*);

// calculate the damage gradient of any give point
Eigen::Vector3d calDamageGradientPoint(Eigen::Vector3d, parametersSim, double, std::map<int, int>*, std::vector<Grid>*);

double calDamageValuePoint(Eigen::Vector3d pos, parametersSim param, double dx, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec);

#endif
