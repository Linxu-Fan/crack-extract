#ifndef DAMAGEGRADIENT_H

#define DAMAGEGRADIENT_H

#include <cmath>

#include "crackExtraction/particles.h"
#include "crackExtraction/utils.h"
#include "crackExtraction/weights.h"

// calculate the damage gradient of all particles and grid nodes.
void calDamageGradient(std::vector<Particle>*, parametersSim, double, std::map<int, int>*, std::vector<Grid>*);

// calculate the damage gradient of any give point
Eigen::Vector3d calDamageGradientPoint(Eigen::Vector3d, parametersSim, double, std::map<int, int>*, std::vector<Grid>*);

double calDamageValuePoint(Eigen::Vector3d pos, parametersSim param, double dx, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec);

#endif
