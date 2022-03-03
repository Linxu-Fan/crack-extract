#ifndef WEIGHTS_H

#define WEIGHTS_H

#include <math.h>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

// Struct of weights
struct weightAndDreri {
    Eigen::Vector3i ppIndex = { 0, 0, 0 }; // each particle's weight
    Eigen::Vector3d space = { 0, 0, 0 }; // each particle's weight derivative
    Eigen::MatrixXd weight;
    Eigen::MatrixXd deltaWeight;

    weightAndDreri(Eigen::Vector3i ippIndex, Eigen::Vector3d ispace, Eigen::MatrixXd iweight, Eigen::MatrixXd ideltaWeight)
        : ppIndex(ippIndex)
        , space(ispace)
        , weight(iweight)
        , deltaWeight(ideltaWeight)
    {
    }
};

struct weightAndDreri calWeight(double, Eigen::Vector3d);

#endif
