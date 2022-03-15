#include "crackExtraction/weights.h"
#include <iostream>

// calculate the weight and derivative of particles
struct weightAndDreri calWeight(double dx, Eigen::Vector3d pos)
{
    Eigen::Vector3d base = pos / dx - Eigen::Vector3d::Constant(0.5);
    Eigen::Vector3i ppIndex = base.cast<int>();
    Eigen::Vector3d space = pos / dx - ppIndex.cast<double>();

    /*std::cout << ppIndex[0] << " " << ppIndex[1] << " " << ppIndex[2] << " " << std::endl;*/

    // calculate weight
    Eigen::Vector3d col0 = 0.5 * (Eigen::Vector3d::Constant(1.5) - space).array().square();
    Eigen::Vector3d col1 = 0.75 - (Eigen::Vector3d::Constant(1.0) - space).array().square();
    Eigen::Vector3d col2 = 0.5 * (Eigen::Vector3d::Constant(0.5) - space).array().square();

    // calculate weight derivative
    Eigen::Vector3d deltaWeightcol0 = (space - Eigen::Vector3d::Constant(1.5)) / dx;
    Eigen::Vector3d deltaWeightcol1 = -2 * (space - Eigen::Vector3d::Constant(1.0)) / dx;
    Eigen::Vector3d deltaWeightcol2 = (space - Eigen::Vector3d::Constant(0.5)) / dx;

    Eigen::MatrixXd weight(3, 3);
    weight << col0, col1, col2;

    Eigen::MatrixXd deltaWeight(3, 3);
    deltaWeight << deltaWeightcol0, deltaWeightcol1, deltaWeightcol2;

    struct weightAndDreri res(ppIndex, space, weight, deltaWeight);

    return res;
};
