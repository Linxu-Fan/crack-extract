#ifndef PARTICLES_H

#define PARTICLES_H

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iostream>
#include <set>
#include <iomanip>
#include <Eigen/Core>





// Struct of particles
struct Particle {
	Eigen::Vector3d pos = { 0 , 0 , 0}, vel = { 0 , 0 , 0}; // each particle's x and y position , velocity, and momentum
	double m = 0; // each particle's mass
	Eigen::Vector3i ppIndex = { 0 , 0 , 0}; // particle base index
	Eigen::MatrixXd weight;
	Eigen::MatrixXd deltaWeight;
	Eigen::Matrix3d F = Eigen::Matrix3d::Identity(); // each particle's deformation gradient
	Eigen::Matrix3d affine = Eigen::Matrix3d::Zero(); // each particle's affine term

	// non-local damage field method: deformation gradient in the previous timestep
	Eigen::Matrix3d pF = Eigen::Matrix3d::Identity(); // each particle's deformation gradient

	// local damage field method parameters
	double Dp = 0; // particle's scalar damage value
	Eigen::Vector3d deltaD = {0 , 0 , 0}; // particle's damage gradient

	// a scalar value that indicates the time of crack initialization
	double Dg = 0;

	// NACC parameters
	double alpha = -0.01;

	bool crackBoundary = false;
	int color = 0;// boundary that traction is applied

	int nearestPoint = -1000;
	Eigen::Vector3d crackSurfaceNormal = {0, 0, 0};


	// node index in the support radius of this particle. The order of the vector is important
	std::vector<int> supportNodes; // store the position of the node in vector "nodeVec" 
	std::vector<double> supportNodeWeight; // store the weight of node to the particle


	Eigen::Matrix3d cauchyStress = Eigen::Matrix3d::Zero(); // each particle's internal cauchy stress



	Particle(Eigen::Vector3d ipos, Eigen::Vector3d ivel, double im, int ic, double iDp) :
		pos(ipos),
		vel(ivel),
		m(im),
		color(ic),
		Dp(iDp){}

};


#endif













