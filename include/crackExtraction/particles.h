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
	double Dp = 0; // particle's scalar damage value
	Eigen::Vector3d deltaD = {0 , 0 , 0}; // particle's damage gradient
	int color = 0;// boundary that traction is applied

	Particle(Eigen::Vector3d ipos, Eigen::Vector3d ivel, double im, int ic, double iDp) :
		pos(ipos),
		vel(ivel),
		m(im),
		color(ic),
		Dp(iDp){}

};


#endif













