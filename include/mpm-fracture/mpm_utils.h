#pragma once
#ifdef _WIN32
// https://docs.microsoft.com/en-us/cpp/c-runtime-library/math-constants?redirectedfrom=MSDN&view=vs-2019
#define _USE_MATH_DEFINES
#include <cmath>
#endif




// eigen
#include <Eigen/Core>
#include <vector>

struct parametersSim {

	int numOfMpmThreads = 1;

	// computational domain
	Eigen::Vector3d length = { 1, 1, 1 }; // computation cube lengths of three dimensions (x , y , z). The origin point is (0 , 0 , 0)
	Eigen::Vector3d minCoordinate = { 0, 0, 0 }; // the minimum coordinate of the computation domain

	// Bcakground Eulerian grid
	double dx = 2E-2; // grid space
	double DP = dx * dx / 4;

	// openvdb voxel size
	double vdbVoxelSize = 0.0005;

	// material properties
	double dt = 1E-7; // timestep
	double dpx = dx / 2; // particle space
	double density = 2450; // particle density
	double vol = dpx * dpx; // each particle's initial volume
	double particle_mass = vol * density; // particle's mass
	double E = 3.2e10; // Young's modulus
	double nu = 0.2; //Poisson ratio
	double mu = E / (2 * (1 + nu)); // lame parameter mu / shear modulus  ::::only for isotropic material
	double lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter lambda  ::::only for isotropic material
	double K = 2 / 3 * mu + lambda; // bulk modulus  ::::only for isotropic material


	// applied force
	std::vector< std::pair<Eigen::Vector3d, Eigen::Vector3d> > appliedForce;



	// local damage field parameters
	double thetaf = 3.7e6;
	double Gf = 3;
	double lch = std::sqrt(3) * dx;
	double HsBar = thetaf * thetaf / 2 / E / Gf;
	double Hs = HsBar * lch / (1 - HsBar * lch);
	double damageThreshold = 0.97; // after this threshold, 
	double sigmoidK = 5; // this parameter control the curevature of the sigmoid function. It is recommend that it is bigger than 5



	// crack extraction algorithm parameter
	double dcx = dx; // crack grid spacing

	// the refined mesh resolution
	double drx = dx / 4.0; // must be samller than dpx


	// chevron marks parameters
	int pathLength = 2000; // the maximum path length
	double stepSize = drx; // stepsize
	double divertAngle = 0.3; // cos(angle) where angle is the angle between this step and next step
	double influenceRadius = 2 * drx; // the radius during which no vertices should be visited
	double smoothRadius = 2 * drx; // background grid spacing for smooth height value
	double pertubationMagitude = 4 * drx; // pertubation magnitude for each vertex
	double ratioOfDeepenMarks = 0.1; // ratio of chevron marks that should be deepened
	double deepenMagnitude = 2.0; // magnitude of deepened chevron marks


};

// return a matrix and a double value
struct matrixAndValue {
	Eigen::Matrix3d F; //matrix
	double S; //double value

	matrixAndValue(Eigen::Matrix3d iF, double iS)
		: F(iF)
		, S(iS)
	{
	}
};

// return a matrix and two double values
struct matrixAndTwoValues {
	Eigen::Matrix3d F; //matrix
	double Dp; //double value
	double Dg; //double value

	matrixAndTwoValues(Eigen::Matrix3d iF, double iDp, double iDg)
		: F(iF)
		, Dp(iDp)
		, Dg(iDg)
	{
	}
};


// return two matricies and a double value
struct TwoMatrixAndValue {
	Eigen::Matrix3d F; //matrix
	Eigen::Matrix3d Sigma; //matrix
	double S; //double value

	TwoMatrixAndValue(Eigen::Matrix3d iF, Eigen::Matrix3d iSi, double iS)
		: F(iF)
		, Sigma(iSi)
		, S(iS)
	{
	}
};


static int calculateID(int x, int y, int z, Eigen::Vector3d len, double dx) // coordinate of x and y, length in three dimensions of the cube, grid space
{
	Eigen::Vector3i length = (len / dx).cast<int>() + Eigen::Vector3i::Constant(1);
	int ID = z * (length(0) * length(1)) + y * length(0) + x;
	return ID;
};


