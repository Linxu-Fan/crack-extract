#include "mpm-fracture/damageGradient.h"
#include "mpm-fracture/weights.h"
#include "mpm-fracture/mpm_utils.h"

// calculate the damage gradient of all particles and grid nodes.
void calDamageGradient(std::vector<Particle>* particles, parametersSim param, double dx, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec)
{
    for (int f = 0; f < particles->size(); f++) {
        struct weightAndDreri WD = calWeight(dx, (*particles)[f].pos);

        (*particles)[f].ppIndex = WD.ppIndex;
        (*particles)[f].weight = WD.weight;
        (*particles)[f].deltaWeight = WD.deltaWeight;
    };

    int count = -1; // count the number of active grid node
    // number of grid nodes per edge

    // calculate node damage value
    for (int f = 0; f < particles->size(); f++) {
        (*particles)[f].deltaD = { 0, 0, 0 };

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    int ID = calculateID((*particles)[f].ppIndex[0] + i, (*particles)[f].ppIndex[1] + j, (*particles)[f].ppIndex[2] + k, param.length, dx);
                    double weight = (*particles)[f].weight(0, i) * (*particles)[f].weight(1, j) * (*particles)[f].weight(2, k);

                    if (weight != 0) {
                        if ((*gridMap).find(ID) == (*gridMap).end()) {
                            count += 1;
                            (*gridMap)[ID] = count;
                            (*nodesVec).push_back(Grid(0.0));

                            (*nodesVec)[count].posIndex = { (*particles)[f].ppIndex[0] + i, (*particles)[f].ppIndex[1] + j, (*particles)[f].ppIndex[2] + k };
                            (*nodesVec)[count].Di += (*particles)[f].Dp * weight;
                            (*nodesVec)[count].sw += weight;
                        } else {
                            int eid = (*gridMap)[ID];

                            (*nodesVec)[eid].Di += (*particles)[f].Dp * weight;
                            (*nodesVec)[eid].sw += weight;
                        };
                    };
                };
            };
        };
    };

    // calculate particle damage gradient
    for (int f = 0; f < particles->size(); f++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    int ID = calculateID((*particles)[f].ppIndex[0] + i, (*particles)[f].ppIndex[1] + j, (*particles)[f].ppIndex[2] + k, param.length, dx);
                    double weight = (*particles)[f].weight(0, i) * (*particles)[f].weight(1, j) * (*particles)[f].weight(2, k);

                    if (weight != 0) {
                        int eid = (*gridMap)[ID];

                        Eigen::Vector3d posD = (*particles)[f].pos - (*nodesVec)[eid].posIndex.cast<double>() * dx;
                        (*particles)[f].deltaD += weight / (dx * dx / 4) * (*nodesVec)[eid].Di / (*nodesVec)[eid].sw * posD;
                    };
                };
            };
        };
    };

    // calculate grid node's damage gradient. This gives the exact value.
    for (int g = 0; g < (*nodesVec).size(); g++) {
		Eigen::Vector3d weightVec = { 0.125, 0.75, 0.125 };
		Eigen::Vector3d posVec = { dx, 0, -dx };

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    int ID = calculateID((int)(*nodesVec)[g].posIndex[0] + i - 1, (int)(*nodesVec)[g].posIndex[1] + j - 1, (int)(*nodesVec)[g].posIndex[2] + k - 1, param.length, dx);
                    double weight = weightVec[i] * weightVec[j] * weightVec[k];

                    if ((*gridMap).find(ID) != (*gridMap).end()) {
                        int eid = (*gridMap)[ID];
						Eigen::Vector3d posD = { posVec[i], posVec[j], posVec[k] };
                        (*nodesVec)[g].deltaDi += weight / (dx * dx / 4) * (*nodesVec)[eid].Di / (*nodesVec)[eid].sw * posD;
                    };
                };
            };
        };

        (*nodesVec)[g].Di = (*nodesVec)[g].Di / (*nodesVec)[g].sw;
    };
};

// calculate the damage gradient of any give point
Eigen::Vector3d calDamageGradientPoint(Eigen::Vector3d pos, parametersSim param, double dx, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec)
{
    Eigen::Vector3d deltaPoint = { 0, 0, 0 };
    Eigen::Vector3d base = (pos) / dx - Eigen::Vector3d::Constant(0.5);
    Eigen::Vector3i ppIndex = base.cast<int>();

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                int ID = calculateID(ppIndex[0] + i, ppIndex[1] + j, ppIndex[2] + k, param.length, dx);

                struct weightAndDreri WD = calWeight(dx, (pos));
                Eigen::MatrixXd weightPoint = WD.weight;
                double weight = weightPoint(0, i) * weightPoint(1, j) * weightPoint(2, k);

                if ((*gridMap).find(ID) != (*gridMap).end()) {
                    int eid = (*gridMap)[ID];
                    Eigen::Vector3d posD = (pos) - (*nodesVec)[eid].posIndex.cast<double>() * dx;
                    deltaPoint += weight / (dx * dx / 4) * (*nodesVec)[eid].Di / (*nodesVec)[eid].sw * posD;
                };
            };
        };
    };

    return deltaPoint;
};


// calculate the damage value of any give point
double calDamageValuePoint(Eigen::Vector3d pos, parametersSim param, double dx, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec)
{
	double dpValue = 0;
	Eigen::Vector3d base = (pos) / dx - Eigen::Vector3d::Constant(0.5);
	Eigen::Vector3i ppIndex = base.cast<int>();

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				int ID = calculateID(ppIndex[0] + i, ppIndex[1] + j, ppIndex[2] + k, param.length, dx);

				struct weightAndDreri WD = calWeight(dx, (pos));
				Eigen::MatrixXd weightPoint = WD.weight;
				double weight = weightPoint(0, i) * weightPoint(1, j) * weightPoint(2, k);

				if ((*gridMap).find(ID) != (*gridMap).end()) {
					int eid = (*gridMap)[ID];
					Eigen::Vector3d posD = (pos)-(*nodesVec)[eid].posIndex.cast<double>() * dx;
					dpValue += weight  * (*nodesVec)[eid].Di;
				};
			};
		};
	};

	return dpValue;
};
