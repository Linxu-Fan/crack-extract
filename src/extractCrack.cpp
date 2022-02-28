#include "mpm-fracture/extractCrack.h"

using namespace voro;
using namespace Eigen;
using namespace std;


// read obj file
struct meshObjFormat readObj(std::string path)
{
	meshObjFormat result;

	// input mesh vertices and faces
	std::vector<Eigen::Vector3d> vertices;
	std::vector<std::vector<int>> faces;

	std::ifstream in;
	in.open(path);
	std::string line;
	while (getline(in, line))
	{
		if (line.size() > 0)
		{
			std::vector<std::string> vecCoor = split(line, " ");
            if(vecCoor.size() == 0)
            {
                std::cout<< "Obj mesh read error!"<<std::endl;
            }
			if (vecCoor[0] == "v")
			{
				Eigen::Vector3d vertex = { std::stod(vecCoor[1]) ,std::stod(vecCoor[2]) ,std::stod(vecCoor[3]) };
				vertices.push_back(vertex);
			}
			if (vecCoor[0] == "f")
			{
				std::vector<int> face;
                
				for (int k = 1; k < vecCoor.size(); k++)
				{
                    //std::cout<<k<<": " <<vecCoor[k]<<std::endl;
					face.push_back(std::stoi(vecCoor[k]) - 1);
				}
				faces.push_back(face);
			}
		}
	}
	in.close();



	result.faces = faces;
	result.vertices = vertices;


	return result;

}


// split a line from a text file
std::vector<std::string> split(const std::string& s, const std::string& seperator)
{
    std::vector<std::string> result;
    typedef std::string::size_type string_size;
    string_size i = 0;

    while (i != s.size()) {
        int flag = 0;
        while (i != s.size() && flag == 0) {
            flag = 1;
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[i] == seperator[x]) {
                    ++i;
                    flag = 0;
                    break;
                }
        }

        flag = 0;
        string_size j = i;
        while (j != s.size() && flag == 0) {
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[j] == seperator[x]) {
                    flag = 1;
                    break;
                }
            if (flag == 0)
                ++j;
        }
        if (i != j) {
            result.push_back(s.substr(i, j - i));
            i = j;
        }
    }
    return result;
}

// Set the damage phase of a grid node vector into a specific value
void setNodeValue(std::vector<Grid>* nodesVec, int va)
{
    int numDamage = (*nodesVec).size();
    for (int k = 0; k < numDamage; k++) {
        (*nodesVec)[k].Di = va;
    }
}

// Find the bounding box boundary nodes and set its damage phase into a specific value
void findBoundaryNodes(std::vector<Particle>* particles, std::vector<Grid>* nodesVec, std::map<int, int>* gridMap, struct parametersSim parti, int va)
{

    int count1 = (*nodesVec).size();
    int numDamage = (*nodesVec).size();
    for (int m = 0; m < numDamage; m++) {
        for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
                for (int k = -1; k < 2; k++) {

                    int ID = calculateID((*nodesVec)[m].posIndex[0] + i, (*nodesVec)[m].posIndex[1] + j, (*nodesVec)[m].posIndex[2] + k, parti.length, parti.dx);
                    if ((*gridMap).find(ID) == (*gridMap).end()) {
                        (*gridMap)[ID] = count1;
                        (*nodesVec).push_back(Grid(0.0));

                        (*nodesVec)[count1].posIndex = { (*nodesVec)[m].posIndex[0] + i, (*nodesVec)[m].posIndex[1] + j, (*nodesVec)[m].posIndex[2] + k };
                        (*nodesVec)[count1].Di = va;

                        count1 += 1;
                    }
                }
            }
        }
    }
}

// Read particles' positions and damage phases
void readParticles( std::vector<Particle>* particlesRaw, std::vector<Particle>* particles, bool ifFully, struct parametersSim param)
{

    for(int i =0; i < (*particlesRaw).size(); i++)
    {
        if (ifFully == true) 
        {
            if ((*particlesRaw)[i].Dp >= param.damageThreshold) 
            {
                (*particles).push_back(Particle((*particlesRaw)[i].pos, (*particlesRaw)[i].vel, 0, 0, 1));
            }
        }
        else 
        {
            (*particles).push_back(Particle((*particlesRaw)[i].pos, (*particlesRaw)[i].vel, 0, 0, 1));
        }
    }
}

// Calculate the damage value of any point and return the value
double ifFullyDamaged(Eigen::Vector3d pos, parametersSim param, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec)
{
    double damageValue = 0;
    Eigen::Vector3d base = pos / param.dx - Eigen::Vector3d::Constant(0.5);
    Eigen::Vector3i ppIndex = base.cast<int>();

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                int ID = calculateID(ppIndex[0] + i, ppIndex[1] + j, ppIndex[2] + k, param.length, param.dx);
                struct weightAndDreri WD = calWeight(param.dx, pos);
                MatrixXd weightPoint = WD.weight;
                double weight = weightPoint(0, i) * weightPoint(1, j) * weightPoint(2, k);

                if ((*gridMap).find(ID) != (*gridMap).end()) {
                    int eid = (*gridMap)[ID];
                    damageValue += weight * (*nodesVec)[eid].Di;
                }
            };
        };
    };

    return damageValue;
}

// Store the index of each vertex and its position
int findIndexVertex(Eigen::Vector3d pos, std::vector<Eigen::Vector3d>* vertexIndex)
{
    int index = vertexIndex->size();
    (*vertexIndex).push_back(pos);

    return index;
}

// Read voronoi tesselation
void readPoints(std::vector<Point>* points)
{
    ifstream in;
    in.open("crack.vol");
    std::string line;
    while (getline(in, line)) {
        stringstream ss(line);
        std::string element;
        std::vector<std::string> eachLine;

        while (ss >> element) {
            eachLine.push_back(element);
        }

        if (eachLine.size() > 0) {
            int index = std::stoi(eachLine[0]);
            Eigen::Vector3d pos = { std::atof(eachLine[1].c_str()), std::atof(eachLine[2].c_str()), std::atof(eachLine[3].c_str()) };
            int numVertices = std::stoi(eachLine[4]);

            // read all vertices in the cell of this point
            std::vector<Eigen::Vector3d> verticsCoor;
            for (int i = 5; i < 5 + numVertices; i++) {
                std::vector<std::string> vecCoor = split(eachLine[i], ",()");

                Eigen::Vector3d coor = { std::atof(vecCoor[0].c_str()), std::atof(vecCoor[1].c_str()), std::atof(vecCoor[2].c_str()) };
                verticsCoor.push_back(coor);
            }

            int indnum = 5 + numVertices;
            int numFaces = std::stoi(eachLine[indnum]);

            // read all faces in the cell of this point
            std::vector<std::vector<int>> verticesFace;
            for (int i = 6 + numVertices; i < 6 + numVertices + numFaces; i++) {
                std::vector<std::string> vecFace = split(eachLine[i], ",()");

                std::vector<int> eachFace;
                for (int j = 0; j < vecFace.size(); j++) {
                    eachFace.push_back(std::stoi(vecFace[j]));
                }
                verticesFace.push_back(eachFace);
            }

            // read all surfaceNormal in the cell of this point
            std::vector<Eigen::Vector3d> surfaceNormal;
            for (int i = 6 + numVertices + numFaces; i < 6 + numVertices + numFaces * 2; i++) {
                std::vector<std::string> vecFaceNormal = split(eachLine[i], ",()");

                Eigen::Vector3d eachFaceNormal = { std::atof(vecFaceNormal[0].c_str()), std::atof(vecFaceNormal[1].c_str()), std::atof(vecFaceNormal[2].c_str()) };
                surfaceNormal.push_back(eachFaceNormal);
            }

            // read all neighbours
            std::vector<int> neighbour;
            for (int k = 6 + numVertices + numFaces * 2; k < eachLine.size(); k++) {
                neighbour.push_back(std::stoi(eachLine[k]));
            }

            (*points).push_back(Point(index, pos, numVertices, verticsCoor, numFaces, verticesFace, surfaceNormal, neighbour));

        } else {
            Eigen::Vector3d pos = { 0, 0, 0 }; // each point's position
            std::vector<Eigen::Vector3d> verticsCoor; // vertices' coordinate
            std::vector<std::vector<int>> verticesFace; // vertices of each face
            std::vector<Eigen::Vector3d> surfaceNormal; // vertices' coordinate
            std::vector<int> neighbour; // neighbour points that share common faces with this point
            std::vector<int> neighbourCalculated; // neighbour points that have already find shared face
            (*points).push_back(Point(-999, pos, 0, verticsCoor, 0, verticesFace, surfaceNormal, neighbour));
        }
    }
    in.close();
}

// Read all structured nodes and calculate the damage gradient
void readParticlesAndCalGradient(std::vector<Grid>* fullyDamagedParticlesNodesVec, std::vector<Particle>* particles, parametersSim param, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec)
{
    for (int i = 0; i < (*fullyDamagedParticlesNodesVec).size(); i++) {
        if ((*fullyDamagedParticlesNodesVec)[i].Di == 1) 
        {
            Eigen::Vector3d ipos = { (*fullyDamagedParticlesNodesVec)[i].posIndex[0] * param.dx, (*fullyDamagedParticlesNodesVec)[i].posIndex[1] * param.dx, (*fullyDamagedParticlesNodesVec)[i].posIndex[2] * param.dx };
            Eigen::Vector3d ivel = { 0, 0, 0 };
            (*particles).push_back(Particle(ipos, ivel, 0, 0, 1));
        }

        if ((*fullyDamagedParticlesNodesVec)[i].Di == 2) {
            Eigen::Vector3d ipos = { (*fullyDamagedParticlesNodesVec)[i].posIndex[0] * param.dx, (*fullyDamagedParticlesNodesVec)[i].posIndex[1] * param.dx, (*fullyDamagedParticlesNodesVec)[i].posIndex[2] * param.dx };
            Eigen::Vector3d ivel = { 0, 0, 0 };
            (*particles).push_back(Particle(ipos, ivel, 0, 0, 0));
        }
    }


    calDamageGradient(particles, param, param.dx, gridMap, nodesVec);
}

// Find paths between two nodes
bool findPath(Eigen::Vector3d startNode, Eigen::Vector3d stopNode, parametersSim param, std::vector<Grid>* fullyDamagedParticlesNodesVec, std::map<int, int>* fullyDamagedParticlesGridMap, std::vector<int>* surfaceNodesID)
{

    Eigen::Vector3d startNodePos = startNode / param.dx;
    Eigen::Vector3i startNodeIndex;
    startNodeIndex[0] = round(startNodePos[0]);
    startNodeIndex[1] = round(startNodePos[1]);
    startNodeIndex[2] = round(startNodePos[2]);
    int startNodeID = calculateID(startNodeIndex[0], startNodeIndex[1], startNodeIndex[2], param.length, param.dx);

    Eigen::Vector3d stopNodePos = stopNode / param.dx;
    Eigen::Vector3i stopNodeIndex;
    stopNodeIndex[0] = round(stopNodePos[0]);
    stopNodeIndex[1] = round(stopNodePos[1]);
    stopNodeIndex[2] = round(stopNodePos[2]);
    int stopNodeID = calculateID(stopNodeIndex[0], stopNodeIndex[1], stopNodeIndex[2], param.length, param.dx);

    Eigen::Vector3i boundingBox = stopNodeIndex - startNodeIndex;
    int maxAllowedSteps = abs(boundingBox[0]) + abs(boundingBox[1]) + abs(boundingBox[2]);

    // if a particle has no neighbours
    bool noNeighbours = false;

    vector<Eigen::Vector3i> visitQueueIndex; // stores node posIndexs that should be visited
    std::map<int, int> IDMap; // the key is the ID of a node, the value is the node's position in std::vector nodeNeighbour
    std::vector<std::vector<int>> nodeParent; // stores parent-layer neighbouring information of nodes

    // Initialize all vectors with starting node
    visitQueueIndex.push_back(startNodeIndex);
    IDMap[startNodeID] = 0;
    std::vector<int> startNodeParent;
    startNodeParent.push_back(-99);
    nodeParent.push_back(startNodeParent); // start node has no parent-layer, so give a negative value

    bool reachMaxmiumBoundingBox = false; // if the expansion layer reach the maximum bounding box
    bool reachStopNode = false;
    int lastLayerPos = 0, lastLayerLength = 1;

    //////////////
    int depth = 0; // the number of step used

    do {

        std::vector<Eigen::Vector3i> visitQueueIndexLayer;
        std::vector<std::vector<int>> nodeParentLayer;
        int lengthOfALayer = 0;

        for (int k = lastLayerPos; k < lastLayerPos + lastLayerLength; k++) {
            Eigen::Vector3i parentNodeIndex = visitQueueIndex[k]; // the node that is going to be visited
            int parentNodeID = calculateID(parentNodeIndex[0], parentNodeIndex[1], parentNodeIndex[2], param.length, param.dx);

            for (int axis = 0; axis < 3; axis++) // three axis directions
            {
                for (int pn = 0; pn < 2; pn++) {
                    Eigen::Vector3i normal = { 0, 0, 0 };
                    normal[axis] = 2 * pn - 1;
                    Eigen::Vector3i neighbourNodePosIndex = parentNodeIndex + normal;

                    int neighbourNodeID = calculateID(neighbourNodePosIndex[0], neighbourNodePosIndex[1], neighbourNodePosIndex[2], param.length, param.dx);
                    if (IDMap.find(neighbourNodeID) == IDMap.end()) // this node is not in the queue
                    {
                        if ((*fullyDamagedParticlesGridMap).find(neighbourNodeID) != (*fullyDamagedParticlesGridMap).end()) {
                            int eid = (*fullyDamagedParticlesGridMap)[neighbourNodeID];
                            if ((*fullyDamagedParticlesNodesVec)[eid].Di == 2) // if this node is a boundary shell
                            {
                                visitQueueIndexLayer.push_back(neighbourNodePosIndex);
                                IDMap[neighbourNodeID] = lengthOfALayer + lastLayerPos + lastLayerLength;

                                // find the neghbours of this node in the last layer
                                std::vector<int> aSingleNodeParent;
                                for (int nt = lastLayerPos; nt < lastLayerPos + lastLayerLength; nt++) {
                                    Eigen::Vector3i lastLayerNode = visitQueueIndex[nt];
                                    Eigen::Vector3i diffIndex = lastLayerNode - neighbourNodePosIndex;
                                    int sumDiff = abs(diffIndex[0]) + abs(diffIndex[1]) + abs(diffIndex[2]);
                                    if (sumDiff == 1) {
                                        int lastLayerNodeID = calculateID(lastLayerNode[0], lastLayerNode[1], lastLayerNode[2], param.length, param.dx);
                                        aSingleNodeParent.push_back(lastLayerNodeID);
                                    }
                                }
                                nodeParentLayer.push_back(aSingleNodeParent);
                                lengthOfALayer += 1;

                                if (neighbourNodeID == stopNodeID) {
                                    reachStopNode = true;
                                }
                            }
                        }
                    }
                }
            }
        }
        lastLayerPos = visitQueueIndex.size();
        lastLayerLength = lengthOfALayer;
        visitQueueIndex.insert(visitQueueIndex.end(), visitQueueIndexLayer.begin(), visitQueueIndexLayer.end());
        nodeParent.insert(nodeParent.end(), nodeParentLayer.begin(), nodeParentLayer.end());

        depth += 1;

        if (lastLayerLength == 0) {
            return true;
        }

        if (depth >= round(1.3 * (double)maxAllowedSteps) && reachStopNode == false) {
            return true;
        }

    } while (reachStopNode == false && noNeighbours == false);

    //// Finding all possible paths
    std::vector<std::vector<int>> paths;
    bool reachStartNode = false;

    std::vector<int> initializeStop;
    initializeStop.push_back(stopNodeID);
    paths.push_back(initializeStop);
    int pathSize = paths.size();

    do {
        int nextLayerSize = 0;
        for (int k = 0; k < pathSize; k++) {
            std::vector<int> partialPath = paths[k];
            std::vector<int> lastNodeParents = nodeParent[IDMap[paths[k].back()]];

            for (int m = 0; m < lastNodeParents.size(); m++) {

                if (count((*surfaceNodesID).begin(), (*surfaceNodesID).end(), lastNodeParents[m]) == 0) {
                    //std::vector<int> partialPathExtend = partialPath;
                    std::vector<int> partialPathExtend;
                    partialPathExtend.push_back(lastNodeParents[m]);

                    if (lastNodeParents[m] == startNodeID) {
                        reachStartNode = true;
                        return false;
                    }

                    paths.push_back(partialPathExtend);
                    nextLayerSize += 1;
                }
            }
        }

        paths.erase(paths.begin(), paths.begin() + pathSize);
        pathSize = nextLayerSize;
        //cout <<"pathSize = "<< pathSize << endl;

    } while (reachStartNode == false && pathSize != 0);

    if (reachStartNode == false) {
        return true;
    }
}

// Find if a pair of nodes belong to critical nodes. The function return true if one node is a critical node
bool ifCriticalNode(Eigen::Vector3d node1, parametersSim param, std::vector<int>* criticalNodeIndex)
{

    Eigen::Vector3d node1Pos = node1 / param.dx;
    Eigen::Vector3i node1Index;
    node1Index[0] = round(node1Pos[0]);
    node1Index[1] = round(node1Pos[1]);
    node1Index[2] = round(node1Pos[2]);
    int node1ID = calculateID(node1Index[0], node1Index[1], node1Index[2], param.length, param.dx);

    if (count((*criticalNodeIndex).begin(), (*criticalNodeIndex).end(), node1ID) != 0) {
        return true;
    } else {
        return false;
    }
}

// Find the nearest boundary node of a critical node
Eigen::Vector3i findNearestBoundaryNode(int nodeIDPoint, std::vector<Point>* points, std::vector<int>* boundaryNodesID, std::vector<Eigen::Vector3i>* boundaryNodesPosIndex, std::map<int, int>* pointIndexFind, parametersSim param, std::vector<int>* criticalNodeIndex)
{
    Eigen::Vector3i nodeIndex = (*boundaryNodesPosIndex)[nodeIDPoint];
    int nodeID = (*boundaryNodesID)[nodeIDPoint];

    bool reachNearest = false;

    std::vector<Eigen::Vector3i> visitQueueIndex; // stores node posIndexs that should be visited
    std::map<int, int> IDMap; // the key is the ID of a node, the value is the node's position in std::vector nodeNeighbour
    std::vector<std::vector<int>> nodeParent; // stores parent-layer neighbouring information of nodes

    // Initialize all vectors with starting node
    visitQueueIndex.push_back(nodeIndex);
    IDMap[nodeID] = 0;
    std::vector<int> startNodeParent;
    startNodeParent.push_back(-99);
    nodeParent.push_back(startNodeParent); // start node has no parent-layer, so give a negative value

    bool reachMaxmiumBoundingBox = false; // if the expansion layer reach the maximum bounding box
    bool reachStopNode = false;
    int lastLayerPos = 0, lastLayerLength = 1;

    do {

        std::vector<Eigen::Vector3i> visitQueueIndexLayer;
        std::vector<std::vector<int>> nodeParentLayer;
        int lengthOfALayer = 0;

        for (int k = lastLayerPos; k < lastLayerPos + lastLayerLength; k++) {
            Eigen::Vector3i parentNodeIndex = visitQueueIndex[k]; // the node that is going to be visited
            int parentNodeID = calculateID(parentNodeIndex[0], parentNodeIndex[1], parentNodeIndex[2], param.length, param.dx);

            for (int axis = 0; axis < 3; axis++) // three axis directions
            {
                for (int pn = 0; pn < 2; pn++) {
                    Eigen::Vector3i normal = { 0, 0, 0 };
                    normal[axis] = 2 * pn - 1;
                    Eigen::Vector3i neighbourNodePosIndex = parentNodeIndex + normal;

                    int neighbourNodeID = calculateID(neighbourNodePosIndex[0], neighbourNodePosIndex[1], neighbourNodePosIndex[2], param.length, param.dx);
                    if (IDMap.find(neighbourNodeID) == IDMap.end()) // this node is not in the queue
                    {

                        if (count((*boundaryNodesID).begin(), (*boundaryNodesID).end(), neighbourNodeID) != 0) {
                            visitQueueIndexLayer.push_back(neighbourNodePosIndex);
                            IDMap[neighbourNodeID] = lengthOfALayer + lastLayerPos + lastLayerLength;

                            // find the neghbours of this node in the last layer
                            std::vector<int> aSingleNodeParent;
                            for (int nt = lastLayerPos; nt < lastLayerPos + lastLayerLength; nt++) {
                                Eigen::Vector3i lastLayerNode = visitQueueIndex[nt];
                                Eigen::Vector3i diffIndex = lastLayerNode - neighbourNodePosIndex;
                                int sumDiff = abs(diffIndex[0]) + abs(diffIndex[1]) + abs(diffIndex[2]);
                                if (sumDiff == 1) {
                                    int lastLayerNodeID = calculateID(lastLayerNode[0], lastLayerNode[1], lastLayerNode[2], param.length, param.dx);
                                    aSingleNodeParent.push_back(lastLayerNodeID);
                                }
                            }
                            nodeParentLayer.push_back(aSingleNodeParent);
                            lengthOfALayer += 1;

                            if (ifCriticalNode(neighbourNodePosIndex.cast<double>() * param.dx, param, criticalNodeIndex) == false) {
                                return neighbourNodePosIndex;
                            }
                        }
                    }
                }
            }
        }
        lastLayerPos = visitQueueIndex.size();
        lastLayerLength = lengthOfALayer;
        visitQueueIndex.insert(visitQueueIndex.end(), visitQueueIndexLayer.begin(), visitQueueIndexLayer.end());
        nodeParent.insert(nodeParent.end(), nodeParentLayer.begin(), nodeParentLayer.end());

        // if this point is a desolate critical point, return a negative value
        if (lastLayerLength == 0) {
            Eigen::Vector3i noConnection = { -999, 0, 0 };
            return noConnection;
        }

    } while (reachNearest == false);
}

// Judge if a pair of points are on different sides of a crack
bool ifTwoSides(int startNode, int stopNode, std::vector<Point>* points, std::vector<int>* boundaryNodesID, std::vector<Eigen::Vector3i>* boundaryNodesPosIndex, std::map<int, int>* pointIndexFind, parametersSim param, std::vector<int>* criticalNodeIndex)
{

    Eigen::Vector3i startNodeIndex = (*boundaryNodesPosIndex)[startNode];
    startNodeIndex = findNearestBoundaryNode(startNode, points, boundaryNodesID, boundaryNodesPosIndex, pointIndexFind, param, criticalNodeIndex);
    // if this point is a desolate critical point, return true. Keep this face though it may become a tooth.
    if (startNodeIndex[0] < 0) {
        return true;
    }
    int startNodeID = calculateID(startNodeIndex[0], startNodeIndex[1], startNodeIndex[2], param.length, param.dx);
    int startNodePointIndex = (*pointIndexFind)[startNodeID];

    Eigen::Vector3i stopNodeIndex = (*boundaryNodesPosIndex)[stopNode];
    if (ifCriticalNode(stopNodeIndex.cast<double>() * param.dx, param, criticalNodeIndex) == true) {
        stopNodeIndex = findNearestBoundaryNode(stopNode, points, boundaryNodesID, boundaryNodesPosIndex, pointIndexFind, param, criticalNodeIndex);
        // if this point is a desolate critical point, return true. Keep this face though it may become a tooth.
        if (stopNodeIndex[0] < 0) {
            return true;
        }
    }
    int stopNodeID = calculateID(stopNodeIndex[0], stopNodeIndex[1], stopNodeIndex[2], param.length, param.dx);
    int stopNodePointIndex = (*pointIndexFind)[stopNodeID];

    //cout << "Real start and stop points are: " << endl;
    //cout << "startNodeIndex = "<< startNodePointIndex <<" pos "<< startNodeIndex[0] << " " << startNodeIndex[1] << " " << startNodeIndex[2] << " " << endl;
    //cout << "stopNodeIndex = " << stopNodePointIndex << " pos " << stopNodeIndex[0] << " " << stopNodeIndex[1] << " " << stopNodeIndex[2] << " " << endl<<endl;

    std::vector<int> sameSideNeighbours; // the same side neighbours
    std::vector<int> otherSideNeighbours; // the other side neighbours

    // find initial cube bounding neighbours of startNode
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                if (i + j + k != 0) {
                    Eigen::Vector3i normal = { i, j, k };
                    Eigen::Vector3i neighbourNodePosIndex = startNodeIndex + normal;

                    int neighbourNodeID = calculateID(neighbourNodePosIndex[0], neighbourNodePosIndex[1], neighbourNodePosIndex[2], param.length, param.dx);

                    if (count((*boundaryNodesID).begin(), (*boundaryNodesID).end(), neighbourNodeID) != 0) {
                        int neighbourNodePointIndex = (*pointIndexFind)[neighbourNodeID];
                        sameSideNeighbours.push_back(neighbourNodePointIndex);
                    }
                }
            }
        }
    }

    // find the same side neighbours of this point
    for (int f = 0; f < (*points)[startNodePointIndex].neighbourSameSide.size(); f++) {
        int sameSideNeighbourPoint = (*points)[startNodePointIndex].neighbourSameSide[f];
        if (count(sameSideNeighbours.begin(), sameSideNeighbours.end(), sameSideNeighbourPoint) == 0) {
            sameSideNeighbours.push_back(sameSideNeighbourPoint);
        }
    }

    // find the other side neighbours of this point
    for (int f = 0; f < (*points)[startNodePointIndex].neighbourOtherSide.size(); f++) {
        int otherSideNeighbourPoint = (*points)[startNodePointIndex].neighbourOtherSide[f];
        if (count(otherSideNeighbours.begin(), otherSideNeighbours.end(), otherSideNeighbourPoint) == 0) {
            otherSideNeighbours.push_back(otherSideNeighbourPoint);
        }
    }

    //cout << "Initial same side neighbours are: ";
    //for (int i = 0; i < sameSideNeighbours.size(); i++)
    //{
    //	cout << sameSideNeighbours[i] << " ";
    //}

    //cout <<endl<< "Initial the other side neighbours are: ";
    //for (int i = 0; i < otherSideNeighbours.size(); i++)
    //{
    //	cout << otherSideNeighbours[i] << " ";
    //}
    //cout << endl << endl << endl << endl;

    bool findStopNode = false;
    int sameSideALayerStart = 0, sameSideALayerLength = sameSideNeighbours.size();
    int otherSideALayerStart = 0, otherSideALayerLength = otherSideNeighbours.size();
    do {

        int sameSideNextLayerLength = 0;
        int otherSideNextLayerLength = 0;

        for (int s = sameSideALayerStart; s < sameSideALayerStart + sameSideALayerLength; s++) {
            int nodePointIndex = sameSideNeighbours[s];

            // find the same side neighbours
            for (int k = 0; k < (*points)[nodePointIndex].neighbourSameSide.size(); k++) {
                int sameSideNeighbourID = (*points)[nodePointIndex].neighbourSameSide[k];

                if (sameSideNeighbourID == stopNodePointIndex) {
                    return false;
                }
                if (count(sameSideNeighbours.begin(), sameSideNeighbours.end(), sameSideNeighbourID) == 0) {
                    sameSideNeighbours.push_back(sameSideNeighbourID);
                    sameSideNextLayerLength += 1;
                }
            }

            // find the other side neighbours
            for (int m = 0; m < (*points)[nodePointIndex].neighbourOtherSide.size(); m++) {
                int otherSideNeighbourID = (*points)[nodePointIndex].neighbourOtherSide[m];

                if (otherSideNeighbourID == stopNodePointIndex) {

                    return true;
                }
                if (count(otherSideNeighbours.begin(), otherSideNeighbours.end(), otherSideNeighbourID) == 0) {
                    //cout << "nodePointIndex = " << nodePointIndex << "; otherSideNeighbourID = " << otherSideNeighbourID << endl;
                    otherSideNeighbours.push_back(otherSideNeighbourID);
                    //otherSideNextLayerLength += 1;
                    otherSideALayerLength += 1;
                }
            }
        }

        for (int s = otherSideALayerStart; s < otherSideALayerStart + otherSideALayerLength; s++) {

            int nodePointIndex = otherSideNeighbours[s];
            // find the same side neighbours; "The same" means in the other crack side
            for (int k = 0; k < (*points)[nodePointIndex].neighbourSameSide.size(); k++) {
                int sameSideNeighbourID = (*points)[nodePointIndex].neighbourSameSide[k];
                if (sameSideNeighbourID == stopNodePointIndex) {
                    return true;
                }
                if (count(otherSideNeighbours.begin(), otherSideNeighbours.end(), sameSideNeighbourID) == 0) {
                    otherSideNeighbours.push_back(sameSideNeighbourID);
                    otherSideNextLayerLength += 1;
                }
            }

            //// find the other side neighbours; "The other" means in the same side of startNode
            //for (int m = 0; m < (*points)[nodePointIndex].neighbourOtherSide.size(); m++)
            //{
            //	int otherSideNeighbourID = (*points)[nodePointIndex].neighbourOtherSide[m];
            //	if (otherSideNeighbourID == stopNodePointIndex)
            //	{
            //		return false;
            //	}
            //	if (count(sameSideNeighbours.begin(), sameSideNeighbours.end(), otherSideNeighbourID) == 0)
            //	{
            //		sameSideNeighbours.push_back(otherSideNeighbourID);
            //		sameSideNextLayerLength += 1;
            //	}

            //}
        }

        sameSideALayerStart = sameSideALayerStart + sameSideALayerLength;
        sameSideALayerLength = sameSideNextLayerLength;
        otherSideALayerStart = otherSideALayerStart + otherSideALayerLength;
        otherSideALayerLength = otherSideNextLayerLength;

        if (sameSideNextLayerLength + sameSideNextLayerLength == 0) {
            return true;
        }

    } while (findStopNode == false);
}

// Extract the crack surface
std::tuple<bool, meshObjFormat, meshObjFormat, std::vector<meshObjFormat> >  extractCrackSurface(std::vector<Particle>* particlesRaw, struct parametersSim param)
{
    SCOPED_TIMER(__FUNCTION__);
    
    clock_t startTime, endTime;
    startTime = clock();

    cout << "Start extracting" << endl;

    //*********Read fully damaged particles***********//
    std::vector<Particle> fullyDamagedParticles;
    std::map<int, int> fullyDamagedParticlesGridMap;
    std::vector<Grid> fullyDamagedParticlesNodesVec;

    readParticles(particlesRaw, &fullyDamagedParticles, true, param);
    calDamageGradient(&fullyDamagedParticles, param, param.dx, &fullyDamagedParticlesGridMap, &fullyDamagedParticlesNodesVec);
    setNodeValue(&fullyDamagedParticlesNodesVec, 1);
    findBoundaryNodes(&fullyDamagedParticles, &fullyDamagedParticlesNodesVec, &fullyDamagedParticlesGridMap, param, 2); // This bounding box is used to generate boundary nodes
    findBoundaryNodes(&fullyDamagedParticles, &fullyDamagedParticlesNodesVec, &fullyDamagedParticlesGridMap, param, 3); // This bounding box is used to generate clip surface mesh
    //*********Read fully damaged particles***********//


	cout << "fullyDamagedParticles.size() = " << fullyDamagedParticles.size() << endl;
	cout << "fullyDamagedParticlesNodesVec.size() = " << fullyDamagedParticlesNodesVec.size() << endl;



    //*********Read all particles***********//
    std::vector<Particle> allParticles;
    std::map<int, int> allParticlesGridMap;
    std::vector<Grid> allParticlesNodesVec;

    readParticles(particlesRaw, &allParticles, false, param);
    calDamageGradient(&allParticles, param, param.dx, &allParticlesGridMap, &allParticlesNodesVec);
    setNodeValue(&allParticlesNodesVec, 1);
    findBoundaryNodes(&allParticles, &allParticlesNodesVec, &allParticlesGridMap, param, 2); // This bounding box is used to generate boundary nodes
    findBoundaryNodes(&allParticles, &allParticlesNodesVec, &allParticlesGridMap, param, 3); // This bounding box is used to generate clip surface mesh
    //*********Read all particles***********//


	cout << "allParticles.size() = " << allParticles.size() << endl;
	cout << "allParticlesNodesVec.size() = " << allParticlesNodesVec.size() << endl;


    // find boundary nodes
    std::vector<int> allParticlesNodeBoundaryIndex; // store the index or ID of allParticles boundary nodes
    for (int i = 0; i < allParticlesNodesVec.size(); i++) {
        if (allParticlesNodesVec[i].Di == 2) {
            int ID = calculateID(allParticlesNodesVec[i].posIndex[0], allParticlesNodesVec[i].posIndex[1], allParticlesNodesVec[i].posIndex[2], param.length, param.dx);
            allParticlesNodeBoundaryIndex.push_back(ID);
        }
    }


    //************find surface nodes//
    std::vector<int> surfaceNodesID; // store surface nodes ID
    std::vector<Eigen::Vector3i> surfaceNodesPosIndex; // store surface position index
    std::vector<int> boundaryNodesIDNoClean; // store boundary nodes ID
    std::vector<Eigen::Vector3i> boundaryNodesPosIndexNoClean; // store boundary position index
    for (int i = 0; i < fullyDamagedParticlesNodesVec.size(); i++) {
        if (fullyDamagedParticlesNodesVec[i].Di == 2) {
            int ID = calculateID(fullyDamagedParticlesNodesVec[i].posIndex[0], fullyDamagedParticlesNodesVec[i].posIndex[1], fullyDamagedParticlesNodesVec[i].posIndex[2], param.length, param.dx);
            if (count(allParticlesNodeBoundaryIndex.begin(), allParticlesNodeBoundaryIndex.end(), ID) == 0) {
                boundaryNodesPosIndexNoClean.push_back(fullyDamagedParticlesNodesVec[i].posIndex);
                boundaryNodesIDNoClean.push_back(ID);
            }

            //************find surface nodes//
            if (count(allParticlesNodeBoundaryIndex.begin(), allParticlesNodeBoundaryIndex.end(), ID) != 0) {
                int surfaceNodeID = calculateID(fullyDamagedParticlesNodesVec[i].posIndex[0], fullyDamagedParticlesNodesVec[i].posIndex[1], fullyDamagedParticlesNodesVec[i].posIndex[2], param.length, param.dx);
                surfaceNodesID.push_back(surfaceNodeID);
                surfaceNodesPosIndex.push_back(fullyDamagedParticlesNodesVec[i].posIndex);
            }
            //************find surface nodes//
        }
    }

    // clean isolate nodes and sharp nodes
    std::vector<int> boundaryNodesID; // store boundary nodes ID
    std::vector<Eigen::Vector3i> boundaryNodesPosIndex; // store boundary position index
    for (int m = 0; m < boundaryNodesPosIndexNoClean.size(); m++) 
    {
        bool stopIsolate = false;
        for (int axis = 0; axis < 3; axis++) // three axis directions
        {
            if(stopIsolate == false)
            {
                for (int pn = 0; pn < 2; pn++) 
                {
                    if(stopIsolate == false)
                    {                    
                        Eigen::Vector3i normal = { 0, 0, 0 };
                        normal[axis] = 2 * pn - 1;
                        Eigen::Vector3i nodePosIndex = boundaryNodesPosIndexNoClean[m] + normal;
                        int nodeID = calculateID(nodePosIndex[0], nodePosIndex[1], nodePosIndex[2], param.length, param.dx);
                        if (count(boundaryNodesIDNoClean.begin(), boundaryNodesIDNoClean.end(), nodeID) != 0) 
                        {
                            boundaryNodesID.push_back(boundaryNodesIDNoClean[m]);
                            boundaryNodesPosIndex.push_back(boundaryNodesPosIndexNoClean[m]);
                            stopIsolate = true;
                        }

                    }
                }
            }
            
        }
    }



    // define a std::map that can find the index of a point in the point vector
    std::map<int, int> pointIndexFind;
    for (int m = 0; m < boundaryNodesID.size(); m++) {
        pointIndexFind[boundaryNodesID[m]] = m;
    }



    ofstream outfile45("boundaryNodes.txt", ios::trunc);
    for (int m = 0; m < boundaryNodesPosIndex.size(); m++) {
        outfile45 << m << " " << scientific << setprecision(8) << boundaryNodesPosIndex[m][0] * param.dx << " " << boundaryNodesPosIndex[m][1] * param.dx << " " << boundaryNodesPosIndex[m][2] * param.dx << endl;
    }
    outfile45.close();

    cout << "The number of boundary nodes is " << boundaryNodesID.size() << endl;

    // find critical nodes
    std::vector<int> criticalNodeIndex;
    int criticalNodeVolumeLength = 1;
    for (int m = 0; m < surfaceNodesPosIndex.size(); m++) {

        for (int i = -criticalNodeVolumeLength; i < criticalNodeVolumeLength + 1; i++) {
            for (int j = -criticalNodeVolumeLength; j < criticalNodeVolumeLength + 1; j++) {
                for (int k = -criticalNodeVolumeLength; k < criticalNodeVolumeLength + 1; k++) {
                    Eigen::Vector3i increment = { i, j, k };
                    Eigen::Vector3i nodePosIndex = surfaceNodesPosIndex[m] + increment;
                    int ID = calculateID(nodePosIndex[0], nodePosIndex[1], nodePosIndex[2], param.length, param.dx);
                    if (count(boundaryNodesID.begin(), boundaryNodesID.end(), ID) != 0) {
                        if (count(criticalNodeIndex.begin(), criticalNodeIndex.end(), ID) == 0) {
                            criticalNodeIndex.push_back(ID);
                            
                        }
                    }
                }
            }
        }
    }



    cout << "Start voro++" << endl;



    double x_min = 1.0E10, x_max = -1.0E10;
    double y_min = 1.0E10, y_max = -1.0E10;
    double z_min = 1.0E10, z_max = -1.0E10;
    int n_x = 10, n_y = 10, n_z = 10;

    for (int m = 0; m <  (*particlesRaw).size(); m++) 
    {
        x_min = std::min(x_min, (*particlesRaw)[m].pos[0]);
        y_min = std::min(y_min, (*particlesRaw)[m].pos[1]);
        z_min = std::min(z_min, (*particlesRaw)[m].pos[2]);

        x_max = std::max(x_max, (*particlesRaw)[m].pos[0]);
        y_max = std::max(y_max, (*particlesRaw)[m].pos[1]);
        z_max = std::max(z_max, (*particlesRaw)[m].pos[2]);
    }
    x_min = x_min - 4 * param.dx;
    x_max = x_max + 4 * param.dx;
    y_min = y_min - 4 * param.dx;
    y_max = y_max + 4 * param.dx;
    z_min = z_min - 4 * param.dx;
    z_max = z_max + 4 * param.dx;



	pre_container pcon(x_min, x_max, y_min, y_max, z_min, z_max, false, false, false);
	// Import the particles from a file
	pcon.import("boundaryNodes.txt");
	pcon.guess_optimal(n_x, n_y, n_z);
	container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, false, false, false, 8);
	pcon.setup(con);


	std::vector<Point> points;
	for (int i = 0; i < boundaryNodesPosIndex.size(); i++)
	{
		Eigen::Vector3d pos = { 0 , 0 , 0 }; // each point's position
		std::vector<Eigen::Vector3d> verticsCoor; // vertices' coordinate
		std::vector<std::vector<int>> verticesFace; // vertices of each face
		std::vector<Eigen::Vector3d> surfaceNormal; // vertices' coordinate
		std::vector<int> neighbour; // neighbour points that share common faces with this point
		std::vector<int> neighbourCalculated; // neighbour points that have already find shared face
		points.push_back(Point(-999, pos, 0, verticsCoor, 0, verticesFace, surfaceNormal, neighbour));
	}





	c_loop_all cl(con);
	if (cl.start())
	{
        
		double x, y, z;
		int index;
		voronoicell_neighbor c;
		std::vector<int> neighbour, verticesFace;
		std::vector<double> vertices;
		std::vector<double> normals;
		do if (con.compute_cell(c, cl))
		{

			int index = cl.pid();
			points[index].index = index;

			cl.pos(x, y, z);
			points[index].pos = { x, y, z };

			c.vertices(x, y, z, vertices);
			points[index].numVertices = vertices.size() / 3;
			for (int m = 0; m < vertices.size() / 3; m++)
			{
				Eigen::Vector3d vert = { vertices[m * 3],  vertices[m * 3 + 1],  vertices[m * 3 + 2] };
				points[index].verticsCoor.push_back(vert);
			}

			c.neighbors(neighbour);
			points[index].neighbour = neighbour;
			points[index].numFaces = neighbour.size();



			c.face_vertices(verticesFace);
			int start = 1, end = verticesFace[0] + 1;
			do
			{
				std::vector<int> faceVert;
				for (int m = start; m < end; m++)
				{
					faceVert.push_back(verticesFace[m]);
				}
				points[index].verticesFace.push_back(faceVert);
				start = end + 1;
				end += verticesFace[end] + 1;
			} while (points[index].verticesFace.size() != neighbour.size());

			c.normals(normals);
			for (int m = 0; m < normals.size() / 3; m++)
			{
				Eigen::Vector3d normal = { normals[m * 3],  normals[m * 3 + 1],  normals[m * 3 + 2] };
				points[index].surfaceNormal.push_back(normal);
			}

			

		} while (cl.inc());
	}




	cout << "Voro++ finished" << endl;







    //***********Read all particles and calculate the damage gradient*************//
    std::vector<Particle> particles;
    std::map<int, int> gridMap;
    std::vector<Grid> nodesVec;
    readParticlesAndCalGradient(&fullyDamagedParticlesNodesVec, &particles, param, &gridMap, &nodesVec);

    double pi = 3.141592653;
    double radius = param.dx * sqrt(3); // support radius of two points
    //radius = 0;

   



	vector<Vector3d> verticesTmp; // vertex index
	vector<vector<int>> facesTmp;



    cout << "Start extracting interior faces" << endl;


    std::vector<set<int>> sharedFace(points.size()); // store the point index which shares a crack surface with the other point

    // find faces that are in the interior and store neighbour information
    for (int i = 0; i < points.size(); i++) {
        Eigen::Vector3d pos = points[i].pos;
        if (ifCriticalNode(pos, param, &criticalNodeIndex) == false) {
            for (int k = 0; k < points[i].neighbour.size(); k++) {
                int neighbourIndex = points[i].neighbour[k];
                if (neighbourIndex > 0) // remove bounding box faces
                {

                    Eigen::Vector3d posNeig = points[neighbourIndex].pos;
                    Eigen::Vector3d posDiff = pos - posNeig;
                    double distancePair = posDiff.norm();

                    if (ifCriticalNode(posNeig, param, &criticalNodeIndex) == false) {
                        if (count(points[i].neighbourCalculated.begin(), points[i].neighbourCalculated.end(), neighbourIndex) == 0) // if these two neighbours have not yet been compared
                        {

                            if (distancePair > radius) // if their distance is larger than the threshold
                            {

                                double existInCrack = ifFullyDamaged((pos + posNeig) / 2.0, param, &gridMap, &nodesVec);
                                if (existInCrack >= param.damageThreshold) // if the middle point is located in the crack area
                                {


                                    vector<int> faceVerteice;
									for (int ver = 0; ver < points[i].verticesFace[k].size(); ver++)
									{
										int indexVertex = points[i].verticesFace[k][ver];
										Vector3d vertex = points[i].verticsCoor[indexVertex];
										int vertexIndex = findIndexVertex(vertex, &verticesTmp);
										faceVerteice.push_back(vertexIndex);
									}
									facesTmp.push_back(faceVerteice);


                                    sharedFace[i].insert(neighbourIndex);
									sharedFace[neighbourIndex].insert(i);

                                    points[i].neighbourOtherSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourOtherSide.push_back(i);

                                } else {
                                    points[i].neighbourSameSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourSameSide.push_back(i);
                                }

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);

                            } else {
                                points[i].neighbourSameSide.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourSameSide.push_back(i);

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);
                            }
                        }
                    }
                }
            }
        }
    }

    cout << "Start extracting boundary faces" << endl;
    // find faces that are in the interior and store neighbour information
    for (int i = 0; i < points.size(); i++) {
        Eigen::Vector3d pos = points[i].pos;

        if (ifCriticalNode(pos, param, &criticalNodeIndex) == false) {
            for (int k = 0; k < points[i].neighbour.size(); k++) {
                int neighbourIndex = points[i].neighbour[k];
                if (neighbourIndex > 0) // remove bounding box faces
                {

                    Eigen::Vector3d posNeig = points[neighbourIndex].pos;
                    Eigen::Vector3d posDiff = pos - posNeig;
                    double distancePair = posDiff.norm();

                    if (ifCriticalNode(posNeig, param, &criticalNodeIndex) == true) {
                        if (count(points[i].neighbourCalculated.begin(), points[i].neighbourCalculated.end(), neighbourIndex) == 0) // if these two neighbours have not yet been compared
                        {

                            if (distancePair > radius) // if their distance is larger than the threshold
                            {
                                //cout << "i = " << i << "; neighbourIndex = " << neighbourIndex << endl;
                                bool twoSide1 = ifTwoSides(i, neighbourIndex, &points, &boundaryNodesID, &boundaryNodesPosIndex, &pointIndexFind, param, &criticalNodeIndex);

                                //twoSide1 = true;
                                if (twoSide1 == true)
                                {



                                    vector<int> faceVerteice;
									for (int ver = 0; ver < points[i].verticesFace[k].size(); ver++)
									{
										int indexVertex = points[i].verticesFace[k][ver];
										Vector3d vertex = points[i].verticsCoor[indexVertex];
										int vertexIndex = findIndexVertex(vertex, &verticesTmp);
										faceVerteice.push_back(vertexIndex);
									}
									facesTmp.push_back(faceVerteice);

                                    sharedFace[i].insert(neighbourIndex);
									sharedFace[neighbourIndex].insert(i);

                                    points[i].neighbourOtherSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourOtherSide.push_back(i);
                                } else {
                                    points[i].neighbourSameSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourSameSide.push_back(i);
                                }

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);

                            } else {
                                points[i].neighbourSameSide.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourSameSide.push_back(i);

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);
                            }
                        }
                    }
                }
            }
        }
    }

    // find faces that defined by critical point
    for (int i = 0; i < points.size(); i++) {

        Eigen::Vector3d pos = points[i].pos;
        if (ifCriticalNode(pos, param, &criticalNodeIndex) == true) {

            for (int k = 0; k < points[i].neighbour.size(); k++) {
                int neighbourIndex = points[i].neighbour[k];
                if (neighbourIndex > 0) // remove bounding box faces
                {
                    Eigen::Vector3d posNeig = points[neighbourIndex].pos;
                    Eigen::Vector3d posDiff = pos - posNeig;
                    double distancePair = posDiff.norm();

                    if (count(points[i].neighbourCalculated.begin(), points[i].neighbourCalculated.end(), neighbourIndex) == 0) // if these two neighbours have not yet been compared
                    {

                        if (distancePair > radius) // if their distance is larger than the threshold
                        {

                            bool twoSide1 = ifTwoSides(i, neighbourIndex, &points, &boundaryNodesID, &boundaryNodesPosIndex, &pointIndexFind, param, &criticalNodeIndex);

                            if (twoSide1 == true)
                            {

                                vector<int> faceVerteice;
                                for (int ver = 0; ver < points[i].verticesFace[k].size(); ver++)
                                {
                                    int indexVertex = points[i].verticesFace[k][ver];
                                    Vector3d vertex = points[i].verticsCoor[indexVertex];
                                    int vertexIndex = findIndexVertex(vertex, &verticesTmp);
                                    faceVerteice.push_back(vertexIndex);
                                }
                                facesTmp.push_back(faceVerteice);



                                sharedFace[i].insert(neighbourIndex);
                                sharedFace[neighbourIndex].insert(i);
                            }

                            // store computation information.
                            points[i].neighbourCalculated.push_back(neighbourIndex);
                            points[neighbourIndex].neighbourCalculated.push_back(i);
                        }
                    }
                }
            }
        }
    }



    endTime = clock();
    cout << "Totle Time is : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;



	////////////////////////////////
	// calculate all fragments
	////////////////////////////////
	std::set<int> remainingPoints;
	for (int i = 0; i < points.size(); i++)
	{
		remainingPoints.insert(i);
	}

	std::vector<std::vector<int>> allFragments;

	do
	{

		std::vector<int> fragment;
		int firstElement = *remainingPoints.begin();
		fragment.push_back(firstElement);
		remainingPoints.erase(firstElement);
		int start = 0; // the index of starting pointof a single layer
		int countLayer = 1;
		do
		{
			int countLayerUpdate = 0;
			for (int i = start; i < start + countLayer; i++) // parse each candidate
			{
				int currentPoint = fragment[i];
				for (int j = 0; j < points[currentPoint].neighbour.size(); j++)
				{
					int candPoint = points[currentPoint].neighbour[j];
					if (candPoint >= 0)
					{
						if (sharedFace[currentPoint].find(candPoint) == sharedFace[currentPoint].end())
						{
							if (std::find(fragment.begin(), fragment.end(), candPoint) == fragment.end())
							{
								fragment.push_back(candPoint);
								remainingPoints.erase(candPoint);
								countLayerUpdate += 1;
							}

						}
					}

				}
			}
			start = start + countLayer;
			countLayer = countLayerUpdate;

		} while (countLayer != 0);

		allFragments.push_back(fragment);



	} while (remainingPoints.size() != 0);

	std::cout << "Number of raw fragments is =" << allFragments.size() << endl;




	// add interior volume to the surrounding fragment
	std::vector<std::vector<int>> allFragmentsRemoveInterior; // store each fragment where interior volume is added to surrounding volume
	std::set<int> interiorFragmentIndex; // store the index of fragment which is an interior volume
	for (int i = 0; i < allFragments.size(); i++)
	{
		std::vector<int> fragment = allFragments[i];
		bool interior = true;
		for (int j = 0; j < fragment.size(); j++)
		{
			for (int h = 0; h < points[fragment[j]].neighbour.size(); h++)
			{
				int neigPoint = points[fragment[j]].neighbour[h];
				if (neigPoint < 0)
				{
					interior = false;
					break;
				}
			}

			if (interior == false)
			{
				break;
			}
		}

		if (interior == true)
		{
			interiorFragmentIndex.insert(i);
		}


	}



	if (interiorFragmentIndex.size() == 0) // all fargments are not interior
	{
		allFragmentsRemoveInterior = allFragments;
	}
	else
	{
		std::map<int, int> mergeInterior; // key and value are the index of interior volume and surrounding volume respectively

		// find the surrounding volume that an interior volume should be added into
		std::set<int>::iterator it;
		for (it = interiorFragmentIndex.begin(); it != interiorFragmentIndex.end(); ++it)
		{
			// find neighbour point in the surrounding volume
			int surroundPoint = -99; // surrounding neighbour point of this interior volume
			for (int k = 0; k < allFragments[*it].size(); k++)
			{
				int seedPoint = allFragments[*it][k];
				for (int i = 0; i < points[seedPoint].neighbour.size(); i++)
				{
					int neigPoint = points[seedPoint].neighbour[i];
					if (sharedFace[seedPoint].find(neigPoint) != sharedFace[seedPoint].end())
					{
						surroundPoint = neigPoint;
						break;
					}
				}

				if (surroundPoint >= 0)
				{
					break;
				}
			}


			// find surrounding volume
			int surroundingVolume = -99;
			for (int i = 0; i < allFragments.size(); i++)
			{
				if (std::find(allFragments[i].begin(), allFragments[i].end(), surroundPoint) != allFragments[i].end())
				{
					surroundingVolume = i;
					break;
				}
			}
			mergeInterior[*it] = surroundingVolume;
            if(surroundingVolume < 0)
            {
                bool findCrackSurface = false;
                meshObjFormat nocrackmesh;
                std::vector<meshObjFormat> noFragments;
                std::tuple<bool, meshObjFormat, meshObjFormat, std::vector<meshObjFormat> > resultReturn(findCrackSurface, nocrackmesh, nocrackmesh, noFragments);
                return resultReturn;
            }
            ASSERT(surroundingVolume >= 0);

		}


		// add interior volume into the surrounding volume
		std::map<int, int>::iterator itMap;
		for (itMap = mergeInterior.begin(); itMap != mergeInterior.end(); itMap++)
		{
			int interiorFragment = itMap->first;
			int surroundingFragment = itMap->second;
			allFragments[surroundingFragment].insert(allFragments[surroundingFragment].end(), allFragments[interiorFragment].begin(), allFragments[interiorFragment].end());
		}


		for (int i = 0; i < allFragments.size(); i++)
		{
			if (mergeInterior.find(i) == mergeInterior.end())
			{
				allFragmentsRemoveInterior.push_back(allFragments[i]);
			}
		}


	}


	std::cout << "Number of final fragments is =" << allFragmentsRemoveInterior.size() << endl;




		
	// remove duplicated vertices
	std::vector<meshObjFormat> allFragmentsObj;
	for (int i = 0; i < allFragmentsRemoveInterior.size(); i++)
	{
		std::vector<int> fragment = allFragmentsRemoveInterior[i];
	
		vector<Vector3d> verticesEachFrag;
		vector<vector<int>> facesEachFrag;
		for (int k = 0; k < fragment.size(); k++)
		{
			// find faces of each voronoi cell
			std::set<int> voroCellVertices;
			std::vector<std::vector<int>> voroCellFaces;
			for (int h = 0; h < points[fragment[k]].verticesFace.size(); h++)
			{
				int opponentPoint = points[fragment[k]].neighbour[h];
				if (std::find(fragment.begin(), fragment.end(), opponentPoint) == fragment.end())
				{
					voroCellFaces.push_back(points[fragment[k]].verticesFace[h]);
					for (int f = 0; f < points[fragment[k]].verticesFace[h].size(); f++)
					{
						voroCellVertices.insert(points[fragment[k]].verticesFace[h][f]);
					}
				}
			}
	
			// remove duplicated vertices
			std::map<int, int> verticesMapping;
			std::set<int>::iterator it;
			for (it = voroCellVertices.begin(); it != voroCellVertices.end(); ++it)
			{
				Eigen::Vector3d candiVert = points[fragment[k]].verticsCoor[*it];
	
				if (verticesEachFrag.size() == 0)
				{
					verticesMapping[*it] = 0;
					verticesEachFrag.push_back(candiVert);
				}
				else 
				{
					// remove duplicated vertices
					int veryIndex = -1;
					bool findDup = false;
					do
					{
						veryIndex += 1;
						Eigen::Vector3d existVert = verticesEachFrag[veryIndex];
						if (candiVert[0] == existVert[0] && candiVert[1] == existVert[1] && candiVert[2] == existVert[2])
						{
							findDup = true;
						}
					} while (findDup == false && veryIndex != verticesEachFrag.size() - 1);
	
	
					if (findDup == true)
					{
						verticesMapping[*it] = veryIndex;
					}
					else
					{
						verticesMapping[*it] = verticesEachFrag.size();
						verticesEachFrag.push_back(candiVert);
					}
	
				}
	
	
			}
	
	
			// reorder face vertices
			for (int f = 0; f < voroCellFaces.size(); f++)
			{
				std::vector<int> orderedFace;
				for (int d = 0; d < voroCellFaces[f].size(); d++)
				{
					orderedFace.push_back(verticesMapping[voroCellFaces[f][d]]);
				}
				facesEachFrag.push_back(orderedFace);
					
			}
	
		}

			
		meshObjFormat fragmentObj;
		fragmentObj.faces = facesEachFrag;
		fragmentObj.vertices = verticesEachFrag;
		allFragmentsObj.push_back(fragmentObj);
		
	}
	
	
	
    // find crack surface that does fully cut
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<int>> faces;
    std::set<std::string> faceSampled;
	for (int i = 0; i < allFragmentsRemoveInterior.size(); i++)
	{
		std::vector<int> fragment = allFragmentsRemoveInterior[i];
	
		for (int k = 0; k < fragment.size(); k++)
		{
			// find faces of each voronoi cell
			for (int h = 0; h < points[fragment[k]].verticesFace.size(); h++)
			{
				int opponentPoint = points[fragment[k]].neighbour[h];
                if(opponentPoint >= 0 && std::find(fragment.begin(), fragment.end(), opponentPoint) == fragment.end())
                {
                    std::string facePositive = std::to_string(fragment[k]) + "#" + std::to_string(opponentPoint);
                    std::string faceNegative = std::to_string(opponentPoint) + "#" + std::to_string(fragment[k]);
                    if(faceSampled.find(facePositive) == faceSampled.end() || faceSampled.find(faceNegative) == faceSampled.end() )
                    {
                        faceSampled.insert(facePositive);
                        faceSampled.insert(faceNegative);

                        int numOfVert = (int)vertices.size();
                        std::vector<int> face;
                        int co = 0;
                        for (int f = 0; f < points[fragment[k]].verticesFace[h].size(); f++)
                        {
                            int vertIndex = points[fragment[k]].verticesFace[h][f];
                            Eigen::Vector3d vertPos = points[fragment[k]].verticsCoor[vertIndex];

                            int vertIndexVertices = -999;
                            for(int fg = 0; fg < vertices.size(); fg++)
                            {
                                Eigen::Vector3d existVert = vertices[fg];
                                if (vertPos[0] == existVert[0] && vertPos[1] == existVert[1] && vertPos[2] == existVert[2])
                                {
                                    vertIndexVertices = fg;
                                    break;
                                }
                            }

                            if(vertIndexVertices == -999)
                            {
                                vertices.push_back(vertPos);
                                face.push_back(numOfVert + co);
                                co += 1;
                            }
                            else
                            {
                                face.push_back(vertIndexVertices);
                            }
                            
                        }
                        faces.push_back(face);


                    }

                }

			}
	

		}


	
	}



    bool findCrackSurface = true;
    meshObjFormat crackSurfacePartialCut;
    crackSurfacePartialCut.vertices = verticesTmp;
    crackSurfacePartialCut.faces = facesTmp;  

    meshObjFormat crackSurfaceFullCut;
    crackSurfaceFullCut.vertices = vertices;
    crackSurfaceFullCut.faces = faces;  

    std::cout << "The number of crack faces is " << facesTmp.size() << endl;

    std::tuple<bool, meshObjFormat,  meshObjFormat, std::vector<meshObjFormat> > resultReturn(findCrackSurface, crackSurfacePartialCut, crackSurfaceFullCut, allFragmentsObj);
    return resultReturn; 


}

// clean faces from voro output
tuple<std::vector<Eigen::Vector3d>, std::vector<vector<int>>, std::vector<vector<int>>> cleanMeshAfterVoro(std::vector<Eigen::Vector3d>* vertices, std::vector<std::vector<int>>* faces)
{
    std::vector<Eigen::Vector3d> verticesVec;
    std::vector<std::vector<int>> facesVec;

    for (int k = 0; k < (*vertices).size(); k++) {
        verticesVec.push_back((*vertices)[k]);
    }
    for (int k = 0; k < (*faces).size(); k++) {
        facesVec.push_back((*faces)[k]);
    }

    //**************** build duplicated vertices index ****************//
    vector<int> vertexDuplicate(verticesVec.size()); // store the index where a vertex appears in the first time. vertexDuplicate[4]=0 means the fifth vertex appears in the first place of vector vertices
    for (int i = 0; i < (*vertices).size(); i++) {
        vertexDuplicate[i] = i;
    }

    for (int i = 0; i < (*vertices).size(); i++) {
        Eigen::Vector3d pos1 = (*vertices)[i];

        for (int j = i + 1; j < (*vertices).size(); j++) {

            Eigen::Vector3d pos2 = (*vertices)[j];
            if (pos1[0] == pos2[0] && pos1[1] == pos2[1] && pos1[2] == pos2[2]) {
                int index = vertexDuplicate[i];
                vertexDuplicate[j] = index;
            }
        }
    }
    //**************** build duplicated vertices index ****************//

    //**************** remove duplicated vertices ****************//
    vector<Eigen::Vector3d> verticesClean;
    int count = 0;
    for (int i = 0; i < vertexDuplicate.size(); i++) {
        if (vertexDuplicate[i] == i) {
            verticesClean.push_back((*vertices)[i]);
            vertexDuplicate[i] = count;
            count += 1;
        } else {
            int index = vertexDuplicate[i];
            vertexDuplicate[i] = vertexDuplicate[index];
        }
    }
    //**************** remove duplicated vertices ****************//

    //**************** reorder face (*vertices) ****************//
    for (int i = 0; i < (*faces).size(); i++) {
        for (int j = 0; j < (*faces)[i].size(); j++) {
            int index = (*faces)[i][j];
            facesVec[i][j] = vertexDuplicate[index];
        }
    }
    //**************** reorder face (*vertices) ****************//

    vector<vector<int>> facesClean;
    for (int k = 0; k < facesVec.size(); k++) {
        set<int> faceVert;
        for (int d = 0; d < facesVec[k].size(); d++) {
            faceVert.insert(facesVec[k][d]);
        }

        if (faceVert.size() >= 3) {
            vector<int> oneFaceVertices;
            for (int d = 0; d < facesVec[k].size(); d++) {
                if (faceVert.count(facesVec[k][d]) == 1) {
                    oneFaceVertices.push_back(facesVec[k][d]);
                    faceVert.erase(facesVec[k][d]);
                }
            }
            facesClean.push_back(oneFaceVertices);
        }
    }

    ////**************** remove duplicate vertices in a face ****************//
    //int maxNumberOfVertices = 0; // the maximum number of (*vertices) in all (*faces)
    //for (int k = 0; k < facesVec.size(); k++)
    //{
    //
    //	set<int> temp;
    //	for (int g = 0; g < facesVec[k].size(); g++)
    //	{
    //		temp.insert(facesVec[k][g]);
    //	}

    //	vector<int> vert;
    //	set<int>::iterator it;
    //	for (it = temp.begin(); it != temp.end(); it++)
    //	{
    //		vert.push_back(*it);
    //	}

    //	facesVec[k] = vert;

    //	if (vert.size() > maxNumberOfVertices)
    //	{
    //		maxNumberOfVertices = vert.size();
    //	}
    //}
    ////**************** remove duplicate vertices in a face ****************//

    //// group (*faces) according to its number of (*vertices)
    //vector<int> minVertextIndex(facesVec.size()); // the minimum vertex index of each face
    //vector<vector<int>> facesGroup(maxNumberOfVertices + 1);
    //for (int k = 0; k < facesVec.size(); k++)
    //{
    //	int numVert = facesVec[k].size();
    //	facesGroup[numVert].push_back(k);

    //	int minIndex = *std::min_element(facesVec[k].begin(), facesVec[k].end());
    //	minVertextIndex[k] = minIndex;
    //}

    //set<int> duplicateFace;
    //for (int k = 0; k < facesGroup.size(); k++)
    //{
    //	int groupSize = facesGroup[k].size();
    //	if (groupSize == 1)
    //	{
    //		duplicateFace.insert(facesGroup[k][0]);
    //	}
    //	else if (groupSize > 1)
    //	{
    //		// compare this face with the remaining (*faces)
    //		for (int m = 0; m < groupSize; m++)
    //		{
    //			int currentFaceIndex = facesGroup[k][m];
    //			set<int> sm(facesVec[currentFaceIndex].begin(), facesVec[currentFaceIndex].end());
    //			for (int n = m + 1; n < groupSize; n++)
    //			{
    //				int nextFaceIndex = facesGroup[k][n];
    //				if (minVertextIndex[currentFaceIndex] == minVertextIndex[nextFaceIndex])
    //				{
    //					set<int> sn(facesVec[nextFaceIndex].begin(), facesVec[nextFaceIndex].end());
    //					set<int> smn;
    //					set_intersection(sm.begin(), sm.end(), sn.begin(), sn.end(), inserter(smn, smn.begin()));

    //					if (smn.size() == sm.size())
    //					{
    //						duplicateFace.insert(m);
    //					}

    //				}
    //			}

    //		}

    //	}

    //}

    vector<vector<int>> verticesInfor(verticesClean.size()); // vector<vector<int>>: store the list of faces which has this vertex
    for (int k = 0; k < facesClean.size(); k++) {
        for (int j = 0; j < facesClean[k].size(); j++) {
            verticesInfor[facesClean[k][j]].push_back(k);
        }
    }

    tuple<vector<Eigen::Vector3d>, vector<vector<int>>, vector<vector<int>>> result(verticesClean, facesClean, verticesInfor);

    return result;
}

// insert points between two points
std::vector<Eigen::Vector3d> insertPoints(Eigen::Vector3d startPos, Eigen::Vector3d endPos, double dx)
{
    std::vector<Eigen::Vector3d> points;
    Eigen::Vector3d direction = (endPos - startPos).normalized();
    double length = (endPos - startPos).norm();

    int i = 1;
    if ((dx - length) < -1.0E-10) {
        do {
            Eigen::Vector3d pos = startPos + direction * (double)i * dx;
            points.push_back(pos);
            i += 1;
        } while (((double)i * dx - length) < -1.0E-10); // To prevent the last inserted point overlap with the ending point
    }

    return points;
}

// refine mesh
struct meshFaceRefined meshRefine(std::vector<Eigen::Vector3d>* vertices, std::vector<std::vector<int>>* faces, vector<vector<int>>* verticesInfor, double dx)
{

    // initialize meshFace structurt vector for adjust normal(order of vertices)
    std::set<int> branchingEdgeVertices; // vertices that are in the branching region
    std::vector<meshFace> meshTemp((*faces).size());
    for (int j = 0; j < (*faces).size(); j++) {
        meshTemp[j].index = j; // initialize index

        for (int i = 0; i < (*faces)[j].size(); i++) {
            // initialize vertices information
            meshTemp[j].vertices.push_back((*faces)[j][i]);

            int startVertex = (*faces)[j][i];
            int endVertex = -100; // -100 is a negative value
            if (i == (*faces)[j].size() - 1) {
                endVertex = (*faces)[j][0];
            } else {
                endVertex = (*faces)[j][i + 1];
            }

            // initialize edge information
            string edgeName = std::to_string(startVertex) + "*" + std::to_string(endVertex);
            pair<int, int> edgeVertex(startVertex, endVertex);
            meshTemp[j].edge[edgeName] = edgeVertex;

            // initialize neighbour information
            set<int> faceIntersect;
            std::set_intersection((*verticesInfor)[startVertex].begin(), (*verticesInfor)[startVertex].end(), (*verticesInfor)[endVertex].begin(), (*verticesInfor)[endVertex].end(), std::inserter(faceIntersect, faceIntersect.begin()));
            set<int>::iterator it;
            std::vector<int> neighbour;
            for (it = faceIntersect.begin(); it != faceIntersect.end(); it++) {
                if (*it != j) {
                    neighbour.push_back(*it);
                }
            }
            meshTemp[j].neighbours[edgeName] = neighbour;

            if (faceIntersect.size() >= 3) {
                branchingEdgeVertices.insert(startVertex);
                branchingEdgeVertices.insert(endVertex);
            }
        }
    }

    // adjust the surface normal into the same orientation
    set<int> notYetVisited; // (*faces) that have not been visited yet
    for (int h = 0; h < (*faces).size(); h++) {
        notYetVisited.insert(h);
    }
    vector<int> facesNeedReverse; // (*faces) that need to reverse its order of vertices

    do {
        set<int> facesSet;
        vector<int> facesVector;
        int startFace = *notYetVisited.begin(); // start from face 0
        facesSet.insert(startFace);
        facesVector.push_back(startFace);
        notYetVisited.erase(startFace);

        int layerStart = 0;
        int currentLayerLength = 1;
        int nextLayerLength = 0;
        bool finished = false; // judge if all (*faces) in this layer have been visited

        do {
            //cout << "facesVector.size() = " << facesVector.size()<<";  layerStart = "<< layerStart<< ";  layerEnd = " << layerStart + currentLayerLength << endl;
            nextLayerLength = 0;
            finished = true; // judge if all (*faces) connected to the starting face have been visited
            for (int k = layerStart; k < layerStart + currentLayerLength; k++) {
                int currentFaceIndex = facesVector[k];
                map<string, pair<int, int>>::iterator it;

                for (it = meshTemp[currentFaceIndex].edge.begin(); it != meshTemp[currentFaceIndex].edge.end(); it++) {

                    string edgeNameReverse = std::to_string(it->second.second) + "*" + std::to_string(it->second.first);
                    vector<int> neighFaces = meshTemp[currentFaceIndex].neighbours[it->first];

                    if (branchingEdgeVertices.count(it->second.second) == 0 && branchingEdgeVertices.count(it->second.first) == 0) {

                        for (int m = 0; m < neighFaces.size(); m++) {
                            int neighIndex = neighFaces[m];

                            if (facesSet.count(neighIndex) == 0) {
                                // ambuigity in the branching area
                                for (int n = 0; n < meshTemp[neighIndex].vertices.size(); n++) {

                                    if (meshTemp[neighIndex].edge.find(edgeNameReverse) == meshTemp[neighIndex].edge.end()) {

                                        meshTemp[neighIndex].edge.clear();
                                        vector<int> reorder = meshTemp[neighIndex].vertices;
                                        reverse(std::begin(reorder), std::end(reorder));
                                        for (int i = 0; i < reorder.size(); i++) {
                                            int startVertex = reorder[i];
                                            int endVertex = -100; // -100 is a negative value
                                            if (i == reorder.size() - 1) {
                                                endVertex = reorder[0];
                                            } else {
                                                endVertex = reorder[i + 1];
                                            }

                                            // initialize edge information
                                            string edgeName = std::to_string(startVertex) + "*" + std::to_string(endVertex);
                                            string edgeNameReorder = std::to_string(endVertex) + "*" + std::to_string(startVertex);
                                            pair<int, int> edgeVertex(startVertex, endVertex);
                                            meshTemp[neighIndex].edge[edgeName] = edgeVertex;
                                            meshTemp[neighIndex].neighbours[edgeName] = meshTemp[neighIndex].neighbours[edgeNameReorder];
                                        }

                                        facesNeedReverse.push_back(neighIndex);
                                    }

                                    facesSet.insert(neighIndex);
                                    facesVector.push_back(neighIndex);
                                    notYetVisited.erase(neighIndex);
                                    nextLayerLength += 1;
                                }
                            }
                        }
                    }
                }
            }

            layerStart = layerStart + currentLayerLength;
            currentLayerLength = nextLayerLength;

        } while (nextLayerLength != 0);

    } while (notYetVisited.size() != 0);

    for (int h = 0; h < facesNeedReverse.size(); h++) {
        reverse(std::begin((*faces)[facesNeedReverse[h]]), std::end((*faces)[facesNeedReverse[h]]));
    }

    // initialize meshFace structurt vector
    std::vector<meshFace> mesh((*faces).size());
    for (int j = 0; j < (*faces).size(); j++) {
        mesh[j].index = j; // initialize index
        mesh[j].insertPointsAfterVertex.resize((*faces)[j].size()); // initialize insertPointsAfterVertex vector

        for (int i = 0; i < (*faces)[j].size(); i++) {
            // initialize vertices information
            mesh[j].vertices.push_back((*faces)[j][i]);

            int startVertex = (*faces)[j][i];
            int endVertex = -100; // -100 is a negative value
            if (i == (*faces)[j].size() - 1) {
                endVertex = (*faces)[j][0];
            } else {
                endVertex = (*faces)[j][i + 1];
            }

            // initialize vertex position vector
            mesh[j].vertexPosition[startVertex] = i;

            // initialize edge information
            string edgeName = std::to_string(startVertex) + "*" + std::to_string(endVertex);
            pair<int, int> edgeVertex(startVertex, endVertex);
            mesh[j].edge[edgeName] = edgeVertex;
            double edgeLength = ((*vertices)[startVertex] - (*vertices)[endVertex]).norm();
            mesh[j].edgeLength[edgeName] = edgeLength;
            if (edgeLength <= dx) {
                mesh[j].edgeInsertedOrNot[edgeName] = true;
            } else {
                mesh[j].edgeInsertedOrNot[edgeName] = false;
            }

            // initialize neighbour information
            set<int> faceIntersect;
            std::set_intersection((*verticesInfor)[startVertex].begin(), (*verticesInfor)[startVertex].end(), (*verticesInfor)[endVertex].begin(), (*verticesInfor)[endVertex].end(), std::inserter(faceIntersect, faceIntersect.begin()));
            set<int>::iterator it;
            std::vector<int> neighbour;
            for (it = faceIntersect.begin(); it != faceIntersect.end(); it++) //使用迭代器进行遍历
            {
                if (*it != j) {
                    neighbour.push_back(*it);
                }
            }
            mesh[j].neighbours[edgeName] = neighbour;

            if (faceIntersect.size() >= 3) {
                branchingEdgeVertices.insert(startVertex);
                branchingEdgeVertices.insert(endVertex);
            }
        }

        // calculate face normal and unit vectors
        vector<Eigen::Vector3d> normalVec;
        double longestEdgeLength = -100;
        Eigen::Vector3d longestEdge = { 0, 0, 0 };
        for (int kk = 0; kk < mesh[j].vertices.size(); kk++) {

            int start = mesh[j].vertices[kk];
            int middle = 0;
            int end = 0;
            if (kk == mesh[j].vertices.size() - 2) {
                middle = mesh[j].vertices[kk + 1];
                end = mesh[j].vertices[0];
            } else if (kk == mesh[j].vertices.size() - 1) {
                middle = mesh[j].vertices[0];
                end = mesh[j].vertices[1];
            } else {
                middle = mesh[j].vertices[kk + 1];
                end = mesh[j].vertices[kk + 2];
            }

            Eigen::Vector3d v1 = (*vertices)[middle] - (*vertices)[start];
            Eigen::Vector3d v2 = (*vertices)[end] - (*vertices)[middle];

            if (v1.norm() > longestEdgeLength) {
                longestEdgeLength = v1.norm();
                longestEdge = v1;
            }

            if (v2.norm() > longestEdgeLength) {
                longestEdgeLength = v2.norm();
                longestEdge = v2;
            }

            Eigen::Vector3d v3 = v1.cross(v2).normalized();
            normalVec.push_back(v3);
        }

        Eigen::Vector3d normal = { 0, 0, 0 };
        for (int hh = 0; hh < normalVec.size(); hh++) {
            normal += normalVec[hh];
        }

        mesh[j].normal = normal / (double)normalVec.size();
        mesh[j].u = longestEdge.normalized();
        mesh[j].v = (longestEdge.cross(normal)).normalized();
    }

    // insert points between vertices in an edge
    int countNumOfVertices = (*vertices).size(); // count the number of vertices in the vertices vector
    for (int j = 0; j < mesh.size(); j++) {
        std::map<string, bool>::iterator it;
        for (it = mesh[j].edgeInsertedOrNot.begin(); it != mesh[j].edgeInsertedOrNot.end(); it++) {
            string edgeName = it->first;
            if (it->second == false) {
                // insert points into this face
                pair<int, int> twoEndingVertex = mesh[j].edge[edgeName];
                vector<Eigen::Vector3d> points = insertPoints((*vertices)[twoEndingVertex.first], (*vertices)[twoEndingVertex.second], dx);
                vector<int> insertPointsVertex;
                for (int k = 0; k < points.size(); k++) {
                    insertPointsVertex.push_back(countNumOfVertices);
                    if (mesh[j].neighbours[edgeName].size() > 2) {
                        branchingEdgeVertices.insert(countNumOfVertices);
                    }
                    (*vertices).push_back(points[k]);
                    countNumOfVertices += 1;
                }
                int vertexPos = mesh[j].vertexPosition[twoEndingVertex.first];
                mesh[j].insertPointsAfterVertex[vertexPos] = insertPointsVertex;
                mesh[j].edgeInsertedOrNot[edgeName] = true; // this edge has already been interpolated

                // insert points into its neighbours
                vector<int> neighbourFaceIndexVec = mesh[j].neighbours[edgeName];
                for (int m = 0; m < neighbourFaceIndexVec.size(); m++) {
                    int neighourFaceIndex = neighbourFaceIndexVec[m];
                    string positiveSequenceName = edgeName;
                    string negativeSequenceName = std::to_string(twoEndingVertex.second) + "*" + std::to_string(twoEndingVertex.first);

                    if (mesh[neighourFaceIndex].edge.find(positiveSequenceName) != mesh[neighourFaceIndex].edge.end()) {
                        if (mesh[neighourFaceIndex].edgeInsertedOrNot[positiveSequenceName] == false) {
                            int vertexPos = mesh[neighourFaceIndex].vertexPosition[twoEndingVertex.first];
                            mesh[neighourFaceIndex].insertPointsAfterVertex[vertexPos] = insertPointsVertex;
                            mesh[neighourFaceIndex].edgeInsertedOrNot[positiveSequenceName] = true; // this edge has already been interpolated
                        }

                    } else {
                        if (mesh[neighourFaceIndex].edgeInsertedOrNot[negativeSequenceName] == false) {
                            vector<int> insertPointsVertexReverse = insertPointsVertex;
                            std::reverse(insertPointsVertexReverse.begin(), insertPointsVertexReverse.end());
                            int vertexPos = mesh[neighourFaceIndex].vertexPosition[twoEndingVertex.second];
                            mesh[neighourFaceIndex].insertPointsAfterVertex[vertexPos] = insertPointsVertexReverse;
                            mesh[neighourFaceIndex].edgeInsertedOrNot[negativeSequenceName] = true; // this edge has already been interpolated
                        }
                    }
                }
            }
        }
    }

    // update edge by insertd points
    for (int j = 0; j < mesh.size(); j++) {
        for (int k = 0; k < mesh[j].vertices.size(); k++) {
            mesh[j].verticesUpdated.push_back(mesh[j].vertices[k]);
            vector<int> insertPoints = mesh[j].insertPointsAfterVertex[k];
            if (insertPoints.size() != 0) {
                mesh[j].verticesUpdated.insert(mesh[j].verticesUpdated.end(), insertPoints.begin(), insertPoints.end());
            }
        }
    }

    // sample points inside of each face
    for (int j = 0; j < mesh.size(); j++) {

        // vertices' position in 2D plane
        vector<Eigen::Vector2d> verticesIn2D;
        double minX = 1.0E10, maxX = -1.0E10, minY = 1.0E10, maxY = -1.0E10; // record the min and max coordinate of the plane
        for (int k = 0; k < mesh[j].vertices.size(); k++) {
            Eigen::Vector3d pos = (*vertices)[mesh[j].vertices[k]];
            Eigen::Vector2d posIn2D = { pos.dot(mesh[j].u), pos.dot(mesh[j].v) };
            if (posIn2D[0] < minX) {
                minX = posIn2D[0];
            }
            if (posIn2D[0] > maxX) {
                maxX = posIn2D[0];
            }
            if (posIn2D[1] < minY) {
                minY = posIn2D[1];
            }
            if (posIn2D[1] > maxY) {
                maxY = posIn2D[1];
            }

            verticesIn2D.push_back(posIn2D);
        }

        // calculate the updatEdgeVertices' coordinates in 2D
        for (int k = 0; k < mesh[j].verticesUpdated.size(); k++) {
            Eigen::Vector3d pos = (*vertices)[mesh[j].verticesUpdated[k]];
            Eigen::Vector2d posIn2D = { pos.dot(mesh[j].u), pos.dot(mesh[j].v) };
            mesh[j].verticesUpdatedCoor2D.push_back(posIn2D);
        }

        // calculate the edge vector
        vector<Eigen::Vector2d> edgeVector;
        for (int k = 0; k < verticesIn2D.size(); k++) {
            Eigen::Vector2d egvec = { 0, 0 };
            if (k == verticesIn2D.size() - 1) {
                egvec = verticesIn2D[0] - verticesIn2D[verticesIn2D.size() - 1];
            } else {
                egvec = verticesIn2D[k + 1] - verticesIn2D[k];
            }
            edgeVector.push_back(egvec);
        }

        // calculate matrix and vector for projecting
        double NP = mesh[j].normal.dot((*vertices)[mesh[j].vertices[0]]);
        Eigen::Matrix3d coffI;
        coffI << mesh[j].u[0], mesh[j].u[1], mesh[j].u[2],
            mesh[j].v[0], mesh[j].v[1], mesh[j].v[2],
            mesh[j].normal[0], mesh[j].normal[1], mesh[j].normal[2];
        Eigen::Matrix3d coffInv = coffI.inverse();

        // https://blog.csdn.net/theArcticOcean/article/details/48632391
        // sample points
        Eigen::Vector2d base = { minX - 0.5 * dx, minY - 0.5 * dx }; // local origin of sample points
        for (int ix = 0; ix < round((maxX - minX) / dx) + 2; ix++) {
            for (int iy = 0; iy < round((maxY - minY) / dx) + 2; iy++) {
                Eigen::Vector2d diffVec = { ix * dx, iy * dx };
                Eigen::Vector2d candidate = base + diffVec;

                // judge if this point is inside of the polygon or not
                bool inside = true;
                Eigen::Vector2d firstCandidateVec = candidate - verticesIn2D[0];
                double firstCrossProductVec = edgeVector[0][0] * firstCandidateVec[1] - edgeVector[0][1] * firstCandidateVec[0];

                for (int k = 1; k < verticesIn2D.size(); k++) {
                    Eigen::Vector2d candidateVec = candidate - verticesIn2D[k];
                    double crossProductVec = edgeVector[k][0] * candidateVec[1] - edgeVector[k][1] * candidateVec[0];
                    if (firstCrossProductVec * crossProductVec <= 0) {
                        inside = false;
                        break;
                    }
                    firstCrossProductVec = crossProductVec;
                }

                if (inside == true) {
                    mesh[j].sampledPointsCoor2D.push_back(candidate);
                    //project this point to 3D and calculate the coordinate
                    Eigen::Vector3d pos = { candidate[0], candidate[1], NP };
                    (*vertices).push_back(coffInv * pos);
                    mesh[j].sampledPoints.push_back(countNumOfVertices);
                    countNumOfVertices += 1;
                }
            }
        }
    }


    std::cout<<"trianglate polygons"<<std::endl;

    // triangulate each polygon
    for (int j = 0; j < mesh.size(); j++) {
        // prepare input struct
        triangulateio in;
        int numberofpoints = mesh[j].sampledPointsCoor2D.size() + mesh[j].verticesUpdatedCoor2D.size();
        int numberofsegments = mesh[j].verticesUpdatedCoor2D.size();

        in.pointlist = (double*)calloc(numberofpoints * 2, sizeof(double));
        for (int k = 0; k < mesh[j].verticesUpdatedCoor2D.size(); k++) {
            in.pointlist[2 * k] = mesh[j].verticesUpdatedCoor2D[k][0];
            in.pointlist[2 * k + 1] = mesh[j].verticesUpdatedCoor2D[k][1];
        }
        for (int k = 0; k < mesh[j].sampledPointsCoor2D.size(); k++) {
            int updatedIndex = k + mesh[j].verticesUpdatedCoor2D.size();
            in.pointlist[2 * updatedIndex] = mesh[j].sampledPointsCoor2D[k][0];
            in.pointlist[2 * updatedIndex + 1] = mesh[j].sampledPointsCoor2D[k][1];
        }

        in.segmentlist = (int*)calloc(numberofsegments * 2, sizeof(int));
        for (int k = 0; k < mesh[j].verticesUpdatedCoor2D.size(); k++) {
            // triangle library index starts from 1
            in.segmentlist[2 * k] = k + 1;
            if (k == mesh[j].verticesUpdatedCoor2D.size() - 1) {
                in.segmentlist[2 * k + 1] = 1;
            } else {
                in.segmentlist[2 * k + 1] = k + 2;
            }
        }

        in.numberofsegments = numberofsegments;
        in.numberofpoints = numberofpoints;
        in.pointattributelist = NULL;
        in.pointmarkerlist = NULL;
        in.numberofpointattributes = 0;
        in.trianglelist = NULL;
        in.triangleattributelist = NULL;
        in.trianglearealist = NULL;
        in.neighborlist = NULL;
        in.numberoftriangles = 0;
        in.numberofcorners = 0;
        in.numberoftriangleattributes = 0;
        in.holelist = NULL;
        in.numberofholes = 0;
        in.regionlist = NULL;
        in.numberofregions = 0;
        in.segmentmarkerlist = NULL;

        // prepare output struct
        triangulateio out;
        out.pointlist = NULL;
        out.trianglelist = NULL;
        out.segmentlist = NULL;
        out.segmentmarkerlist = NULL;
        out.pointmarkerlist = NULL;

        // prepare flag
        char* flag = "pQ";
        

        // Call triangle
        ::triangulate(flag, &in, &out, 0);

        // read triangulation results
        std::vector<int> mapPoints; // map the triangulation.h points index to a global index
        for (int k = 0; k < mesh[j].verticesUpdated.size(); k++) {
            mapPoints.push_back(mesh[j].verticesUpdated[k]);
        }
        for (int k = 0; k < mesh[j].sampledPoints.size(); k++) {
            mapPoints.push_back(mesh[j].sampledPoints[k]);
        }

        for (int m = 0; m < out.numberoftriangles; m++) {
            Eigen::Vector3i triangle = { mapPoints[out.trianglelist[3 * m] - 1], mapPoints[out.trianglelist[3 * m + 1] - 1], mapPoints[out.trianglelist[3 * m + 2] - 1] };

            set<int> threeVertices;
            threeVertices.insert(triangle[0]);
            threeVertices.insert(triangle[1]);
            threeVertices.insert(triangle[2]);

            // There may be bugs in "triangle.h" library. Three vertices in the same edge(straight line). The following code removes the case
            bool goodTriangle = true;
            for (int n = 0; n < mesh[j].vertices.size(); n++) {
                set<int> verticesInAnEdge;
                verticesInAnEdge.insert(mesh[j].vertices[n]);
                for (int g = 0; g < mesh[j].insertPointsAfterVertex[n].size(); g++) {
                    verticesInAnEdge.insert(mesh[j].insertPointsAfterVertex[n][g]);
                }

                if (n == mesh[j].vertices.size() - 1) {
                    verticesInAnEdge.insert(mesh[j].vertices[0]);
                } else {
                    verticesInAnEdge.insert(mesh[j].vertices[n + 1]);
                }

                set<int> verticesIntersection;
                std::set_intersection(verticesInAnEdge.begin(), verticesInAnEdge.end(), threeVertices.begin(), threeVertices.end(), std::inserter(verticesIntersection, verticesIntersection.begin()));

                if (verticesIntersection.size() == 3) {
                    goodTriangle = false;
                }
            }

            if (goodTriangle == true) {
                if (triangle[0] >= 0 && triangle[0] < (*vertices).size()) {
                    if (triangle[1] >= 0 && triangle[1] < (*vertices).size()) {
                        if (triangle[2] >= 0 && triangle[2] < (*vertices).size()) {
                            mesh[j].trianglesInThisFace.push_back(triangle);
                        }
                    }
                }
            }
        }
    }


     std::cout<<"Finish trianglating polygons"<<std::endl;

    // refined faces
    std::vector<Eigen::Vector3i> facesRefined;
    for (int k = 0; k < mesh.size(); k++) {
        for (int m = 0; m < mesh[k].trianglesInThisFace.size(); m++) {
            Eigen::Vector3i threeVertices = { mesh[k].trianglesInThisFace[m][0], mesh[k].trianglesInThisFace[m][1], mesh[k].trianglesInThisFace[m][2] };
            facesRefined.push_back(threeVertices);
        }
    }

    // refined vertices
    std::vector<Eigen::Vector3d> verticesRefined;
    for (int k = 0; k < (*vertices).size(); k++) {
        verticesRefined.push_back((*vertices)[k]);
    }

    // calculate the average surface normal
    std::vector<Eigen::Vector3d> surfaceNorml((*vertices).size()); // the surface normal for each vertex

    for (int k = 0; k < mesh.size(); k++) {
        for (int m = 0; m < mesh[k].trianglesInThisFace.size(); m++) {
            for (int n = 0; n < 3; n++) {
                int verIndex = mesh[k].trianglesInThisFace[m][n];
                surfaceNorml[verIndex] = mesh[k].normal;
                //if (verIndex >= 0 && verIndex <= (*vertices).size()) //unknown bug, temporary solution
                //{
                //	surfaceNorml[verIndex] = mesh[k].normal;
                //}
            }
        }
    }

    // neighbouring vertices
    vector<set<int>> neighbouringVertex((*vertices).size());
    vector<vector<Eigen::Vector3d>> VertexNormalAverage((*vertices).size());
    for (int k = 0; k < mesh.size(); k++) {
        for (int m = 0; m < mesh[k].trianglesInThisFace.size(); m++) {

            neighbouringVertex[mesh[k].trianglesInThisFace[m][0]].insert(mesh[k].trianglesInThisFace[m][1]);
            neighbouringVertex[mesh[k].trianglesInThisFace[m][0]].insert(mesh[k].trianglesInThisFace[m][2]);

            neighbouringVertex[mesh[k].trianglesInThisFace[m][1]].insert(mesh[k].trianglesInThisFace[m][0]);
            neighbouringVertex[mesh[k].trianglesInThisFace[m][1]].insert(mesh[k].trianglesInThisFace[m][2]);

            neighbouringVertex[mesh[k].trianglesInThisFace[m][2]].insert(mesh[k].trianglesInThisFace[m][0]);
            neighbouringVertex[mesh[k].trianglesInThisFace[m][2]].insert(mesh[k].trianglesInThisFace[m][1]);

            Eigen::Vector3d normal = mesh[k].normal;
            VertexNormalAverage[mesh[k].trianglesInThisFace[m][0]].push_back(normal);
            VertexNormalAverage[mesh[k].trianglesInThisFace[m][1]].push_back(normal);
            VertexNormalAverage[mesh[k].trianglesInThisFace[m][2]].push_back(normal);
        }
    }

    vector<Eigen::Vector3d> VertexNormal;
    for (int k = 0; k < (*vertices).size(); k++) {
        Eigen::Vector3d normal = { 0, 0, 0 };
        for (int m = 0; m < VertexNormalAverage[k].size(); m++) {
            normal += VertexNormalAverage[k][m];
        }
        normal = normal / (double)VertexNormalAverage[k].size();

        VertexNormal.push_back(normal);
    }

    struct meshFaceRefined result;
    result.branchingEdgeVertices = branchingEdgeVertices;
    result.facesRefined = facesRefined;
    result.neighbouringVertex = neighbouringVertex;
    result.VertexNormal = VertexNormal;
    result.verticesRefined = verticesRefined;

    return result;
}





////////////////////////////////
// function for multi-velocity field
////////////////////////////////
crackSurfaceInfor extractCrackSurface_MultiVel(std::vector<Particle>* particlesRaw, std::vector<Eigen::Vector3d>* mpmParticlePositions, struct parametersSim param)
{
    SCOPED_TIMER(__FUNCTION__);
    //*********Read fully damaged particles***********//
    std::vector<Particle> fullyDamagedParticles;
    std::map<int, int> fullyDamagedParticlesGridMap;
    std::vector<Grid> fullyDamagedParticlesNodesVec;

    readParticles(particlesRaw, &fullyDamagedParticles, true, param);
    calDamageGradient(&fullyDamagedParticles, param, param.dx, &fullyDamagedParticlesGridMap, &fullyDamagedParticlesNodesVec);
    setNodeValue(&fullyDamagedParticlesNodesVec, 1);
    findBoundaryNodes(&fullyDamagedParticles, &fullyDamagedParticlesNodesVec, &fullyDamagedParticlesGridMap, param, 2); // This bounding box is used to generate boundary nodes
    findBoundaryNodes(&fullyDamagedParticles, &fullyDamagedParticlesNodesVec, &fullyDamagedParticlesGridMap, param, 3); // This bounding box is used to generate clip surface mesh
    //*********Read fully damaged particles***********//




    //*********Read all particles***********//
    std::vector<Particle> allParticles;
    std::map<int, int> allParticlesGridMap;
    std::vector<Grid> allParticlesNodesVec;

    readParticles(particlesRaw, &allParticles, false, param);
    calDamageGradient(&allParticles, param, param.dx, &allParticlesGridMap, &allParticlesNodesVec);
    setNodeValue(&allParticlesNodesVec, 1);
    findBoundaryNodes(&allParticles, &allParticlesNodesVec, &allParticlesGridMap, param, 2); // This bounding box is used to generate boundary nodes
    findBoundaryNodes(&allParticles, &allParticlesNodesVec, &allParticlesGridMap, param, 3); // This bounding box is used to generate clip surface mesh
    //*********Read all particles***********//


    // find boundary nodes
    std::vector<int> allParticlesNodeBoundaryIndex; // store the index or ID of allParticles boundary nodes
    for (int i = 0; i < allParticlesNodesVec.size(); i++) {
        if (allParticlesNodesVec[i].Di == 2) {
            int ID = calculateID(allParticlesNodesVec[i].posIndex[0], allParticlesNodesVec[i].posIndex[1], allParticlesNodesVec[i].posIndex[2], param.length, param.dx);
            allParticlesNodeBoundaryIndex.push_back(ID);
        }
    }


    //************find surface nodes//
    std::vector<int> surfaceNodesID; // store surface nodes ID
    std::vector<Eigen::Vector3i> surfaceNodesPosIndex; // store surface position index
    std::vector<int> boundaryNodesIDNoClean; // store boundary nodes ID
    std::vector<Eigen::Vector3i> boundaryNodesPosIndexNoClean; // store boundary position index
    for (int i = 0; i < fullyDamagedParticlesNodesVec.size(); i++) {
        if (fullyDamagedParticlesNodesVec[i].Di == 2) {
            int ID = calculateID(fullyDamagedParticlesNodesVec[i].posIndex[0], fullyDamagedParticlesNodesVec[i].posIndex[1], fullyDamagedParticlesNodesVec[i].posIndex[2], param.length, param.dx);
            if (count(allParticlesNodeBoundaryIndex.begin(), allParticlesNodeBoundaryIndex.end(), ID) == 0) {
                boundaryNodesPosIndexNoClean.push_back(fullyDamagedParticlesNodesVec[i].posIndex);
                boundaryNodesIDNoClean.push_back(ID);
            }

            //************find surface nodes//
            if (count(allParticlesNodeBoundaryIndex.begin(), allParticlesNodeBoundaryIndex.end(), ID) != 0) {
                int surfaceNodeID = calculateID(fullyDamagedParticlesNodesVec[i].posIndex[0], fullyDamagedParticlesNodesVec[i].posIndex[1], fullyDamagedParticlesNodesVec[i].posIndex[2], param.length, param.dx);
                surfaceNodesID.push_back(surfaceNodeID);
                surfaceNodesPosIndex.push_back(fullyDamagedParticlesNodesVec[i].posIndex);
            }
            //************find surface nodes//
        }
    }

    // clean isolate nodes and sharp nodes
    std::vector<int> boundaryNodesID; // store boundary nodes ID
    std::vector<Eigen::Vector3i> boundaryNodesPosIndex; // store boundary position index
    for (int m = 0; m < boundaryNodesPosIndexNoClean.size(); m++) {
        for (int axis = 0; axis < 3; axis++) // three axis directions
        {
            for (int pn = 0; pn < 2; pn++) {
                Eigen::Vector3i normal = { 0, 0, 0 };
                normal[axis] = 2 * pn - 1;
                Eigen::Vector3i nodePosIndex = boundaryNodesPosIndexNoClean[m] + normal;
                int nodeID = calculateID(nodePosIndex[0], nodePosIndex[1], nodePosIndex[2], param.length, param.dx);
                if (count(boundaryNodesIDNoClean.begin(), boundaryNodesIDNoClean.end(), nodeID) != 0) {
                    boundaryNodesID.push_back(boundaryNodesIDNoClean[m]);
                    boundaryNodesPosIndex.push_back(boundaryNodesPosIndexNoClean[m]);
                    goto stopNotIsolate;
                }
            }
        }

    stopNotIsolate:
        int a;
    }



    // define a std::map that can find the index of a point in the point vector
    std::map<int, int> pointIndexFind;
    for (int m = 0; m < boundaryNodesID.size(); m++) {
        pointIndexFind[boundaryNodesID[m]] = m;
    }



    ofstream outfile45("boundaryNodes.txt", ios::trunc);
    for (int m = 0; m < boundaryNodesPosIndex.size(); m++) {
        outfile45 << m << " " << scientific << setprecision(8) << boundaryNodesPosIndex[m][0] * param.dx << " " << boundaryNodesPosIndex[m][1] * param.dx << " " << boundaryNodesPosIndex[m][2] * param.dx << endl;
    }
    outfile45.close();



    // find critical nodes
    std::vector<int> criticalNodeIndex;
    int criticalNodeVolumeLength = 1;
    for (int m = 0; m < surfaceNodesPosIndex.size(); m++) {

        for (int i = -criticalNodeVolumeLength; i < criticalNodeVolumeLength + 1; i++) {
            for (int j = -criticalNodeVolumeLength; j < criticalNodeVolumeLength + 1; j++) {
                for (int k = -criticalNodeVolumeLength; k < criticalNodeVolumeLength + 1; k++) {
                    Eigen::Vector3i increment = { i, j, k };
                    Eigen::Vector3i nodePosIndex = surfaceNodesPosIndex[m] + increment;
                    int ID = calculateID(nodePosIndex[0], nodePosIndex[1], nodePosIndex[2], param.length, param.dx);
                    if (count(boundaryNodesID.begin(), boundaryNodesID.end(), ID) != 0) {
                        if (count(criticalNodeIndex.begin(), criticalNodeIndex.end(), ID) == 0) {
                            criticalNodeIndex.push_back(ID);
                            
                        }
                    }
                }
            }
        }
    }




    double x_min = 1.0E10, x_max = -1.0E10;
    double y_min = 1.0E10, y_max = -1.0E10;
    double z_min = 1.0E10, z_max = -1.0E10;
    int n_x = 10, n_y = 10, n_z = 10;

    for (int m = 0; m <  (*mpmParticlePositions).size(); m++) 
    {
        x_min = std::min(x_min, (*mpmParticlePositions)[m][0]);
        y_min = std::min(y_min, (*mpmParticlePositions)[m][1]);
        z_min = std::min(z_min, (*mpmParticlePositions)[m][2]);

        x_max = std::max(x_max, (*mpmParticlePositions)[m][0]);
        y_max = std::max(y_max, (*mpmParticlePositions)[m][1]);
        z_max = std::max(z_max, (*mpmParticlePositions)[m][2]);
    }
    x_min = x_min - 4 * param.dx;
    x_max = x_max + 4 * param.dx;
    y_min = y_min - 4 * param.dx;
    y_max = y_max + 4 * param.dx;
    z_min = z_min - 4 * param.dx;
    z_max = z_max + 4 * param.dx;



	pre_container pcon(x_min, x_max, y_min, y_max, z_min, z_max, false, false, false);
	// Import the particles from a file
	pcon.import("boundaryNodes.txt");
	pcon.guess_optimal(n_x, n_y, n_z);
	container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, false, false, false, 8);
	pcon.setup(con);


	std::vector<Point> points;
	for (int i = 0; i < boundaryNodesPosIndex.size(); i++)
	{
		Eigen::Vector3d pos = { 0 , 0 , 0 }; // each point's position
		std::vector<Eigen::Vector3d> verticsCoor; // vertices' coordinate
		std::vector<std::vector<int>> verticesFace; // vertices of each face
		std::vector<Eigen::Vector3d> surfaceNormal; // vertices' coordinate
		std::vector<int> neighbour; // neighbour points that share common faces with this point
		std::vector<int> neighbourCalculated; // neighbour points that have already find shared face
		points.push_back(Point(-999, pos, 0, verticsCoor, 0, verticesFace, surfaceNormal, neighbour));
	}



    // wallShell wallFace(&fullyDamagedParticlesGridMap, &fullyDamagedParticlesNodesVec, param, -99);
    // con.add_wall(wallFace);
	c_loop_all cl(con);
	if (cl.start())
	{
        
		double x, y, z;
		int index;
		voronoicell_neighbor c;
		std::vector<int> neighbour, verticesFace;
		std::vector<double> vertices;
		std::vector<double> normals;
		do if (con.compute_cell(c, cl))
		{

			int index = cl.pid();
			points[index].index = index;

			cl.pos(x, y, z);
			points[index].pos = { x, y, z };

			c.vertices(x, y, z, vertices);
			points[index].numVertices = vertices.size() / 3;
			for (int m = 0; m < vertices.size() / 3; m++)
			{
				Eigen::Vector3d vert = { vertices[m * 3],  vertices[m * 3 + 1],  vertices[m * 3 + 2] };
				points[index].verticsCoor.push_back(vert);
			}

			c.neighbors(neighbour);
			points[index].neighbour = neighbour;
			points[index].numFaces = neighbour.size();



			c.face_vertices(verticesFace);
			int start = 1, end = verticesFace[0] + 1;
			do
			{
				std::vector<int> faceVert;
				for (int m = start; m < end; m++)
				{
					faceVert.push_back(verticesFace[m]);
				}
				points[index].verticesFace.push_back(faceVert);
				start = end + 1;
				end += verticesFace[end] + 1;
			} while (points[index].verticesFace.size() != neighbour.size());

			c.normals(normals);
			for (int m = 0; m < normals.size() / 3; m++)
			{
				Eigen::Vector3d normal = { normals[m * 3],  normals[m * 3 + 1],  normals[m * 3 + 2] };
				points[index].surfaceNormal.push_back(normal);
			}

			

		} while (cl.inc());
	}









    //***********Read all particles and calculate the damage gradient*************//
    std::vector<Particle> particles;
    std::map<int, int> gridMap;
    std::vector<Grid> nodesVec;
    readParticlesAndCalGradient(&fullyDamagedParticlesNodesVec, &particles, param, &gridMap, &nodesVec);

    double pi = 3.141592653;
    double radius = param.dx * sqrt(3); // support radius of two points



    // find crack surface that does fully cut
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<int>> faces;
    std::vector<set<int>> sharedFace(points.size()); // store the point index which shares a crack surface with the other point
    // find faces that are in the interior and store neighbour information
    for (int i = 0; i < points.size(); i++) {
        Eigen::Vector3d pos = points[i].pos;
        if (ifCriticalNode(pos, param, &criticalNodeIndex) == false) {
            for (int k = 0; k < points[i].neighbour.size(); k++) {
                int neighbourIndex = points[i].neighbour[k];
                if (neighbourIndex > 0) // remove bounding box faces
                {

                    Eigen::Vector3d posNeig = points[neighbourIndex].pos;
                    Eigen::Vector3d posDiff = pos - posNeig;
                    double distancePair = posDiff.norm();

                    if (ifCriticalNode(posNeig, param, &criticalNodeIndex) == false) {
                        if (count(points[i].neighbourCalculated.begin(), points[i].neighbourCalculated.end(), neighbourIndex) == 0) // if these two neighbours have not yet been compared
                        {

                            if (distancePair > radius) // if their distance is larger than the threshold
                            {

                                double existInCrack = ifFullyDamaged((pos + posNeig) / 2.0, param, &gridMap, &nodesVec);
                                if (existInCrack >= param.damageThreshold) // if the middle point is located in the crack area
                                {
                                    std::vector<int> faceVerteice;
                                    for (int ver = 0; ver < points[i].verticesFace[k].size(); ver++) {
                                        int indexVertex = points[i].verticesFace[k][ver];
                                        Eigen::Vector3d vertex = points[i].verticsCoor[indexVertex];
                                        int vertexIndex = findIndexVertex(vertex, &vertices);
                                        faceVerteice.push_back(vertexIndex);
                                    }
                                    faces.push_back(faceVerteice);

                                    sharedFace[i].insert(neighbourIndex);
									sharedFace[neighbourIndex].insert(i);

                                    points[i].neighbourOtherSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourOtherSide.push_back(i);

                                } else {
                                    points[i].neighbourSameSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourSameSide.push_back(i);
                                }

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);

                            } else {
                                points[i].neighbourSameSide.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourSameSide.push_back(i);

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);
                            }
                        }
                    }
                }
            }
        }
    }

    // find faces that are in the interior and store neighbour information
    for (int i = 0; i < points.size(); i++) {
        Eigen::Vector3d pos = points[i].pos;

        if (ifCriticalNode(pos, param, &criticalNodeIndex) == false) {
            for (int k = 0; k < points[i].neighbour.size(); k++) {
                int neighbourIndex = points[i].neighbour[k];
                if (neighbourIndex > 0) // remove bounding box faces
                {

                    Eigen::Vector3d posNeig = points[neighbourIndex].pos;
                    Eigen::Vector3d posDiff = pos - posNeig;
                    double distancePair = posDiff.norm();

                    if (ifCriticalNode(posNeig, param, &criticalNodeIndex) == true) {
                        if (count(points[i].neighbourCalculated.begin(), points[i].neighbourCalculated.end(), neighbourIndex) == 0) // if these two neighbours have not yet been compared
                        {

                            if (distancePair > radius) // if their distance is larger than the threshold
                            {
                                //cout << "i = " << i << "; neighbourIndex = " << neighbourIndex << endl;
                                bool twoSide1 = ifTwoSides(i, neighbourIndex, &points, &boundaryNodesID, &boundaryNodesPosIndex, &pointIndexFind, param, &criticalNodeIndex);

                                //twoSide1 = true;
                                if (twoSide1 == true)
                                {
                                    std::vector<int> faceVerteice;
                                    for (int ver = 0; ver < points[i].verticesFace[k].size(); ver++) {
                                        int indexVertex = points[i].verticesFace[k][ver];
                                        Eigen::Vector3d vertex = points[i].verticsCoor[indexVertex];
                                        int vertexIndex = findIndexVertex(vertex, &vertices);
                                        faceVerteice.push_back(vertexIndex);
                                    }
                                    faces.push_back(faceVerteice);

                                    sharedFace[i].insert(neighbourIndex);
									sharedFace[neighbourIndex].insert(i);

                                    points[i].neighbourOtherSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourOtherSide.push_back(i);
                                } else {
                                    points[i].neighbourSameSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourSameSide.push_back(i);
                                }

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);

                            } else {
                                points[i].neighbourSameSide.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourSameSide.push_back(i);

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);
                            }
                        }
                    }
                }
            }
        }
    }

    // find faces that defined by critical point
    for (int i = 0; i < points.size(); i++) {

        Eigen::Vector3d pos = points[i].pos;
        if (ifCriticalNode(pos, param, &criticalNodeIndex) == true) {

            for (int k = 0; k < points[i].neighbour.size(); k++) {
                int neighbourIndex = points[i].neighbour[k];
                if (neighbourIndex > 0) // remove bounding box faces
                {
                    Eigen::Vector3d posNeig = points[neighbourIndex].pos;
                    Eigen::Vector3d posDiff = pos - posNeig;
                    double distancePair = posDiff.norm();

                    if (count(points[i].neighbourCalculated.begin(), points[i].neighbourCalculated.end(), neighbourIndex) == 0) // if these two neighbours have not yet been compared
                    {

                        if (distancePair > radius) // if their distance is larger than the threshold
                        {

                            bool twoSide1 = ifTwoSides(i, neighbourIndex, &points, &boundaryNodesID, &boundaryNodesPosIndex, &pointIndexFind, param, &criticalNodeIndex);

                            if (twoSide1 == true)
                            {
                                std::vector<int> faceVerteice;
                                for (int ver = 0; ver < points[i].verticesFace[k].size(); ver++) {
                                    int indexVertex = points[i].verticesFace[k][ver];
                                    Eigen::Vector3d vertex = points[i].verticsCoor[indexVertex];
                                    int vertexIndex = findIndexVertex(vertex, &vertices);
                                    faceVerteice.push_back(vertexIndex);
                                }
                                faces.push_back(faceVerteice);


                                sharedFace[i].insert(neighbourIndex);
                                sharedFace[neighbourIndex].insert(i);

                            }

                            // store computation information.
                            points[i].neighbourCalculated.push_back(neighbourIndex);
                            points[neighbourIndex].neighbourCalculated.push_back(i);
                        }
                    }
                }
            }
        }
    }




	////////////////////////////////
	// calculate all fragments
	////////////////////////////////
	std::set<int> remainingPoints;
	for (int i = 0; i < points.size(); i++)
	{
		remainingPoints.insert(i);
	}

	std::vector<std::vector<int>> allFragments;

	do
	{

		std::vector<int> fragment;
		int firstElement = *remainingPoints.begin();
		fragment.push_back(firstElement);
		remainingPoints.erase(firstElement);
		int start = 0; // the index of starting pointof a single layer
		int countLayer = 1;
		do
		{
			int countLayerUpdate = 0;
			for (int i = start; i < start + countLayer; i++) // parse each candidate
			{
				int currentPoint = fragment[i];
				for (int j = 0; j < points[currentPoint].neighbour.size(); j++)
				{
					int candPoint = points[currentPoint].neighbour[j];
					if (candPoint >= 0)
					{
						if (sharedFace[currentPoint].find(candPoint) == sharedFace[currentPoint].end())
						{
							if (std::find(fragment.begin(), fragment.end(), candPoint) == fragment.end())
							{
								fragment.push_back(candPoint);
								remainingPoints.erase(candPoint);
								countLayerUpdate += 1;
							}

						}
					}

				}
			}
			start = start + countLayer;
			countLayer = countLayerUpdate;

		} while (countLayer != 0);

		allFragments.push_back(fragment);



	} while (remainingPoints.size() != 0);






	// add interior volume to the surrounding fragment
	std::vector<std::vector<int>> allFragmentsRemoveInterior; // store each fragment where interior volume is added to surrounding volume
	std::set<int> interiorFragmentIndex; // store the index of fragment which is an interior volume
	for (int i = 0; i < allFragments.size(); i++)
	{
		std::vector<int> fragment = allFragments[i];
		bool interior = true;
		for (int j = 0; j < fragment.size(); j++)
		{
			for (int h = 0; h < points[fragment[j]].neighbour.size(); h++)
			{
				int neigPoint = points[fragment[j]].neighbour[h];
				if (neigPoint < 0)
				{
					interior = false;
					break;
				}
			}

			if (interior == false)
			{
				break;
			}
		}

		if (interior == true)
		{
			interiorFragmentIndex.insert(i);
		}


	}




	if (interiorFragmentIndex.size() == 0) // all fargments are not interior
	{
		allFragmentsRemoveInterior = allFragments;
	}
	else
	{
		std::map<int, int> mergeInterior; // key and value are the index of interior volume and surrounding volume respectively

		// find the surrounding volume that an interior volume should be added into
		std::set<int>::iterator it;
		for (it = interiorFragmentIndex.begin(); it != interiorFragmentIndex.end(); ++it)
		{
			// find neighbour point in the surrounding volume
			int surroundPoint = -99; // surrounding neighbour point of this interior volume
			for (int k = 0; k < allFragments[*it].size(); k++)
			{
				int seedPoint = allFragments[*it][k];
				for (int i = 0; i < points[seedPoint].neighbour.size(); i++)
				{
					int neigPoint = points[seedPoint].neighbour[i];
					if (sharedFace[seedPoint].find(neigPoint) != sharedFace[seedPoint].end())
					{
						surroundPoint = neigPoint;
						break;
					}
				}

				if (surroundPoint >= 0)
				{
					break;
				}
			}


			// find surrounding volume
			int surroundingVolume = -99;
			for (int i = 0; i < allFragments.size(); i++)
			{
				if (std::find(allFragments[i].begin(), allFragments[i].end(), surroundPoint) != allFragments[i].end())
				{
					surroundingVolume = i;
					break;
				}
			}
			mergeInterior[*it] = surroundingVolume;
            if(surroundingVolume < 0)
            {       
                crackSurfaceInfor refinedMeshResult;
                return refinedMeshResult;

            }
            ASSERT(surroundingVolume >= 0);

		}


		// add interior volume into the surrounding volume
		std::map<int, int>::iterator itMap;
		for (itMap = mergeInterior.begin(); itMap != mergeInterior.end(); itMap++)
		{
			int interiorFragment = itMap->first;
			int surroundingFragment = itMap->second;
			allFragments[surroundingFragment].insert(allFragments[surroundingFragment].end(), allFragments[interiorFragment].begin(), allFragments[interiorFragment].end());
		}


		for (int i = 0; i < allFragments.size(); i++)
		{
			if (mergeInterior.find(i) == mergeInterior.end())
			{
				allFragmentsRemoveInterior.push_back(allFragments[i]);
			}
		}


	}






    crackSurfaceInfor refinedMeshResult;
    refinedMeshResult.numOfFragments = allFragmentsRemoveInterior.size();
    refinedMeshResult.vertices = vertices;
    refinedMeshResult.faces = faces;
    for(int i = 0; i < faces.size(); i++)
    {
        Eigen::Vector3d v1 = vertices[faces[i][0]];
        Eigen::Vector3d v2 = vertices[faces[i][1]];
        Eigen::Vector3d v3 = vertices[faces[i][2]];
        Eigen::Vector3d normal = (v1 - v2).cross(v2 - v3).normalized();
        refinedMeshResult.facesNormal.push_back(normal);

    }

    return refinedMeshResult;


}




////////////////////////////////
// function for getting partial cracks
////////////////////////////////
meshFaceRefined extractCrackSurface_PartialCracks(std::vector<Particle>* particlesRaw, std::vector<Eigen::Vector3d>* mpmParticlePositions, struct parametersSim param)
{
    SCOPED_TIMER(__FUNCTION__);

    printf("%s\n", __FUNCTION__);
    //*********Read fully damaged particles***********//
    std::vector<Particle> fullyDamagedParticles;
    std::map<int, int> fullyDamagedParticlesGridMap;
    std::vector<Grid> fullyDamagedParticlesNodesVec;

    readParticles(particlesRaw, &fullyDamagedParticles, true, param);
    calDamageGradient(&fullyDamagedParticles, param, param.dx, &fullyDamagedParticlesGridMap, &fullyDamagedParticlesNodesVec);
    setNodeValue(&fullyDamagedParticlesNodesVec, 1);
    findBoundaryNodes(&fullyDamagedParticles, &fullyDamagedParticlesNodesVec, &fullyDamagedParticlesGridMap, param, 2); // This bounding box is used to generate boundary nodes
    findBoundaryNodes(&fullyDamagedParticles, &fullyDamagedParticlesNodesVec, &fullyDamagedParticlesGridMap, param, 3); // This bounding box is used to generate clip surface mesh
    //*********Read fully damaged particles***********//



    //*********Read all particles***********//
    std::vector<Particle> allParticles;
    std::map<int, int> allParticlesGridMap;
    std::vector<Grid> allParticlesNodesVec;

    readParticles(particlesRaw, &allParticles, false, param);
    calDamageGradient(&allParticles, param, param.dx, &allParticlesGridMap, &allParticlesNodesVec);
    setNodeValue(&allParticlesNodesVec, 1);
    findBoundaryNodes(&allParticles, &allParticlesNodesVec, &allParticlesGridMap, param, 2); // This bounding box is used to generate boundary nodes
    findBoundaryNodes(&allParticles, &allParticlesNodesVec, &allParticlesGridMap, param, 3); // This bounding box is used to generate clip surface mesh
    //*********Read all particles***********//


    // find boundary nodes
    std::vector<int> allParticlesNodeBoundaryIndex; // store the index or ID of allParticles boundary nodes
    for (int i = 0; i < allParticlesNodesVec.size(); i++) {
        if (allParticlesNodesVec[i].Di == 2) {
            int ID = calculateID(allParticlesNodesVec[i].posIndex[0], allParticlesNodesVec[i].posIndex[1], allParticlesNodesVec[i].posIndex[2], param.length, param.dx);
            allParticlesNodeBoundaryIndex.push_back(ID);
        }
    }


    //************find surface nodes//
    std::vector<int> surfaceNodesID; // store surface nodes ID
    std::vector<Eigen::Vector3i> surfaceNodesPosIndex; // store surface position index
    std::vector<int> boundaryNodesIDNoClean; // store boundary nodes ID
    std::vector<Eigen::Vector3i> boundaryNodesPosIndexNoClean; // store boundary position index
    for (int i = 0; i < fullyDamagedParticlesNodesVec.size(); i++) {
        if (fullyDamagedParticlesNodesVec[i].Di == 2) {
            int ID = calculateID(fullyDamagedParticlesNodesVec[i].posIndex[0], fullyDamagedParticlesNodesVec[i].posIndex[1], fullyDamagedParticlesNodesVec[i].posIndex[2], param.length, param.dx);
            if (count(allParticlesNodeBoundaryIndex.begin(), allParticlesNodeBoundaryIndex.end(), ID) == 0) {
                boundaryNodesPosIndexNoClean.push_back(fullyDamagedParticlesNodesVec[i].posIndex);
                boundaryNodesIDNoClean.push_back(ID);
            }

            //************find surface nodes//
            if (count(allParticlesNodeBoundaryIndex.begin(), allParticlesNodeBoundaryIndex.end(), ID) != 0) {
                int surfaceNodeID = calculateID(fullyDamagedParticlesNodesVec[i].posIndex[0], fullyDamagedParticlesNodesVec[i].posIndex[1], fullyDamagedParticlesNodesVec[i].posIndex[2], param.length, param.dx);
                surfaceNodesID.push_back(surfaceNodeID);
                surfaceNodesPosIndex.push_back(fullyDamagedParticlesNodesVec[i].posIndex);
            }
            //************find surface nodes//
        }
    }

    // clean isolate nodes and sharp nodes
    std::vector<int> boundaryNodesID; // store boundary nodes ID
    std::vector<Eigen::Vector3i> boundaryNodesPosIndex; // store boundary position index
    for (int m = 0; m < boundaryNodesPosIndexNoClean.size(); m++) {
        for (int axis = 0; axis < 3; axis++) // three axis directions
        {
            for (int pn = 0; pn < 2; pn++) {
                Eigen::Vector3i normal = { 0, 0, 0 };
                normal[axis] = 2 * pn - 1;
                Eigen::Vector3i nodePosIndex = boundaryNodesPosIndexNoClean[m] + normal;
                int nodeID = calculateID(nodePosIndex[0], nodePosIndex[1], nodePosIndex[2], param.length, param.dx);
                if (count(boundaryNodesIDNoClean.begin(), boundaryNodesIDNoClean.end(), nodeID) != 0) {
                    boundaryNodesID.push_back(boundaryNodesIDNoClean[m]);
                    boundaryNodesPosIndex.push_back(boundaryNodesPosIndexNoClean[m]);
                    goto stopNotIsolate;
                }
            }
        }

    stopNotIsolate:
        int a;
    }



    // define a std::map that can find the index of a point in the point vector
    std::map<int, int> pointIndexFind;
    for (int m = 0; m < boundaryNodesID.size(); m++) {
        pointIndexFind[boundaryNodesID[m]] = m;
    }



    ofstream outfile45("boundaryNodes.txt", ios::trunc);
    for (int m = 0; m < boundaryNodesPosIndex.size(); m++) {
        outfile45 << m << " " << scientific << setprecision(8) << boundaryNodesPosIndex[m][0] * param.dx << " " << boundaryNodesPosIndex[m][1] * param.dx << " " << boundaryNodesPosIndex[m][2] * param.dx << endl;
    }
    outfile45.close();



    // find critical nodes
    std::vector<int> criticalNodeIndex;
    int criticalNodeVolumeLength = 1;
    for (int m = 0; m < surfaceNodesPosIndex.size(); m++) {

        for (int i = -criticalNodeVolumeLength; i < criticalNodeVolumeLength + 1; i++) {
            for (int j = -criticalNodeVolumeLength; j < criticalNodeVolumeLength + 1; j++) {
                for (int k = -criticalNodeVolumeLength; k < criticalNodeVolumeLength + 1; k++) {
                    Eigen::Vector3i increment = { i, j, k };
                    Eigen::Vector3i nodePosIndex = surfaceNodesPosIndex[m] + increment;
                    int ID = calculateID(nodePosIndex[0], nodePosIndex[1], nodePosIndex[2], param.length, param.dx);
                    if (count(boundaryNodesID.begin(), boundaryNodesID.end(), ID) != 0) {
                        if (count(criticalNodeIndex.begin(), criticalNodeIndex.end(), ID) == 0) {
                            criticalNodeIndex.push_back(ID);
                        }
                    }
                }
            }
        }
    }



    double x_min = 1.0E10, x_max = -1.0E10;
    double y_min = 1.0E10, y_max = -1.0E10;
    double z_min = 1.0E10, z_max = -1.0E10;
    int n_x = 10, n_y = 10, n_z = 10;

    for (int m = 0; m <  (*mpmParticlePositions).size(); m++) 
    {
        x_min = std::min(x_min, (*mpmParticlePositions)[m][0]);
        y_min = std::min(y_min, (*mpmParticlePositions)[m][1]);
        z_min = std::min(z_min, (*mpmParticlePositions)[m][2]);

        x_max = std::max(x_max, (*mpmParticlePositions)[m][0]);
        y_max = std::max(y_max, (*mpmParticlePositions)[m][1]);
        z_max = std::max(z_max, (*mpmParticlePositions)[m][2]);
    }
    x_min = x_min - 40 * param.dx;
    x_max = x_max + 40 * param.dx;
    y_min = y_min - 40 * param.dx;
    y_max = y_max + 40 * param.dx;
    z_min = z_min - 40 * param.dx;
    z_max = z_max + 40 * param.dx;



	pre_container pcon(x_min, x_max, y_min, y_max, z_min, z_max, false, false, false);
	// Import the particles from a file
	pcon.import("boundaryNodes.txt");
	pcon.guess_optimal(n_x, n_y, n_z);
	container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, false, false, false, 8);
	pcon.setup(con);


	std::vector<Point> points;
	for (int i = 0; i < boundaryNodesPosIndex.size(); i++)
	{
		Eigen::Vector3d pos = { 0 , 0 , 0 }; // each point's position
		std::vector<Eigen::Vector3d> verticsCoor; // vertices' coordinate
		std::vector<std::vector<int>> verticesFace; // vertices of each face
		std::vector<Eigen::Vector3d> surfaceNormal; // vertices' coordinate
		std::vector<int> neighbour; // neighbour points that share common faces with this point
		std::vector<int> neighbourCalculated; // neighbour points that have already find shared face
		points.push_back(Point(-999, pos, 0, verticsCoor, 0, verticesFace, surfaceNormal, neighbour));
	}


    wallShell wallFace(&fullyDamagedParticlesGridMap, &fullyDamagedParticlesNodesVec, param, -99);
    con.add_wall(wallFace);
	c_loop_all cl(con);
	if (cl.start())
	{
        
		double x, y, z;
		int index;
		voronoicell_neighbor c;
		std::vector<int> neighbour, verticesFace;
		std::vector<double> vertices;
		std::vector<double> normals;
		do if (con.compute_cell(c, cl))
		{

			int index = cl.pid();
			points[index].index = index;

			cl.pos(x, y, z);
			points[index].pos = { x, y, z };

			c.vertices(x, y, z, vertices);
			points[index].numVertices = vertices.size() / 3;
			for (int m = 0; m < vertices.size() / 3; m++)
			{
				Eigen::Vector3d vert = { vertices[m * 3],  vertices[m * 3 + 1],  vertices[m * 3 + 2] };
				points[index].verticsCoor.push_back(vert);
			}

			c.neighbors(neighbour);
			points[index].neighbour = neighbour;
			points[index].numFaces = neighbour.size();



			c.face_vertices(verticesFace);
			int start = 1, end = verticesFace[0] + 1;
			do
			{
				std::vector<int> faceVert;
				for (int m = start; m < end; m++)
				{
					faceVert.push_back(verticesFace[m]);
				}
				points[index].verticesFace.push_back(faceVert);
				start = end + 1;
				end += verticesFace[end] + 1;
			} while (points[index].verticesFace.size() != neighbour.size());

			c.normals(normals);
			for (int m = 0; m < normals.size() / 3; m++)
			{
				Eigen::Vector3d normal = { normals[m * 3],  normals[m * 3 + 1],  normals[m * 3 + 2] };
				points[index].surfaceNormal.push_back(normal);
			}

			

		} while (cl.inc());
	}






    //***********Read all particles and calculate the damage gradient*************//
    std::vector<Particle> particles;
    std::map<int, int> gridMap;
    std::vector<Grid> nodesVec;
    readParticlesAndCalGradient(&fullyDamagedParticlesNodesVec, &particles, param, &gridMap, &nodesVec);

    double pi = 3.141592653;
    double radius = param.dx * sqrt(3); // support radius of two points



    std::vector<Eigen::Vector3d> vertices; // vertex index

    cout << "Start extracting interior faces" << endl;

    std::vector<std::vector<int>> faces;
    // find faces that are in the interior and store neighbour information
    for (int i = 0; i < points.size(); i++) {
        Eigen::Vector3d pos = points[i].pos;
        if (ifCriticalNode(pos, param, &criticalNodeIndex) == false) {
            for (int k = 0; k < points[i].neighbour.size(); k++) {
                int neighbourIndex = points[i].neighbour[k];
                if (neighbourIndex > 0) // remove bounding box faces
                {

                    Eigen::Vector3d posNeig = points[neighbourIndex].pos;
                    Eigen::Vector3d posDiff = pos - posNeig;
                    double distancePair = posDiff.norm();

                    if (ifCriticalNode(posNeig, param, &criticalNodeIndex) == false) {
                        if (count(points[i].neighbourCalculated.begin(), points[i].neighbourCalculated.end(), neighbourIndex) == 0) // if these two neighbours have not yet been compared
                        {

                            if (distancePair > radius) // if their distance is larger than the threshold
                            {

                                double existInCrack = ifFullyDamaged((pos + posNeig) / 2.0, param, &gridMap, &nodesVec);
                                if (existInCrack >= param.damageThreshold) // if the middle point is located in the crack area
                                {

                                    std::vector<int> faceVerteice;
                                    for (int ver = 0; ver < points[i].verticesFace[k].size(); ver++) {
                                        int indexVertex = points[i].verticesFace[k][ver];
                                        Eigen::Vector3d vertex = points[i].verticsCoor[indexVertex];
                                        int vertexIndex = findIndexVertex(vertex, &vertices);
                                        faceVerteice.push_back(vertexIndex);
                                    }
                                    faces.push_back(faceVerteice);

                                    points[i].neighbourOtherSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourOtherSide.push_back(i);

                                } else {
                                    points[i].neighbourSameSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourSameSide.push_back(i);
                                }

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);

                            } else {
                                points[i].neighbourSameSide.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourSameSide.push_back(i);

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);
                            }
                        }
                    }
                }
            }
        }
    }

    cout << "Finish extracting interior faces" << endl;

    cout << "Start extracting boundary faces" << endl;
    // find faces that are in the interior and store neighbour information
    for (int i = 0; i < points.size(); i++) {
        Eigen::Vector3d pos = points[i].pos;

        if (ifCriticalNode(pos, param, &criticalNodeIndex) == false) {
            for (int k = 0; k < points[i].neighbour.size(); k++) {
                int neighbourIndex = points[i].neighbour[k];
                if (neighbourIndex > 0) // remove bounding box faces
                {

                    Eigen::Vector3d posNeig = points[neighbourIndex].pos;
                    Eigen::Vector3d posDiff = pos - posNeig;
                    double distancePair = posDiff.norm();

                    if (ifCriticalNode(posNeig, param, &criticalNodeIndex) == true) {
                        if (count(points[i].neighbourCalculated.begin(), points[i].neighbourCalculated.end(), neighbourIndex) == 0) // if these two neighbours have not yet been compared
                        {

                            if (distancePair > radius) // if their distance is larger than the threshold
                            {
                                //cout << "i = " << i << "; neighbourIndex = " << neighbourIndex << endl;
                                bool twoSide1 = ifTwoSides(i, neighbourIndex, &points, &boundaryNodesID, &boundaryNodesPosIndex, &pointIndexFind, param, &criticalNodeIndex);

                                //twoSide1 = true;
                                if (twoSide1 == true) {
                                    std::vector<int> faceVerteice;
                                    for (int ver = 0; ver < points[i].verticesFace[k].size(); ver++) {
                                        int indexVertex = points[i].verticesFace[k][ver];
                                        Eigen::Vector3d vertex = points[i].verticsCoor[indexVertex];
                                        int vertexIndex = findIndexVertex(vertex, &vertices);
                                        faceVerteice.push_back(vertexIndex);
                                    }
                                    faces.push_back(faceVerteice);

                                    points[i].neighbourOtherSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourOtherSide.push_back(i);
                                } else {
                                    points[i].neighbourSameSide.push_back(neighbourIndex);
                                    points[neighbourIndex].neighbourSameSide.push_back(i);
                                }

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);

                            } else {
                                points[i].neighbourSameSide.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourSameSide.push_back(i);

                                // store computation information.
                                points[i].neighbourCalculated.push_back(neighbourIndex);
                                points[neighbourIndex].neighbourCalculated.push_back(i);
                            }
                        }
                    }
                }
            }
        }
    }

    // find faces that defined by critical point
    for (int i = 0; i < points.size(); i++) {

        Eigen::Vector3d pos = points[i].pos;
        if (ifCriticalNode(pos, param, &criticalNodeIndex) == true) {

            for (int k = 0; k < points[i].neighbour.size(); k++) {
                int neighbourIndex = points[i].neighbour[k];
                if (neighbourIndex > 0) // remove bounding box faces
                {
                    Eigen::Vector3d posNeig = points[neighbourIndex].pos;
                    Eigen::Vector3d posDiff = pos - posNeig;
                    double distancePair = posDiff.norm();

                    if (count(points[i].neighbourCalculated.begin(), points[i].neighbourCalculated.end(), neighbourIndex) == 0) // if these two neighbours have not yet been compared
                    {

                        if (distancePair > radius) // if their distance is larger than the threshold
                        {

                            bool twoSide1 = ifTwoSides(i, neighbourIndex, &points, &boundaryNodesID, &boundaryNodesPosIndex, &pointIndexFind, param, &criticalNodeIndex);

                            if (twoSide1 == true) {
                                std::vector<int> faceVerteice;
                                for (int ver = 0; ver < points[i].verticesFace[k].size(); ver++) {
                                    int indexVertex = points[i].verticesFace[k][ver];
                                    Eigen::Vector3d vertex = points[i].verticsCoor[indexVertex];
                                    int vertexIndex = findIndexVertex(vertex, &vertices);
                                    faceVerteice.push_back(vertexIndex);
                                }
                                faces.push_back(faceVerteice);
                            }

                            // store computation information.
                            points[i].neighbourCalculated.push_back(neighbourIndex);
                            points[neighbourIndex].neighbourCalculated.push_back(i);
                        }
                    }
                }
            }
        }
    }

    cout << "Finish extracting boundary faces" << endl;


    std::cout << "The number of crack faces is " << faces.size() << endl;

 

    for (int h = 0; h < faces.size(); h++) {
        for (int g = 0; g < faces[h].size(); g++) {
            faces[h][g] = faces[h][g] - 1;
        }
    }


   




    struct meshFaceRefined refinedMeshResult;
    if (vertices.size() != 0 && faces.size() != 1) {
        tuple<std::vector<Eigen::Vector3d>, std::vector<vector<int>>, std::vector<vector<int>>> result = cleanMeshAfterVoro(&vertices, &faces);

        std::ofstream outfile21("vertices.txt", std::ios::trunc);
        for (int k = 0; k < std::get<0>(result).size(); k++) {
            Eigen::Vector3d scale = std::get<0>(result)[k];
            outfile21 << std::scientific << std::setprecision(8) << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
        }
        outfile21.close();

        std::ofstream outfile22("faces.txt", std::ios::trunc);
        for (int k = 0; k < std::get<1>(result).size(); k++) {
            outfile22 << std::scientific << std::setprecision(8) << "f ";
            for (int j = 0; j < std::get<1>(result)[k].size(); j++) {
                outfile22 << std::scientific << std::setprecision(8) << std::get<1>(result)[k][j] << "  ";
            }
            outfile22 << std::endl;
        }
        outfile22.close();

        std::ofstream outfile23("verticesInfor.txt", std::ios::trunc);
        for (int k = 0; k < std::get<2>(result).size(); k++) {
            outfile23 << std::scientific << std::setprecision(8) << "f ";
            for (int j = 0; j < std::get<2>(result)[k].size(); j++) {
                outfile23 << std::scientific << std::setprecision(8) << std::get<2>(result)[k][j] << "  ";
            }
            outfile23 << std::endl;
        }
        outfile23.close();

        std::cout << "refine crack mesh" << std::endl;

        

        
        try 
        {
            refinedMeshResult = meshRefine(&std::get<0>(result), &std::get<1>(result), &std::get<2>(result), param.drx);
        }
        catch (...) 
        {
            std::cout << "Failed to run triangle.h library" << std::endl;
        }

  
    }





    return  refinedMeshResult;



}









////////////////////////////////
// cut objects with ftetwild
////////////////////////////////
void cutObject(std::string nameObject, std::vector<meshObjFormat>* fragments, std::vector<meshObjFormat>* fragmentCollision)
{


	for (int i = 0; i < (*fragments).size(); i++)
	{
        std::string frag_name = nameObject + "_Frag_" + std::to_string(i);
        writeOffFile((*fragments)[i].vertices, (*fragments)[i].faces, frag_name);
		string comm = "./FloatTetwild_bin -i " + frag_name +  ".off --log ftetwild.log --is-quiet --level 0 --manifold-surface -o " + nameObject + "_Frag_" + std::to_string(i) + "_Res";
		FILE* pf = popen(comm.c_str(), "r");



	}


    // detect if all fragments are fixed
	for (int i = 0; i < (*fragments).size(); i++)
	{

        std::string fixedFragName = nameObject + "_Frag_" + std::to_string(i) + "_Res__sf.obj";
        bool existFixFrag = false;
        do
        {
            if(FILE* file = fopen(fixedFragName.c_str(), "r"))
            {
                existFixFrag = true;
            }
        } while (existFixFrag == false);
        
	}


    // delete unnecessary files
    for (int i = 0; i < (*fragments).size(); i++)
	{
        // delete unnecessary files
        std::string n0 = nameObject + "_Frag_" + std::to_string(i) + "_Res_.csv";
        std::string n1 = nameObject + "_Frag_" + std::to_string(i) + "_Res_.msh";
        std::string n2 = nameObject + "_Frag_" + std::to_string(i) + "_Res__tracked_surface.stl";
        // int del0 = remove(n0.c_str());
        // int del1 = remove(n1.c_str());
        // int del2 = remove(n2.c_str());
	}



     // delete unnecessary files
    for (int i = 0; i < (*fragments).size(); i++)
	{
        std::string fixedFragName = nameObject + "_Frag_" + std::to_string(i) + "_Res__sf.obj";

		std::ofstream outfile("bool" + std::to_string(i) + ".json", std::ios::trunc); 
		outfile << "{" << std::endl;
		outfile << "    \"operation\": \"intersection\"," << std::endl;
		outfile << "    \"left\": \"" + fixedFragName +"\"," << std::endl;
		outfile << "    \"right\": \"" + nameObject +".obj\"" << std::endl;
		outfile << "}" << std::endl;
		outfile.close();


        std::string boolFile = "bool" + std::to_string(i) + ".json";
        bool existFixFrag = false;
        do
        {
            if(FILE* file = fopen(boolFile.c_str(), "r"))
            {
                existFixFrag = true;
            }
        } while (existFixFrag == false);


		string comm2 = "./FloatTetwild_bin --csg bool" + std::to_string(i) + ".json" +  " --log ftetwild.log --manifold-surface --is-quiet --level 0 -o cuttedFrag" + std::to_string(i);
	
		FILE* pf2 = popen(comm2.c_str(), "r");
	}




    // detect if all fragments are cutted
	for (int i = 0; i < (*fragments).size(); i++)
	{

        std::string fixedFragName = "cuttedFrag" + std::to_string(i) + "__sf.obj";
        bool existFixFrag = false;
        do
        {
            if(FILE* file = fopen(fixedFragName.c_str(), "r"))
            {
                existFixFrag = true;
            }
        } while (existFixFrag == false);
        
	}


    for (int i = 0; i < (*fragments).size(); i++)
    {
        std::string bo = "cuttedFrag" + std::to_string(i) + "__sf.obj";
        meshObjFormat res = readObj(bo);
        (*fragmentCollision).push_back(res);
    }


    // delete unnecessary files
    for (int i = 0; i < (*fragments).size(); i++)
	{
        std::string n_sf = nameObject + "_Frag_" + std::to_string(i) + "_Res__sf.obj";
        //int del_n_sf = remove(n_sf.c_str());

        // delete unnecessary files
        std::string n0 = "cuttedFrag" + std::to_string(i) + "_.csv";
        std::string n1 = "cuttedFrag" + std::to_string(i) + "_.msh";
        std::string n2 = "cuttedFrag" + std::to_string(i) + "__tracked_surface.stl";
        std::string n3 = "cuttedFrag" + std::to_string(i) + "__sf.obj";
        // int del0 = remove(n0.c_str());
        // int del1 = remove(n1.c_str());
        // int del2 = remove(n2.c_str());
        // int del3 = remove(n3.c_str());


        std::string n_json = "bool" + std::to_string(i) + ".json";
        int del_n_json = remove(n_json.c_str());

	}

}


#define mcCheckError(errCode) mcCheckError_(errCode, __FILE__, __LINE__)

void mcDebugOutput(McDebugSource source,
    McDebugType type,
    unsigned int id,
    McDebugSeverity severity,
    size_t length,
    const char* message,
    const void* userParam)
{
    printf("MCUT LOG: %s\n", message);
#if 0
    printf("Debug message ( %d ), length=%zu\n%s\n--\n", id, length, message);
    printf("userParam=%p\n", userParam);

    switch (source) {
    case MC_DEBUG_SOURCE_API:
        printf("Source: API");
        break;
    case MC_DEBUG_SOURCE_KERNEL:
        printf("Source: Kernel");
        break;
    case MC_DEBUG_SOURCE_ALL:
        break;
    }

    printf("\n");

    switch (type) {
    case MC_DEBUG_TYPE_ERROR:
        printf("Type: Error");
        break;
    case MC_DEBUG_TYPE_DEPRECATED_BEHAVIOR:
        printf("Type: Deprecated Behaviour");
        break;
    case MC_DEBUG_TYPE_OTHER:
        printf("Type: Other");
        break;
    case MC_DEBUG_TYPE_ALL:
        break;
    }

    printf("\n");

    switch (severity) {
    case MC_DEBUG_SEVERITY_HIGH:
        printf("Severity: high");
        break;
    case MC_DEBUG_SEVERITY_MEDIUM:
        printf("Severity: medium");
        break;
    case MC_DEBUG_SEVERITY_LOW:
        printf("Severity: low");
        break;
    case MC_DEBUG_SEVERITY_NOTIFICATION:
        printf("Severity: notification");
        break;
    case MC_DEBUG_SEVERITY_ALL:
        break;
    }

    printf("\n\n");
    #endif
}

void mcCheckError_(McResult err, const char* file, int line)
{
    std::string error;
    switch (err) {
    case MC_OUT_OF_MEMORY:
        error = "MC_OUT_OF_MEMORY";
        break;
    case MC_INVALID_VALUE:
        error = "MC_INVALID_VALUE";
        break;
    case MC_INVALID_OPERATION:
        error = "MC_INVALID_OPERATION";
        break;
    case MC_NO_ERROR:
        error = "MC_NO_ERROR";
        break;
        case MC_RESULT_MAX_ENUM:
        error = "UNKNOWN";
        break;
    }
    if (err) {
        std::cout << error << " | " << file << " (" << line << ")" << std::endl;
    }
}

////////////////////////////////
// intersection of two objects with mcut
////////////////////////////////
std::vector<std::pair<meshObjFormat, std::pair<std::set<int> , std::map<int, int> > > > intersection_MCUT(meshObjFormat* collisionMeshParent, meshObjFormat* fragmentVolume)
{
    std::vector<std::pair<meshObjFormat, std::pair<std::set<int> , std::map<int, int> > > > results; // 1. fragment; 2. vertex index of those who are on the crack surface; 3. vertex index of those who are on the src mesh


    /////////////////////////////////
    // read parent collision mesh
    /////////////////////////////////
    InputMesh srcMesh;

    // copy vertices
    for (int i = 0; i < (int)collisionMeshParent->vertices.size(); ++i)
    {
        srcMesh.vertexCoordsArray.push_back(collisionMeshParent->vertices[i][0]);
        srcMesh.vertexCoordsArray.push_back(collisionMeshParent->vertices[i][1]);
        srcMesh.vertexCoordsArray.push_back(collisionMeshParent->vertices[i][2]);
    }

    // copy faces
    for (int i = 0; i < (int)collisionMeshParent->faces.size(); ++i)
    {
        for (int j = 0; j < (int)collisionMeshParent->faces[i].size(); ++j)
        {
            srcMesh.faceIndicesArray.push_back(collisionMeshParent->faces[i][j]);
        }

        srcMesh.faceSizesArray.push_back((uint32_t)collisionMeshParent->faces[i].size());
    }

    printf("src mesh:\n\tvertices=%d\n\tfaces=%d\n", (int)srcMesh.vertexCoordsArray.size(), (int)srcMesh.faceIndicesArray.size());



    /////////////////////////////////
    // read fragment mesh
    /////////////////////////////////
    InputMesh cutMesh;

    // copy vertices
    for (int i = 0; i < (int)fragmentVolume->vertices.size(); ++i)
    {
        cutMesh.vertexCoordsArray.push_back(fragmentVolume->vertices[i][0]);
        cutMesh.vertexCoordsArray.push_back(fragmentVolume->vertices[i][1]);
        cutMesh.vertexCoordsArray.push_back(fragmentVolume->vertices[i][2]);
    }

    // copy faces
    for (int i = 0; i < (int)fragmentVolume->faces.size(); ++i)
    {
        for (int j = 0; j < (int)fragmentVolume->faces[i].size(); ++j)
        {
            cutMesh.faceIndicesArray.push_back(fragmentVolume->faces[i][j]);
        }

        cutMesh.faceSizesArray.push_back((uint32_t)fragmentVolume->faces[i].size());
    }


    printf("cut mesh:\n\tvertices=%d\n\tfaces=%d\n", (int)cutMesh.vertexCoordsArray.size(), (int)cutMesh.faceIndicesArray.size());

    // create a context
    // -------------------
    McContext context = MC_NULL_HANDLE;
    McResult err = mcCreateContext(&context, MC_DEBUG);

    ASSERT(err == MC_NO_ERROR);
    // config debug output
    // -----------------------
    uint64_t numBytes = 0;
    McFlags contextFlags;
    err = mcGetInfo(context, MC_CONTEXT_FLAGS, 0, nullptr, &numBytes);
    mcCheckError(err);

    ASSERT(sizeof(McFlags) == numBytes);

    err = mcGetInfo(context, MC_CONTEXT_FLAGS, numBytes, &contextFlags, nullptr);
    mcCheckError(err);

    if (contextFlags & MC_DEBUG) {
        mcDebugMessageCallback(context, mcDebugOutput, nullptr);
        mcDebugMessageControl(context, McDebugSource::MC_DEBUG_SOURCE_ALL, McDebugType::MC_DEBUG_TYPE_ALL, McDebugSeverity::MC_DEBUG_SEVERITY_ALL, true);
    }



    //  do the cutting (boolean ops)
    // -------------------------------
    // printf("\nInputs: \n\tShape A = 'cube.obj'.\n\tShape B = 'torus.obj'\n\n");

    // We can either let MCUT compute all possible meshes (including patches etc.), or we can
    // constrain the library to compute exactly the boolean op mesh we want. This 'constrained' case
    // is done with the following flags.
    // NOTE: you can extend these flags by bitwise ORing with additional flags (see `McDispatchFlags' in mcut.h)
    const std::map<std::string, McFlags> booleanOps = {{"INTERSECTION", MC_DISPATCH_FILTER_FRAGMENT_SEALING_INSIDE | MC_DISPATCH_FILTER_FRAGMENT_LOCATION_BELOW}};


    const McFlags boolOpFlags = MC_DISPATCH_FILTER_FRAGMENT_SEALING_INSIDE | MC_DISPATCH_FILTER_FRAGMENT_LOCATION_BELOW;
    const std::string boolOpName = "INTERSECTION";


    //std::cout<<"inner mcut = "<< 0<<std::endl;

    err = mcDispatch(
        context,
        MC_DISPATCH_VERTEX_ARRAY_DOUBLE |          // vertices are in array of doubles
            MC_DISPATCH_ENFORCE_GENERAL_POSITION | // perturb if necessary
            MC_DISPATCH_INCLUDE_VERTEX_MAP | //
            MC_DISPATCH_FILTER_FRAGMENT_SEALING_INSIDE | MC_DISPATCH_FILTER_FRAGMENT_LOCATION_BELOW,  // INTERSECTION                         // filter flags which specify the type of output we want
        // source mesh
        reinterpret_cast<const void *>(srcMesh.vertexCoordsArray.data()),
        reinterpret_cast<const uint32_t *>(srcMesh.faceIndicesArray.data()),
        srcMesh.faceSizesArray.data(),
        static_cast<uint32_t>(srcMesh.vertexCoordsArray.size() / 3),
        static_cast<uint32_t>(srcMesh.faceSizesArray.size()),
        // cut mesh
        reinterpret_cast<const void *>(cutMesh.vertexCoordsArray.data()),
        cutMesh.faceIndicesArray.data(),
        cutMesh.faceSizesArray.data(),
        static_cast<uint32_t>(cutMesh.vertexCoordsArray.size() / 3),
        static_cast<uint32_t>(cutMesh.faceSizesArray.size()));

    ASSERT(err == MC_NO_ERROR);



    //std::cout<<"inner mcut = "<< 1<<std::endl;

    // query the number of available connected component
    // --------------------------------------------------
    uint32_t numConnComps;

    err = mcGetConnectedComponents(context, MC_CONNECTED_COMPONENT_TYPE_ALL, 0, NULL, &numConnComps);
    ASSERT(err == MC_NO_ERROR);

    printf("MCUT available connected components: %u\n", numConnComps);

    err = mcGetConnectedComponents(context, MC_CONNECTED_COMPONENT_TYPE_FRAGMENT, 0, NULL, &numConnComps);
    ASSERT(err == MC_NO_ERROR);


    if (numConnComps == 0)
    {
        fprintf(stdout, "No connected components found by MCUT, which is strange\n");
        std::getc(stdin);
        return results;
    }


    //ASSERT(numConnComps == 1); // exactly 1 result (for this example)

    std::vector<McConnectedComponent> connectedComponents(numConnComps, MC_NULL_HANDLE);
    connectedComponents.resize(numConnComps);
    err = mcGetConnectedComponents(context, MC_CONNECTED_COMPONENT_TYPE_FRAGMENT, (uint32_t)connectedComponents.size(), connectedComponents.data(), NULL);

    ASSERT(err == MC_NO_ERROR);

    // query the data of each connected component from MCUT
    // -------------------------------------------------------

    
    for(uint32_t ccn = 0; ccn < numConnComps; ccn++)
    {
        McConnectedComponent connComp = connectedComponents[ccn];

        // query the vertices
        // ----------------------
        uint64_t numBytes = 0;
        err = mcGetConnectedComponentData(context, connComp, MC_CONNECTED_COMPONENT_DATA_VERTEX_DOUBLE, 0, NULL, &numBytes);
        ASSERT(err == MC_NO_ERROR);
        uint32_t ccVertexCount = (uint32_t)(numBytes / (sizeof(double) * 3));
        std::vector<double> ccVertices((uint64_t)ccVertexCount * 3u, 0);
        err = mcGetConnectedComponentData(context, connComp, MC_CONNECTED_COMPONENT_DATA_VERTEX_DOUBLE, numBytes, (void *)ccVertices.data(), NULL);
        ASSERT(err == MC_NO_ERROR);

        // query the faces
        // -------------------
        numBytes = 0;
        err = mcGetConnectedComponentData(context, connComp, MC_CONNECTED_COMPONENT_DATA_FACE, 0, NULL, &numBytes);
        ASSERT(err == MC_NO_ERROR);
        std::vector<uint32_t> ccFaceIndices(numBytes / sizeof(uint32_t), 0);
        err = mcGetConnectedComponentData(context, connComp, MC_CONNECTED_COMPONENT_DATA_FACE, numBytes, ccFaceIndices.data(), NULL);
        ASSERT(err == MC_NO_ERROR);
        // query the face sizes
        // ------------------------
        numBytes = 0;
        err = mcGetConnectedComponentData(context, connComp, MC_CONNECTED_COMPONENT_DATA_FACE_SIZE, 0, NULL, &numBytes);
        ASSERT(err == MC_NO_ERROR);
        std::vector<uint32_t> faceSizes(numBytes / sizeof(uint32_t), 0);
        err = mcGetConnectedComponentData(context, connComp, MC_CONNECTED_COMPONENT_DATA_FACE_SIZE, numBytes, faceSizes.data(), NULL);
        ASSERT(err == MC_NO_ERROR);
        // query the face map
        // ------------------------
        const uint32_t ccFaceCount = static_cast<uint32_t>(faceSizes.size());

        /// ------------------------------------------------------------------------------------

        // Here we show, how to know when connected components, pertain particular boolean operations.

        McPatchLocation patchLocation = (McPatchLocation)0;

        err = mcGetConnectedComponentData(context, connComp, MC_CONNECTED_COMPONENT_DATA_PATCH_LOCATION, sizeof(McPatchLocation), &patchLocation, NULL);
        ASSERT(err == MC_NO_ERROR);
        McFragmentLocation fragmentLocation = (McFragmentLocation)0;
        err = mcGetConnectedComponentData(context, connComp, MC_CONNECTED_COMPONENT_DATA_FRAGMENT_LOCATION, sizeof(McFragmentLocation), &fragmentLocation, NULL);
        ASSERT(err == MC_NO_ERROR);

        err = mcGetConnectedComponentData(context, connComp, MC_CONNECTED_COMPONENT_DATA_VERTEX_MAP, 0, NULL, &numBytes);
        ASSERT(err == MC_NO_ERROR);
        std::vector<uint32_t> ccVertexMap;
        ccVertexMap.resize(numBytes / sizeof(uint32_t));
        err = mcGetConnectedComponentData(context, connComp, MC_CONNECTED_COMPONENT_DATA_VERTEX_MAP, numBytes, ccVertexMap.data(), NULL);
        ASSERT(err == MC_NO_ERROR);



        // save cc mesh to .obj file
        // -------------------------

        meshObjFormat intersectionVol;
        std::set<int> crackSurfaceVertices; // vertices which are located on the crack surface
        std::map<int, int> parentSurfaceVertices; // vertices which are located on the src surface
        // write vertices and normals
        for (uint32_t i = 0; i < ccVertexCount; ++i)
        {
            double x = ccVertices[(uint64_t)i * 3 + 0];
            double y = ccVertices[(uint64_t)i * 3 + 1];
            double z = ccVertices[(uint64_t)i * 3 + 2];
            Eigen::Vector3d vert = {x, y, z};
            intersectionVol.vertices.push_back(vert);

            const int ccVertexIdx = i;
            // input mesh (source mesh or cut mesh) vertex index (which may be offsetted)
            const uint32_t imVertexIdxRaw = ccVertexMap.at(ccVertexIdx);
            bool vertexIsFromSrcMesh = (imVertexIdxRaw < (std::uint32_t)srcMesh.V.size());
            const bool isSeamVertex = (imVertexIdxRaw == MC_UNDEFINED_VALUE);
            uint32_t imVertexIdx = imVertexIdxRaw; // actual index value, accounting for offset


            if (!vertexIsFromSrcMesh && !isSeamVertex)
            {
                imVertexIdx = (imVertexIdxRaw - (std::uint32_t)srcMesh.V.size()); // account for offset
                crackSurfaceVertices.insert(i);
            }
            else
            {
                parentSurfaceVertices[i] = imVertexIdxRaw;
            }
            
        }
        int faceVertexOffsetBase = 0;

        std::cout<< " crackSurfaceVertices.size() = "<<crackSurfaceVertices.size()<<std::endl;
        std::cout<< " parentSurfaceVertices.size() = "<<parentSurfaceVertices.size()<<std::endl;


        // for each face in CC
        for (uint32_t f = 0; f < ccFaceCount; ++f)
        {
            bool reverseWindingOrder = (fragmentLocation == MC_FRAGMENT_LOCATION_BELOW) && (patchLocation == MC_PATCH_LOCATION_OUTSIDE);
            int faceSize = faceSizes.at(f);

            std::vector<int> face;
            // for each vertex in face
            for (int v = (reverseWindingOrder ? (faceSize - 1) : 0);
                    (reverseWindingOrder ? (v >= 0) : (v < faceSize));
                    v += (reverseWindingOrder ? -1 : 1))
            {
                face.push_back(ccFaceIndices[(uint64_t)faceVertexOffsetBase + v]);
            } 
            intersectionVol.faces.push_back(face);
            faceVertexOffsetBase += faceSize;
        }
        std::pair<std::set<int> , std::map<int, int> > crack_surface_vertices = std::make_pair(crackSurfaceVertices, parentSurfaceVertices);
        std::pair<meshObjFormat, std::pair<std::set<int> , std::map<int, int> > > res = std::make_pair(intersectionVol, crack_surface_vertices);

        results.push_back(res);

        
    }
        

    // 6. free connected component data
    // --------------------------------
    err = mcReleaseConnectedComponents(context, (uint32_t)connectedComponents.size(), connectedComponents.data());
    ASSERT(err == MC_NO_ERROR);



    // 7. destroy context
    // ------------------
    err = mcReleaseContext(context);

    ASSERT(err == MC_NO_ERROR);

    return results;
}



////////////////////////////////
// triangulate a polygon mesh with libigl
////////////////////////////////
void triangulate_polygonMesh(meshObjFormat* polygonMesh)
{
    std::vector<std::vector<int>> faces;
    
    int numOfFaces = (int)polygonMesh->faces.size();
    // NOTE: "curveMesh" contains only triangulated faces because we only
    // use it to store intersection points. Intersecting faces are removed
    for (int i = 0; i < numOfFaces; ++i)
    {
        std::vector<int> face = polygonMesh->faces[i];
        int faceDim = (int)face.size();
        if (faceDim == 3)
        {
            faces.push_back(face);
        }
        else
        {
            Eigen::MatrixXd V(faceDim, 3);
            Eigen::VectorXi I(faceDim), C(2);
            C(0) = 0;
            C(1) = faceDim;
            std::map<int, int> faceToMeshMap;
            for (int j = 0; j < faceDim; ++j)
            {
                I(j) = j;
                Eigen::Vector3d vert = polygonMesh->vertices[polygonMesh->faces[i][j]];
                V.row(j) = vert;
                faceToMeshMap[j] = polygonMesh->faces[i][j];
            }

            // Convert polygon representation to triangles
            Eigen::MatrixXi F;
            Eigen::VectorXi J;
            igl::polygons_to_triangles(I, C, F, J);


            for (int j = 0; j < F.rows(); ++j)
            {
                std::vector<int> faceAdded;
                faceAdded.push_back(faceToMeshMap[F(j,0)]);
                faceAdded.push_back(faceToMeshMap[F(j,1)]);
                faceAdded.push_back(faceToMeshMap[F(j,2)]);
                faces.push_back(faceAdded);
            }


        }


    }

    polygonMesh->faces = faces;
   

}


////////////////////////////////
// upsample a triangle mesh from cleaned ftetwild output
////////////////////////////////
meshObjFormat upsampleMesh(parametersSim param, meshObjFormat *inputMesh)
{
    meshObjFormat outputMesh;
    
    // approximate the upsample step
    double averageEdgeLength = 0;
    int countNum = 0;
    for(int m = 0; m <(*inputMesh).faces.size(); m++)
    {
        double edge0 = ((*inputMesh).vertices[(*inputMesh).faces[m][0]] - (*inputMesh).vertices[(*inputMesh).faces[m][1]]).norm();
        double edge1 = ((*inputMesh).vertices[(*inputMesh).faces[m][1]] - (*inputMesh).vertices[(*inputMesh).faces[m][2]]).norm(); 
        double edge2 = ((*inputMesh).vertices[(*inputMesh).faces[m][2]] - (*inputMesh).vertices[(*inputMesh).faces[m][0]]).norm(); 
        averageEdgeLength += edge0; 
        averageEdgeLength += edge1; 
        averageEdgeLength += edge2; 
        countNum += 3; 
    }
    averageEdgeLength = averageEdgeLength / countNum;


    if(averageEdgeLength <= param.drx)
    {
        outputMesh.vertices = (*inputMesh).vertices;
        outputMesh.faces = (*inputMesh).faces;
    }
    else
    {
        int sampleStep = trunc(log(averageEdgeLength) / log(2));
        sampleStep = 3;
        // upsample a triangle mesh
        Eigen::MatrixXd V((*inputMesh).vertices.size() ,3);
        Eigen::MatrixXi F((*inputMesh).faces.size(), 3);
        for(int m = 0; m <( *inputMesh).vertices.size(); m++)
        {
            V(m , 0) = (*inputMesh).vertices[m][0];
            V(m , 1) = (*inputMesh).vertices[m][1];
            V(m , 2) = (*inputMesh).vertices[m][2];
        }
        for(int m = 0; m <( *inputMesh).faces.size(); m++)
        {
            F(m , 0) = (*inputMesh).faces[m][0];
            F(m , 1) = (*inputMesh).faces[m][1];
            F(m , 2) = (*inputMesh).faces[m][2];
        }
        igl::upsample(V, F, sampleStep);

        // get the upsampled mesh
        for(int m = 0; m <V.rows(); m++)
        {
            outputMesh.vertices.push_back(V.row(m));
        }
        for(int m = 0; m <F.rows(); m++)
        {
            std::vector<int> face;
            face.push_back(F(m,0));
            face.push_back(F(m,1));
            face.push_back(F(m,2));
            outputMesh.faces.push_back(face);
        }


    }

  
    return outputMesh;

}



////////////////////////////////
// cut objects with mcut
////////////////////////////////
// 1. param: parameter; 2. objectName: name of the object to be cut; 3. collisionMeshParent: the object to be cut
// 4. fragments: each fragment; 5. 
void cutObject_MCUT(parametersSim param, std::string objectName, meshObjFormat* collisionMeshParent, std::vector<meshObjFormat>* fragments, std::vector<meshObjFormat>* finalFragments)
{

    std::cout << "Expected fragments = "<< (*fragments).size() <<std::endl;

    std::cout<<"mcut = "<< "1"<<std::endl;
//#pragma omp parallel for num_threads(param.numOfMpmThreads)
	for (int i = 0; i < (*fragments).size(); i++)
	{
        std::string frag_name = objectName + "_Frag_" + std::to_string(i);
        writeOffFile((*fragments)[i].vertices, (*fragments)[i].faces, frag_name);
		string comm = "./FloatTetwild_bin -i " + frag_name +  ".off --log ftetwild.log --manifold-surface -l 0.025 --is-quiet --level 0 -o " + objectName + "_Frag_" + std::to_string(i) + "_Res > " + objectName + "-tetwild.log.";
		std::cout << "running: " << comm << std::endl;
        int ret = std::system(comm.c_str());
        if(ret != 0)
        {
            printf("failed to run ftetwild");
            std::abort();
        }
        //FILE* pf = popen(comm.c_str(), "w");
        //int res = pclose(pf);
	}
#if 0
    std::cout<<"mcut = "<< "2"<<std::endl;
    // detect if all fragments are fixed
	for (int i = 0; i < (*fragments).size(); i++)
	{

        std::string fixedFragName = objectName + "_Frag_" + std::to_string(i) + "_Res__sf.obj";
        bool existFixFrag = false;
        do
        {
            if(FILE* file = fopen(fixedFragName.c_str(), "r"))
            {
                existFixFrag = true;
            }
        } while (existFixFrag == false);
        
	}

    std::cout<<"mcut = "<< "3"<<std::endl;
    // delete unnecessary files
    for (int i = 0; i < (*fragments).size(); i++)
	{
        // delete unnecessary files
        std::string n0 = objectName + "_Frag_" + std::to_string(i) + "_Res_.csv";
        std::string n1 = objectName + "_Frag_" + std::to_string(i) + "_Res_.msh";
        std::string n2 = objectName + "_Frag_" + std::to_string(i) + "_Res__tracked_surface.stl";
        int del0 = remove(n0.c_str());
        int del1 = remove(n1.c_str());
        int del2 = remove(n2.c_str());

        std::string b0 = objectName + "_Frag_" + std::to_string(i) + ".off";
        int delb0 = remove(b0.c_str());
	}
#endif
    std::cout<<"mcut = "<< "4"<<std::endl;
    for(int i = 0; i < (*fragments).size(); i++)
    {
        std::string pathVol = objectName + "_Frag_" + std::to_string(i) + "_Res__sf.obj";
        //std::cout<<"mcut = "<< "4.1"<<std::endl;
        std::cout<<pathVol <<std::endl;
        meshObjFormat fragmentVolume_ = readObj(pathVol);
        meshObjFormat fragmentVolume_pre =  cleanFtetWild(&fragmentVolume_);
        meshObjFormat fragmentVolume = upsampleMesh(param, &fragmentVolume_pre);

        //std::cout<<"mcut = "<< "4.2"<<std::endl;

        std::string name11 = "srcMesh_" + std::to_string(i);
		writeOffFile((*collisionMeshParent).vertices, (*collisionMeshParent).faces, name11);

        std::string name22 = "cutMesh_" + std::to_string(i) ;
		writeOffFile(fragmentVolume.vertices, fragmentVolume.faces, name22);

        
        std::vector<std::pair<meshObjFormat, std::pair<std::set<int> , std::map<int, int> > > > intersectioVol = intersection_MCUT(collisionMeshParent, &fragmentVolume);
        std::cout<<"mcut = "<< "4.3"<<std::endl;

        if(intersectioVol.size() != 0)
        {
            for(int kj = 0; kj < intersectioVol.size(); kj++)
            {
                (*finalFragments).push_back(intersectioVol[kj].first);
                // triangulate_polygonMesh(&intersectioVol[kj].first);
                // (*fragmentCollision).push_back(intersectioVol[kj].first);
                // (*crackSurVertices).push_back(intersectioVol[kj].second);

            }
        }

    }
    
#if 0
    std::cout<<"mcut = "<< "5"<<std::endl;
     // delete unnecessary files
    for (int i = 0; i < (*fragments).size(); i++)
	{
        std::string fixedFragName = objectName + "_Frag_" + std::to_string(i) + "_Res__sf.obj";
        int del_fixedFrag = remove(fixedFragName.c_str());
	}
#endif

}
















