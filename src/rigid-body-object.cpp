
#include "mpm-fracture/object.h"

#if 0
// MeshDecimation
#ifndef MeshDecimation_h
#define MeshDecimation_h
#include "Simplify.h"
#endif
#endif // MeshDecimation

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>    // level set CSG operations
#include <openvdb/tools/LevelSetUtil.h> // segmentSDF(...)

#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/decimate.h>
#include <igl/edge_flaps.h>
#include <igl/parallel_for.h>
#include <igl/qslim.h>
#include <igl/shortest_edge_and_midpoint.h>

void Object::initVolumeMeshData()
{

    renderingMesh->need_curvatures();
    renderingMesh->need_adjacentfaces();
    renderingMesh->need_normals();
    renderingMesh->need_neighbors();
    renderingMesh->need_across_edge();

    trimeshHighResKDTree = new trimesh::KDtree(renderingMesh->vertices);

    trimeshHighResMaxEdgeLen = renderingMesh->stat(trimesh::TriMesh::STAT_MAX, trimesh::TriMesh::STAT_EDGELEN);

    ASSERT(trimeshHighResMaxEdgeLen > 0);

    for (int i = 0; i < (int)renderingMesh->faces.size(); ++i)
    {
        const trimesh::TriMesh::Face &f = renderingMesh->faces[i];
        float area = trimesh::len(trimesh::trinorm(renderingMesh->vertices[f[0]], renderingMesh->vertices[f[1]], renderingMesh->vertices[f[2]]));
        elemArea[i] = area;
        trimesh::point c = (renderingMesh->vertices[f[0]] + renderingMesh->vertices[f[1]] + renderingMesh->vertices[f[2]]) / 3.0f;
        elemCtr[i] = toEigen(c);
    }
}

void Object::storeData(ColliderData collMesh)
{
    storedData.push_back(collMesh);
}

/* notify the BulletWrapper that this objects has a new collision shape
 * that should be updated, the wrapper must first remove this objects RB
 * from the dynamics world, then call doCollisionUpdate() and finally
 * add this object's RB back into the dynamics world
 */
bool Object::pendingCollisionUpdate() { return collUpdate.shape != NULL; }

void Object::doCollisionUpdate()
{
    if (pendingCollisionUpdate())
    {
        printf("\n%% processing collision shape update ... ");

        if (collInUse.shape)
        {
            collInUse.deleteAll();
        }

        collInUse.set(collUpdate);
        collUpdate.setNull();

        rigidBodyPtr->setCollisionShape(collInUse.shape);
        // update mass and inertia ...
        double mass = updateVolume * m_params->density;
        btVector3 inertia;
        collInUse.shape->calculateLocalInertia(mass, inertia);
        rigidBodyPtr->setMassProps(mass, inertia);

        // correct for the change in COM between the old and new collision shapes
        //Eigen::Vector3f com = updateCOM - vdbCOM;
        Eigen::Vector3f com = updateCOM - centerOfMass;
        btVector3 btCOM(com[0], com[1], com[2]); // assuming the COM was at the origin in the original mesh, this is the COM shift
        rigidBodyPtr->getWorldTransform().setOrigin(
            rigidBodyPtr->getWorldTransform().getOrigin() + rigidBodyPtr->getWorldTransform().getBasis() * btCOM);
        //vdbCOM = updateCOM;
        centerOfMass = updateCOM;
        printf("\n%% collision shape update done ... ");
        printf("\n%% !!! new rest.coeff. is %.3lf (%s)\n", rigidBodyPtr->getRestitution(), ((std::string *)rigidBodyPtr->getUserPointer())->c_str());
        printf("%% new margin is %.3lg\n", rigidBodyPtr->getCollisionShape()->getMargin());
    }
}

float Object::getEovNuSq()
{
    return m_params->youngsModulus / (1.0 - (m_params->poissonsRatio * m_params->poissonsRatio));
}

float Object::getCurvatureAndElementID(const Eigen::Vector3f &point, int &closestFaceIndex)
{
    ASSERT(renderingMesh != nullptr);
    ASSERT(iglRenderingMesh != nullptr);

#if 1
    closestFaceIndex = -1;
    

    typedef Eigen::Matrix<float, 1, 3> RowVector;

    RowVector p (3);
    p(0) = point(0);
    p(1) = point(1);
    p(2) = point(2);
    RowVector closestPointOnSurface(3);
    //printf("iglRenderingMesh->V iglRenderingMesh->F\n", (int)iglRenderingMesh->V.rows(), (int)iglRenderingMesh->F.rows());
    float dist = iglRenderingMesh->aabbTree.squared_distance(iglRenderingMesh->V, iglRenderingMesh->F, p, closestFaceIndex, closestPointOnSurface);
    //printf("closest point dist = %f closestFaceIndex = %d\n", dist, closestFaceIndex);
    ASSERT(closestFaceIndex != -1);

    //ASSERT(closestFaceIndex < renderingMesh->faces.size());
     Eigen::Vector3i face = iglRenderingMesh->F.row(closestFaceIndex);
     //printf("facet = %d %d %d\n", face[0], face[1], face[2]);
    // calculate the mean curviture on the closest triangle
    const float curvature = (renderingMesh->curv1[face[0]] + renderingMesh->curv1[face[1]] + renderingMesh->curv1[face[2]]) / 3.0f;

    //printf("curvature = %f\n", curvature);
    return curvature;
#if 0
    Eigen::MatrixXf L; // barycentric coordinates

    float curvature = 0;

    // compute barycentric coords for each closest point on the closest tri
    {
        Eigen::MatrixXf P(1, 3);
        P.row(0) = point;
        Eigen::Vector3i tri = iglRenderingMesh->F.row(closestFaceIndex);
        Eigen::MatrixXf A(1, 3);
        A.row(0) = iglRenderingMesh->V.row(tri(0));
        Eigen::MatrixXf B(1, 3);
        B.row(0) = iglRenderingMesh->V.row(tri(1));
        Eigen::MatrixXf C(1, 3);
        C.row(0) = iglRenderingMesh->V.row(tri(2));

        igl::barycentric_coordinates(P, A, B, C, L);

        const Eigen::Vector3f b = L.row(0);

        curvature = iglRenderingMesh->meanCurvature(tri(0)) * b(0) + //
                    iglRenderingMesh->meanCurvature(tri(1)) * b(1) + //
                    iglRenderingMesh->meanCurvature(tri(2)) * b(2);

        return curvature;
    }
#endif
#else

    // find the closest triangle to "point"
    //std::cout << "sdgasgdasgdkjgda" << std::endl;
    closestFaceIndex = getClosestFace(renderingMesh, trimeshHighResKDTree, toTrimesh(point), trimeshHighResMaxEdgeLen);
    //std::cout << closestFaceIndex << std::endl;
    ASSERT(closestFaceIndex != -1);

    const trimesh::TriMesh::Face &face = renderingMesh->faces[closestFaceIndex];

    // calculate the mean curviture on the closest triangle
    const float curvature = (renderingMesh->curv1[face[0]] + renderingMesh->curv1[face[1]] + renderingMesh->curv1[face[2]]) / 3.0f;
    return curvature;
#endif
}

int Object::getClosestFace(const trimesh::TriMesh *m, const trimesh::KDtree *kd, const trimesh::point &p, float maxdist2)
{
    // const float *match = kd->closest_to_pt(p, maxdist2);
    // The closest vertex on the triangle might be much further away than the closest
    // point on the surface.  So, we need to be conservative here.

    // "p" comes from the collision-shape mesh which, after decimation, can slightly change shape (may become smaller
    // than rendering mesh).
    // This change of shape can sometimes imply that "p" is more than "maxdist2" (mean edge length in redering mesh) away from the closest point
    // on the rendering mesh. To resolve this we scale "maxdist2" by "factor".
    this->renderingMesh->need_bbox();
    //trimesh::vec3 bbsize = this->renderingMesh->bbox.size();
    float factor = 100;
    // "match" is a pointer to an element
    //std::cout << p << std::endl;
    const float *match = kd->closest_to_pt(p, factor);

    if (!match)
    {
        fprintf(stderr, "failed to find closest point (p=(%f %f %f), maxdist2=%f)\n", p[0], p[1], p[2], maxdist2);
        return -1;
    }
    // compute vertex index using pointer arithmetic ("&(m->vertices[0][0])" is the address
    // of x-comp of first vertex)
    int ind = (match - (const float *)&(m->vertices[0][0])) / 3;
    int nv = m->vertices.size();
    if (ind < 0 || ind >= nv)
    {
        fprintf(stderr, "index out of range (p=(%f %f %f), maxdist2=%f)\n", p[0], p[1], p[2], maxdist2);
        fprintf(stderr, "index = %d, the number of vertices = %d", ind, nv);
        return -1;
    }

    // get adjacent faces
    const std::vector<int> &a = m->adjacentfaces[ind];
    if (a.empty())
    {
        fprintf(stderr, "closest point has zero adj faces (p=(%f %f %f), maxdist2=%f)\n", p[0], p[1], p[2], maxdist2);
        return -1;
    }

    // now find the closest triangle by searching the triangles
    // which are adjacent to the closest point we found

    int closestTriangle = -1;
    float closest_dist2 = factor;

    for (size_t i = 0; i < a.size(); i++)
    {
        trimesh::point c = getClosestPointOnTriangle(m, a[i], p);
        float this_dist2 = trimesh::dist2(c, p);
        if (this_dist2 < closest_dist2)
        {
            closest_dist2 = this_dist2;
            closestTriangle = a[i];
        }
    }

    if (closestTriangle == -1)
    {
        fprintf(stderr, "could not find closest tri from adj list (p=(%f %f %f), maxdist2=%f)\n", p[0], p[1], p[2], maxdist2);
    }

    return closestTriangle;
}

trimesh::point Object::getClosestPointOnTriangle(const trimesh::TriMesh *m, int i, const trimesh::point &p)
{
    const trimesh::TriMesh::Face &f = m->faces[i];
    const trimesh::point &v0 = m->vertices[f[0]];
    const trimesh::point &v1 = m->vertices[f[1]];
    const trimesh::point &v2 = m->vertices[f[2]];
    trimesh::vec3 a = v1 - v0;
    trimesh::vec3 b = v2 - v0;
    trimesh::vec3 p1 = p - v0;
    trimesh::vec3 n = a CROSS b;

    float A[3][3] = {{a[0], b[0], n[0]},
                     {a[1], b[1], n[1]},
                     {a[2], b[2], n[2]}};
    float x[3] = {p1[0], p1[1], p1[2]};
    int indx[3];
    trimesh::ludcmp<float, 3>(A, indx);
    trimesh::lubksb<float, 3>(A, indx, x);

    if (x[0] >= 0.0f && x[1] >= 0.0f && x[0] + x[1] <= 1.0f)
        return v0 + x[0] * a + x[1] * b;

    trimesh::point c01 = closestPointOnSegment(v0, v1, p);
    trimesh::point c12 = closestPointOnSegment(v1, v2, p);
    trimesh::point c20 = closestPointOnSegment(v2, v0, p);
    float d01 = trimesh::dist2(c01, p);
    float d12 = trimesh::dist2(c12, p);
    float d20 = trimesh::dist2(c20, p);
    if (d01 < d12)
    {
        if (d01 < d20)
            return c01;
        else
            return c20;
    }
    else
    {
        if (d12 < d20)
            return c12;
        else
            return c20;
    }
}

// Helper for dist2mesh:
// Find closest point to p on segment from v0 to v1
trimesh::point Object::closestPointOnSegment(const trimesh::point &v0, const trimesh::point &v1, const trimesh::point &p)
{
    trimesh::vec3 v01 = v1 - v0;
    float d = (p - v0) DOT v01;
    d /= trimesh::len2(v01);
    if (d < 0.0f)
        d = 0.0f;
    else if (d > 1.0f)
        d = 1.0f;
    return v0 + d * v01;
}

void Object::clearContacts()
{
    balancedTractions.clear();
    contactTractions.clear();
    contactDuration = 0.0;
}

/* Add a contact point which will contribute a boundary condition for the fracture simulation
 * p is the point of impact in this rigid body's local coordinate system (alternatively a triangle-ID can be specified)
 * d is the direction of the impact (also in local coords), should point towards the inside for collisions
 * the simulation will run for as many time-steps as required by the contact with the longest duration
 */

void Object::addContact(unsigned int tri, Eigen::Vector3f &d, float impulse, float duration)
{
    //printf("%% ... %s contact (%.3lg %.3lg %.3lg), tri %d\n",((string*)rb->getUserPointer())->c_str(),p[0],p[1],p[2],tri);

    if (duration == 0)
    {
        return; // prevent division by zero
    }

    Eigen::Vector3f traction(d);
    traction.normalize(); // d should be normalized by caller, but just to be sure ...
    traction *= impulse / duration;

    traction /= elemArea[tri];

    if (contactTractions.count(tri) == 0)
    { // triangle has not tractions applied to it yet (in current timestep)
        contactTractions[tri].setZero();
    }

    contactTractions[tri] += traction; // accumulation traction forces

    if (duration > contactDuration)
    {
        contactDuration = duration; // store the max. contact duration for the sim-run
    }

    // printf("%% ... contact on tri %6d - traction is (%.3lg %.3lg %.3lg)\n", tri, traction[0], traction[1], traction[2]);
}

/* Check if the current set of contact points is sufficient to run a fracture simulation
 * returns 0 if the longest contact duration is too short for the time-step size
 * otherwise returns the total magnitude of deformational forces
 * the caller should then decide whether the force magnitude warrants running the simulation
 * This function must be called after "calculateBalancedContactTractions"
 */
float Object::getTotalContactForce(float fractureSimTimeStepSize)
{
    int steps = (int)((contactDuration / fractureSimTimeStepSize) + 0.5); // round
    if (steps == 0)
    {
        return 0.0;
    }

    float totalForce = 0.0;

    for (std::map<int, Eigen::Vector3f>::iterator it = balancedTractions.begin(); it != balancedTractions.end(); ++it)
    {
        totalForce += it->second.norm() * elemArea[it->first];
    }

    return totalForce;
}

// balance tractions --> goal is to have sum(forces)=0 and sum(torques)=0 over the surface
// reads from member variable contactTractions
// This function must be called before "getTotalContactForce"
void Object::calculateBalancedContactTractions()
{
    this->balancedTractions.clear();
    Eigen::Vector3f force(0.0, 0.0, 0.0);
    Eigen::Vector3f torque(0.0, 0.0, 0.0);

    // initialze q for all surface (i.e. non-crack) elements
    // also compute the net force and torque due to the contactTractions
    for (int i = 0; i < (int)renderingMesh->faces.size(); ++i)
    {
        if (contactTractions.count(i))
        { // does the face have a contact force/Traction?
            balancedTractions[i] = contactTractions[i];
            force += contactTractions[i] * elemArea[i];
            torque += contactTractions[i].cross(elemCtr[i]) * elemArea[i];
        }
        else
        {
            balancedTractions[i].setZero();
        }
    }

    //from now on iterate over q to get all elements that are not cracks

    // we want to find a traction field that is
    // -- as close as possible to the contactTractions (in L2-norm)
    // -- has 0 net force and
    // -- has 0 net torque (both integrated over the entire surface)
    // we find this by setting q = q_in - q_c (i.e. contactTractions - correction)
    // then do a constrained quadratic minimization on the correction
    // the function we want to minimize is the squared L2-norm of the correction q_c'*q_c
    // this results in a linear system which we solve via a Schur complement system
    // since our target function is so simple, and we have only 2 constraints
    // with 3 dimensions each (so 6 in all), we only need to invert a 6x6 matrix.

    // build the 6x6 Schur complement system S=B(A^-1)B' = -b, with b=[ force ; torque ]
    // actually A is an identity since we use the squared L2-norm so we ignore it from now on
    // but we could weight e.g. by element areas just as easily
    // as long a A is diagonal (with non-zero entries) we're good and S is symmetric
    Eigen::MatrixXf S(6, 6);
    S.setZero();
    float tmp;
    for (std::map<int, Eigen::Vector3f>::iterator it = balancedTractions.begin(); it != balancedTractions.end(); ++it)
    {
        // run once through the elements and compute the entries of S
        // S(0,0) is (force-x, force-x), S(1,1) is (force-y, force-y) etc.
        // S(0,1)=S(0,2)=S(1,2) = 0 -- forces do not mix
        // S(0,3)=S(1,4)=S(2,5) = 0 -- forces never affect torques in the same axis
        // S(1,3) is (force-y, torque-x), S(2,3) is (force-z,torque-x) etc.
        // we'll build S block-wise S=[ Sff  Sft ; Sft', Stt ]
        // actually we'll first only build the upper triangle of S and then transpose-copy them down

        // ok let's start with the forces, Sff = Bf*Bf' where Bf consists of
        // blocks that are a*I_3x3 with a the element area and I the identity
        // consequently Sff = sum(a^2)*I_3x3 i.e. S(0,0)=S(1,1)=S(2,2)=sum(a^2) the sum of the element areas squared
        tmp = elemArea[it->first];
        tmp *= tmp; // tmp is the squared area of the element
        S(0, 0) += tmp;
        S(1, 1) += tmp;
        S(2, 2) += tmp;

        // now do the torques, Stt = Bt*Bt' where Bt consists of blocks
        // that are a*M with M being the cross-product from the right operator
        // with the elements centroid: M*b = (b x elemCtr[i])
        // i.e. Stt=sum(a^2*M*M') with M=[0 mz -my ; -mz 0 mx ; my -mx 0], mx,my,mz being the components of elemCtr[i]
        // and  M*M' = [my^2+mz^2  -mx*my  -mx*mz ; -mx*my  mx^2+mz^2  -my*mz ; -mx*mz  -my*mz  mx^2+my^2]
        // first do the diagonal entries
        S(3, 3) += tmp * (elemCtr[it->first][1] * elemCtr[it->first][1] + elemCtr[it->first][2] * elemCtr[it->first][2]);
        S(4, 4) += tmp * (elemCtr[it->first][0] * elemCtr[it->first][0] + elemCtr[it->first][2] * elemCtr[it->first][2]);
        S(5, 5) += tmp * (elemCtr[it->first][0] * elemCtr[it->first][0] + elemCtr[it->first][1] * elemCtr[it->first][1]);
        // now the upper triangle
        S(3, 4) -= tmp * (elemCtr[it->first][0] * elemCtr[it->first][1]);
        S(3, 5) -= tmp * (elemCtr[it->first][0] * elemCtr[it->first][2]);
        S(4, 5) -= tmp * (elemCtr[it->first][1] * elemCtr[it->first][2]);

        // finally the force-torque interaction block Sft = Bf*Bt'
        // Sft = a*I*a*M' = a^2*M'
        // again we need only the upper triangle of this
        S(0, 4) += tmp * elemCtr[it->first][2];
        S(0, 5) -= tmp * elemCtr[it->first][1]; // note the minus!
        S(1, 5) += tmp * elemCtr[it->first][0];
    }
    // fill in the lower triangle of Stt
    S(4, 3) = S(3, 4);
    S(5, 3) = S(3, 5);
    S(5, 4) = S(4, 5);
    // now the lower triangle of Sft (this is a skew-symmetric block!)
    S(1, 3) = -S(0, 4);
    S(2, 3) = -S(0, 5);
    S(2, 4) = -S(1, 5);
    // and finally the 3x3 lower block of Sft'
    S(3, 1) = S(1, 3);
    S(3, 2) = S(2, 3);
    S(4, 2) = S(2, 4);
    S(4, 0) = S(0, 4);
    S(5, 0) = S(0, 5);
    S(5, 1) = S(1, 5);

    Eigen::VectorXf v(6);
    v[0] = force[0];
    v[1] = force[1];
    v[2] = force[2]; // note: no minus here!
    v[3] = torque[0];
    v[4] = torque[1];
    v[5] = torque[2];
    S.llt().solveInPlace(v);
    //Eigen::LLT<Eigen::MatrixXf> lltSolver(S); // Mathematica could pre-solve this analytically ... but it's a terribly long formula
    //lltSolver.solveInPlace(v);

    // now that we have the Schur complement solution
    // the traction correction is given by A*q_c = -B'*v
    // (note that we've already omitted a minus above so we also ignore this one)
    // again we'll ignore A as it is the identity in the case of the squared L2-norm,
    // but obviously we could use any positive diagonal weight very easily

    // compute and apply q_c
    Eigen::Vector3f q_ce; // per element traction correction
    for (std::map<int, Eigen::Vector3f>::iterator it = balancedTractions.begin(); it != balancedTractions.end(); ++it)
    {
        // as above, B consists of a 3x3 block Bf and another one Bt
        // so B'*v = [ Bf' Bt' ]*v = Bf'*v[0-2] + Bt'*v[3-5]
        // the result is the traction correction for this element

        // let's start with the force block, remember Bf has blocks that are a*I_3x3 per element
        q_ce[0] = elemArea[it->first] * v[0];
        q_ce[1] = elemArea[it->first] * v[1];
        q_ce[2] = elemArea[it->first] * v[2];

        // now we'll add the contribution of the torque block a*M
        // with M=[0 mz -my ; -mz 0 mx ; my -mx 0], mx,my,mz being the components of elemCtr[i]
        // this is skew-symmetric, so we just need to swap the signs to get the transpose
        q_ce[0] += elemArea[it->first] * (-elemCtr[it->first][2] * v[4] + elemCtr[it->first][1] * v[5]);
        q_ce[1] += elemArea[it->first] * (elemCtr[it->first][2] * v[3] - elemCtr[it->first][0] * v[5]);
        q_ce[2] += elemArea[it->first] * (-elemCtr[it->first][1] * v[3] + elemCtr[it->first][0] * v[4]);

        it->second -= q_ce; // apply the correction
    }

    //// testing ...
    //printf("%% ... original force %.3g (%.3lg %.3lg %.3lg)\n%% ... original torque %.3g (%.3lg %.3lg %.3lg)\n",
    //	force.norm(),force[0],force[1],force[2], torque.norm(),torque[0],torque[1],torque[2]);
    //force.setZero(); torque.setZero();
    //for(std::map<int, Eigen::Vector3f>::iterator it = q.begin(); it != q.end(); ++it){
    //	force  += it->second*elemArea[it->first];
    //	torque += it->second.cross(elemCtr[it->first])*elemArea[it->first];
    //}
    //printf("%% ... remaining force %.3g (%.3lg %.3lg %.3lg)\n%% ... remaining torque %.3g (%.3lg %.3lg %.3lg)\n",
    //	force.norm(),force[0],force[1],force[2], torque.norm(),torque[0],torque[1],torque[2]);

    // finally distribute remaining forces from collision elements, since these will receive Diri BCs
    // we're not doing this anymore and solve the Neumann problem with a regularizer instead
    //if(0) for(std::map<int, Eigen::Vector3f>::iterator it_c = contactTractions.begin(); it_c != contactTractions.end(); ++it_c){
    //	//printf("%% ... ... remaining traction on collision (%.3lg %.3lg %.3lg)\n",q[it->first][0],q[it->first][1],q[it->first][2]);
    //	force = -q[it_c->first] / (bem->getElems().size()-contactTractions.size()); // use force vector as buffer -- it's still a traction
    //	for(std::map<int, Eigen::Vector3f>::iterator it = q.begin(); it != q.end(); ++it){
    //		if( contactTractions.count(it->first)==0 ){ // this is not a collision element
    //			it->second += force;
    //		}
    //	}
    //	q[it_c->first].setZero();
    //}

    //
    // remove zero-value entries
    //
    std::map<int, Eigen::Vector3f> cpy;
    for (std::map<int, Eigen::Vector3f>::const_iterator i = balancedTractions.cbegin(); i != balancedTractions.cend(); ++i)
    {
        if (i->second.norm() > 0)
        {
            cpy.emplace(i->first, i->second);
        }
    }

    balancedTractions = cpy;
}

// volume of triangle mesh + its volume center of mass
btVector3 Object::getMeshVolumeAndCentreOfMass(const trimesh::TriMesh *mesh, float &meshVolume)
{
    btVector3 fragColShapeMeshCOM(0, 0, 0);
    meshVolume = 0;

    for (int i = 0; i < mesh->faces.size(); ++i)
    {
        const trimesh::TriMesh::Face &face = mesh->faces.at(i);
        const btVector3 &a = toBullet(mesh->vertices.at(face[0]));
        const btVector3 &b = toBullet(mesh->vertices.at(face[1]));
        const btVector3 &c = toBullet(mesh->vertices.at(face[2]));

        float vol = a.dot(b.cross(c)) / 6.0;
        meshVolume += vol;
        fragColShapeMeshCOM += vol * (a + b + c) / 4.0; // volume weighted centroid of tet (a,b,c,0) (simplex integration)
    }

    fragColShapeMeshCOM /= meshVolume;

    return fragColShapeMeshCOM;
}


float Object::getCenteredMeshFromVBD(
    trimesh::TriMesh *&pMesh,
    Eigen::Vector3f &centreOfMass,
    const openvdb::FloatGrid::Ptr& gridPtr,
    const float vdbVolToMeshAdaptivity)
{
    //openvdb::tools::VolumeToMesh volumeToMeshHandle(*gridPtr, 0, 0.5);
    //volumeToMeshHandle(*gridPtr);
    std::vector<openvdb::Vec3s> out_points;
    std::vector<openvdb::Vec4I > out_quads;
    std::vector<openvdb::Vec3I > out_tris;
    openvdb::FloatGrid::Ptr g = gridPtr->deepCopy();
    openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*g, out_points, out_tris, out_quads, 0, vdbVolToMeshAdaptivity);


    pMesh = new trimesh::TriMesh;

    printf("openvdb::tools::volumeToMesh: out_points = %d out_tris = %d out_quads = %d\n", (int)out_points.size(), (int)out_tris.size(), (int)out_quads.size());
    for(unsigned int i = 0; i < out_points.size(); i++) {
               //pMesh->vertices[i * 3] = out_points[i].x();
              //  pMesh->vertices[i * 3 + 1] = out_points[i].y();
              // pMesh->vertices[i * 3 + 2] = out_points[i].z();
               openvdb::Vec3s &v = (out_points)[i];
        // NOTE: This is a hack! The y coord is negated because of mixup with blender's XZY coordinate system ****************************************************************************************
        pMesh->vertices.push_back(trimesh::vec3(v[0], v[1], v[2]));
       }

       for(unsigned int i = 0; i < out_quads.size(); i++) {
       //        mesh->quads[i * 4] = out_quads[i].x();
      //         mesh->quads[i * 4 + 1] = out_quads[i].y();
      //         mesh->quads[i * 4 + 2] = out_quads[i].z();
       //        mesh->quads[i * 4 + 3] = out_quads[i].w();

                openvdb::Vec4I p = out_quads[i];

                trimesh::TriMesh::Face f0;
                f0[0] = p.z();
                f0[1] = p.y();
                f0[2] = p.x();
                trimesh::TriMesh::Face f1;
                f1[0] = p.w();
                f1[1] = p.z();
                f1[2] = p.x();

                pMesh->faces.push_back(f0);
                pMesh->faces.push_back(f1);
       }

       for(unsigned int i = 0; i < out_tris.size(); i++) {
              // mesh->triangles[i * 3] = out_tris[i].x();
              // mesh->triangles[i * 3 + 1] = out_tris[i].y();
              // mesh->triangles[i * 3 + 2] = out_tris[i].z();

               trimesh::TriMesh::Face f0;
                f0[0] = out_tris[i].z();
                f0[1] = out_tris[i].y();
                f0[2] = out_tris[i].x();

                pMesh->faces.push_back(f0);

                
       }


#if 0
    openvdb::tools::PointList *verts = &volumeToMeshHandle.pointList();
    openvdb::tools::PolygonPoolList *polys = &volumeToMeshHandle.polygonPoolList();

    for (size_t i = 0; i < volumeToMeshHandle.pointListSize(); i++)
    {
        openvdb::Vec3s &v = (*verts)[i];
        // NOTE: This is a hack! The y coord is negated because of mixup with blender's XZY coordinate system ****************************************************************************************
        pMesh->vertices.push_back(trimesh::vec3(v[0], v[1], v[2]));
        //file << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    for (size_t i = 0; i < volumeToMeshHandle.polygonPoolListSize(); i++)
    {

        for (size_t ndx = 0; ndx < (*polys)[i].numTriangles(); ndx++)
        {
            openvdb::Vec3I *p = &((*polys)[i].triangle(ndx));
            //file << "f " << p->x() << " " << p->y() << " " << p->z() << std::endl;
        }

        for (size_t ndx = 0; ndx < (*polys)[i].numQuads(); ndx++)
        {
            openvdb::Vec4I *p = &((*polys)[i].quad(ndx));

            trimesh::TriMesh::Face f0;
            f0[0] = p->z();
            f0[1] = p->y();
            f0[2] = p->x();
            trimesh::TriMesh::Face f1;
            f1[0] = p->w();
            f1[1] = p->z();
            f1[2] = p->x();

            pMesh->faces.push_back(f0);
            pMesh->faces.push_back(f1);

            //file << "f " << p->z() + 1 << " " << p->y() + 1 << " " << p->x() + 1 << std::endl;
            //file << "f " << p->w() + 1 << " " << p->z() + 1 << " " << p->x() + 1 << std::endl;
        }
    }
#endif


    getLargestConnectedComponent(pMesh);

    float meshVolume;
    btVector3 com = getMeshVolumeAndCentreOfMass(pMesh, meshVolume);
    centreOfMass = toEigen(com);

    // recenter the object's mesh around the center of mass
    for (int i = 0; i < pMesh->vertices.size(); ++i)
    {
        pMesh->vertices[i][0] -= centreOfMass[0];
        pMesh->vertices[i][1] -= centreOfMass[1];
        pMesh->vertices[i][2] -= centreOfMass[2];
    }

    return meshVolume;
}

void toLibiglMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, trimesh::TriMesh &tm)
{
    V.resize(tm.vertices.size(), 3);
    for (int i = 0; i < tm.vertices.size(); ++i)
    {
        const trimesh::point &p = tm.vertices[i];
        V(i, 0) = p[0];
        V(i, 1) = p[1];
        V(i, 2) = p[2];
    }

    F.resize(tm.faces.size(), 3);
    for (int i = 0; i < tm.faces.size(); ++i)
    {
        const trimesh::TriMesh::Face &f = tm.faces[i];
        F(i, 0) = f[0];
        F(i, 1) = f[1];
        F(i, 2) = f[2];
    }
}

void toTrimesh(trimesh::TriMesh &tm, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    tm.vertices.clear();
    for (int i = 0; i < V.rows(); ++i)
    {
        trimesh::point3 p(V(i, 0), V(i, 1), V(i, 2));
        tm.vertices.push_back(p);
    }

    tm.faces.clear();
    for (int i = 0; i < F.rows(); ++i)
    {
        trimesh::TriMesh::Face f(F(i, 0), F(i, 1), F(i, 2));
        tm.faces.push_back(f);
    }
}

void reduceTriangles(trimesh::TriMesh *&mesh, int target)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    toLibiglMesh(V, F, *mesh);

    Eigen::MatrixXd U;
    Eigen::MatrixXi G;
    Eigen::VectorXi J;
    Eigen::VectorXi I;
    printf("remesh target: %d\nrun igl::decimate\n", target);
    //igl::qslim(V, F, (std::size_t)target, U, G, J, I);
    igl::decimate(V, F, (std::size_t)target, U, G, J, I);

    toTrimesh(*mesh, U, G);
}

// reduce number of triangles
// if highResMesh->faces.size() <= targetTris, then no-op
void Object::decimateMesh(
    trimesh::TriMesh *&mesh,
    const int targetTris)
{
    SCOPED_TIMER(__FUNCTION__);

    if (mesh->faces.size() <= targetTris)
    {
        printf("%%skip mesh simplification\n");
        return;
    }

    printf("%% do mesh simplification (faces=%d, target=%d)\n", (int)mesh->faces.size(), targetTris);
#ifdef MeshDecimation_h
    
    Simplify::vertices.resize(mesh->vertices.size());

    for (int i = 0; i < mesh->vertices.size(); ++i) {
        const trimesh::vec3& p = mesh->vertices[i];
        Simplify::Vertex v;
        v.p.x = p[0];
        v.p.y = p[1];
        v.p.z = p[2];
        Simplify::vertices[i] = v;
    }

    Simplify::triangles.resize(mesh->faces.size());

    for (int i = 0; i < mesh->faces.size(); ++i) {
        const trimesh::TriMesh::Face& f = mesh->faces[i];
        Simplify::Triangle t;
        t.v[0] = f[0];
        t.v[1] = f[1];
        t.v[2] = f[2];
        Simplify::triangles[i] = t;
    }

    printf("%% input: %d triangles %d vertices\n", (int)Simplify::triangles.size(), (int)Simplify::vertices.size());
    Simplify::simplify_mesh(targetTris, 16, true);

    // Simplification done;
    // Visualize the Result

    printf("%% output: %d triangles %d vertices\n", (int)Simplify::triangles.size(), (int)Simplify::vertices.size());

    ASSERT(Simplify::triangles.size() > 0);
    ASSERT(Simplify::vertices.size() > 0);

    // wipe data
    mesh->vertices.clear();
    mesh->faces.clear();

    // copy vertices
    mesh->vertices.resize(Simplify::vertices.size());

    for (int i = 0; i < Simplify::vertices.size(); ++i) {
        const Simplify::Vertex& v = Simplify::vertices[i];
        mesh->vertices[i] = trimesh::vec3(v.p.x, v.p.y, v.p.z);
    }

    mesh->faces.resize(Simplify::triangles.size());

    for (int i = 0; i < Simplify::triangles.size(); ++i) {
        const Simplify::Triangle& f = Simplify::triangles[i];
        mesh->faces[i] = trimesh::TriMesh::Face(f.v[0], f.v[1], f.v[2]);
    }

    Simplify::triangles.clear();
    Simplify::vertices.clear();
#else
    reduceTriangles(mesh, targetTris);
#endif
    printf("%% new number of faces=%d\n", (int)mesh->faces.size());
    printf("%% decimation complete\n");
}

void Object::setCollMeshFromRenderingMesh(
    ColliderData &coll,
    const trimesh::TriMesh *renderingMesh,
    std::string writeMeshFile,
    unsigned int targetTris,
    const float vdbVolToMeshAdaptivity)
{
#if 0
    trimesh::TriMesh* levelSetTrimesh = nullptr;
    Eigen::Vector3f com;
    double volume = getCenteredMeshFromVBD(levelSetTrimesh, com, grid, vdbVolToMeshAdaptivity);

    int numPoints = levelSetTrimesh->vertices.size();
    ASSERT(numPoints > 0);
    int numTriangles = levelSetTrimesh->faces.size();
    ASSERT(numTriangles > 0);
    ASSERT(levelSetTrimesh != nullptr);
#endif
    coll.trimeshLowRes = new trimesh::TriMesh; // NOTE: thsi memory must be freed
    coll.trimeshLowRes->vertices = renderingMesh->vertices;
    coll.trimeshLowRes->faces = renderingMesh->faces;

    if (renderingMesh->faces.size() > (int)(targetTris * 1.1))
    {
        decimateMesh(coll.trimeshLowRes, targetTris);
    }

    if (!writeMeshFile.empty())
    {
#if 1 //defined(DUMP_DEBUG_MESH_DATA_TO_FILE)
        //coll.trimeshLowRes->write(writeMeshFile.c_str());
        //vcg::tri::io::ExporterOBJ<VcgMesh> writerObj;
        //writerObj.Save(mesh,writeMeshFile.append(".obj").c_str(), vcg::tri::io::Mask::IOM_NONE);
#else
        printf("%% skip logging coll-mesh to file\n");
#endif
    }

    coll.iface = new btTriangleIndexVertexArray(
        (int)coll.trimeshLowRes->faces.size(),
        (int *)&coll.trimeshLowRes->faces[0],
        3 * sizeof(int),
        (int)coll.trimeshLowRes->vertices.size(),
        (btScalar *)&coll.trimeshLowRes->vertices[0],
        3 * sizeof(float));

    btGImpactMeshShape *gimesh = new btGImpactMeshShape(coll.iface);
    gimesh->updateBound();
    coll.shape = gimesh;

    //return com;
}

void Object::setCollMeshFromRenderingMesh(
    ColliderData &coll,
    const trimesh::TriMesh *renderingMesh,
    openvdb::FloatGrid::Ptr& renderingMeshGridPtr,
    std::string writeMeshFile)
{
    coll.trimeshLowRes = new trimesh::TriMesh; // NOTE: thsi memory must be freed
#if 1

    if(renderingMesh->faces.size() < 16384)
    {
        coll.trimeshLowRes->faces = renderingMesh->faces;
        coll.trimeshLowRes->vertices = renderingMesh->vertices;
    }
    else
    {

        //trimesh::TriMesh* levelSetTrimesh = nullptr;
        float v_;
        auto lowReslevelSetGridPtr = meshToVDBLevelSetGrid(renderingMesh, renderingMeshGridPtr->voxelSize()[0] *10, "-", v_);
        //auto g = openvdb::tools::levelSetRebuild(*gridPtr, 0.0); // make sure the narrow-band is properly resolved (also removes most interior cracks)
        //g->setName(name.str());
        //g->setGridClass(openvdb::GRID_LEVEL_SET);

        Eigen::Vector3f com;
        double volume = getCenteredMeshFromVBD(coll.trimeshLowRes, com, lowReslevelSetGridPtr, 0);

        if(coll.trimeshLowRes->faces.empty())
        {
            coll.trimeshLowRes->faces = renderingMesh->faces;
            coll.trimeshLowRes->vertices = renderingMesh->vertices;
        }
        //int numPoints = levelSetTrimesh->vertices.size();
        //ASSERT(numPoints > 0);
        //int numTriangles = levelSetTrimesh->faces.size();
        //ASSERT(numTriangles > 0);
        //ASSERT(levelSetTrimesh != nullptr);
    }
#endif
    
    //coll.trimeshLowRes->vertices = renderingMesh->vertices;
    //coll.trimeshLowRes->faces = renderingMesh->faces;

    ///if (renderingMesh->faces.size() > (int)(targetTris * 1.1))
    //{
   //     decimateMesh(coll.trimeshLowRes, targetTris);
    //}

    if (!writeMeshFile.empty())
    {
#if 0 //defined(DUMP_DEBUG_MESH_DATA_TO_FILE)
        coll.trimeshLowRes->write(writeMeshFile.c_str());
        //vcg::tri::io::ExporterOBJ<VcgMesh> writerObj;
        //writerObj.Save(mesh,writeMeshFile.append(".obj").c_str(), vcg::tri::io::Mask::IOM_NONE);
#else
        printf("%% skip logging coll-mesh to file\n");
#endif
    }

    coll.iface = new btTriangleIndexVertexArray(
        (int)coll.trimeshLowRes->faces.size(),
        (int *)&coll.trimeshLowRes->faces[0],
        3 * sizeof(int),
        (int)coll.trimeshLowRes->vertices.size(),
        (btScalar *)&coll.trimeshLowRes->vertices[0],
        3 * sizeof(float));

    btGImpactMeshShape *gimesh = new btGImpactMeshShape(coll.iface);
    gimesh->updateBound();
    coll.shape = gimesh;

   
}

btRigidBody *Object::createConvexRigidBody(
    openvdb::FloatGrid::Ptr& grid,
    Object &parent,
    std::string writeMeshFile,
    double gridVolume,
    int renderMeshRemeshTarget)
{
    const float vdbVolToMeshAdaptivity = 0.2;
    trimesh::TriMesh *mesh;
    Eigen::Vector3f com;
    double volume = getCenteredMeshFromVBD(mesh, com, grid, vdbVolToMeshAdaptivity);
    //decimateMesh(mesh, renderMeshRemeshTarget);

    std::string fpath = writeMeshFile.append(".obj").c_str();
    printf("save rendering mesh: %s\n", fpath.c_str());
    mesh->write(fpath);
    //vcg::tri::io::ExporterOBJ<VcgMesh> writerObj;
    //writerObj.Save(mesh, writeMeshFile.append(".obj").c_str(), vcg::tri::io::Mask::IOM_NONE);

    // this method is supposed to be used for small fragments
    // a convex hull collision shape should be accurate enough ...
    btConvexHullShape *shape = new btConvexHullShape();
    for (int i = 0; i < mesh->vertices.size(); ++i)
    {
        //shape->addPoint(btVector3(mesh.vert[i].P()[0], mesh.vert[i].P()[1], mesh.vert[i].P()[2]), false); // don't update AABB of shape when inserting points
        shape->addPoint(toBullet(mesh->vertices[i]), false); // don't update AABB of shape when inserting points
    }
    mesh->vertices.clear();
    mesh->faces.clear();
    shape->recalcLocalAabb();

    com -= parent.centerOfMass; // parent.getCOMshift();

    btRigidBody *rb = createFragmentRB(
        parent.rigidBodyPtr, // parent.getRB(),
        shape,
        btVector3(com[0], com[1], com[2]),
        (btScalar)((gridVolume > 0.0 ? gridVolume : volume) * m_params->density /*parent.getMaterial()->getDensity()*/));

    return rb;
}

btRigidBody *Object::createFragmentRB(btRigidBody *parent, btCollisionShape *shape, const btVector3 &btCOM, btScalar mass)
{
    btVector3 inertia;
    shape->calculateLocalInertia(mass, inertia);

    // compute the position of the fragment in world space
    // based on the shift in COM from the parent's material space to it's own
    // as well as the parent's rotation
    btVector3 position = parent->getWorldTransform().getOrigin();
    btQuaternion rotation = parent->getWorldTransform().getRotation();
    position += quatRotate(rotation, btCOM);

    btRigidBody::btRigidBodyConstructionInfo rbci(
        mass,
        new btDefaultMotionState(btTransform(rotation, position)),
        shape,
        inertia);
    printf("%% ... new fragment has inertia (%.1le,%.1le,%.1le), com (%.3lf %.3lf %.3lf)\n", inertia[0], inertia[1], inertia[2], btCOM[0], btCOM[1], btCOM[2]);
    btRigidBody *rb = new btRigidBody(rbci);
    rb->setRestitution(parent->getRestitution());
    rb->setFriction(parent->getFriction());
    rb->setLinearVelocity(parent->getVelocityInLocalPoint(btCOM));
    rb->setAngularVelocity(parent->getAngularVelocity());
    return rb;
}

#if 0
double FractureRB::updateMass() {
  double volume = openvdb::tools::levelSetVolume(*fractSim->getLevelSet().getObjectGrid());
  double mass = volume * mat->getDensity();
  btVector3 inertia;
  rb->getCollisionShape()->calculateLocalInertia(mass, inertia);
  rb->setMassProps(mass, inertia);
  return volume;
}
#endif






#ifdef USE_LEVEL_SET


// 1. segment the level-set, get grids
// 2. check size - can be almost original, large or small
// 3. if we need to update original, do that
// 4. mesh the grids, make RBs
// 5. make breakable RBs for large ones
int Object::splitFragments(
    const Params& simParams,
    const GenericMesh& crackSurfaceMesh,
    const std::string& outDir,
    std::vector<Object*>& largeFragments,
    std::vector<btRigidBody*>& smallFragments)
{
    SCOPED_TIMER(__FUNCTION__);

    ASSERT(crackSurfaceMesh.m.faces.size() > 0);
    printf("\n%% splitFragments\n");

    largeFragments.clear();
    smallFragments.clear();

    // convert fracture points to grid
    const std::string& parentName = *((std::string*)this->rigidBodyPtr->getUserPointer());
    printf("parentName: %s (the object we are breaking)\n", parentName.c_str());

    printf("convert crack-mesh to VDB level set\n");

    openvdb::FloatGrid::Ptr fractureSurfaceLevelSetGridPtr = crackSurfaceMeshToLevelSetGrid(
        crackSurfaceMesh,
        simParams.vdbVoxelSize,
        parentName);

    openvdb::FloatGrid::Ptr copyOfCrackGrid = fractureSurfaceLevelSetGridPtr->deepCopy();
    openvdb::FloatGrid::Ptr copyOfObjGrid = this->levelSetGridPtr->deepCopy();

    // Compute the difference (A / B) of the two level sets.
    openvdb::tools::csgDifference(*copyOfObjGrid, *copyOfCrackGrid);

    openvdb::FloatGrid::Ptr csgSubtractedObjGrid = copyOfObjGrid; // cutted piece
    csgSubtractedObjGrid->setName(csgSubtractedObjGrid->getName() + "Cutted");

#if defined(DUMP_DEBUG_MESH_DATA_TO_FILE)
    saveOpenVDBGrid(csgSubtractedObjGrid, outDir);
#endif

    // list of fragment pieces after cutting (as level set grids)
    std::vector<openvdb::FloatGrid::Ptr> fragmentLevelSetGridPtrList;

    openvdb::tools::segmentSDF(*csgSubtractedObjGrid, fragmentLevelSetGridPtrList);

#if defined(DUMP_DEBUG_MESH_DATA_TO_FILE)
    //openvdb::GridPtrVec fragmentLevelSetGridPtrVec;
    //for (auto& g : fragmentLevelSetGridPtrList) {
    //    fragmentLevelSetGridPtrVec.push_back(g); //copy
    //}
    saveOpenVDBGrids(fragmentLevelSetGridPtrVec, outDir, csgSubtractedObjGrid->getName() + "Frags");
#endif

    printf("\n%% have %d fragments", (int)fragmentLevelSetGridPtrList.size()); // a bit of debug output

    // if there's only one segment we don't need to do anything, if there are
    // no segments, something went wrong and we can't do anything ... should
    // probably print a warning or something
    if (fragmentLevelSetGridPtrList.empty()) {
        return -1;
    }

    printf("\n%% compute volume");
    bool keepParent = false;
    //double parentVolume = openvdb::tools::levelSetVolume(*fractSim->getLevelSet().getObjectGrid());
    double parentVolume = openvdb::tools::levelSetVolume(*this->levelSetGridPtr);
    printf("\n%% parent volume: %f", parentVolume);
    float originalVolume = this->levelSetVolumeOrig;
    unsigned int ignoreCount = 0;

    if (parentVolume < 0.0) {
        printf("\n%% !!! negative parent volume !!!");
        parentVolume *= -1.0;
    }

    // for each fragment level
    for (unsigned int i = 0; i < fragmentLevelSetGridPtrList.size(); ++i) {
        printf("\n%% ------ process fragment %d", i);
        // get pointer/iterator to the current fragment
        openvdb::FloatGrid::Ptr fragmentLevelSetGridPtr = fragmentLevelSetGridPtrList[i];

        // calculate volume of fragment
        double segmentVolume = openvdb::tools::levelSetVolume(*fragmentLevelSetGridPtr);

        if (segmentVolume < 0.0) {
            printf("\n%% !!! negative segment volume (%.2lg) !!!", segmentVolume);
            segmentVolume *= -1.0; // something odd happened
        }

        printf("\n%% fragment volume %f", segmentVolume);

        // init the to-be name of the fragment level set

        printf(" vol ratio %f\n", segmentVolume / originalVolume);

        if (segmentVolume <= m_params->ignoredFragmentVolumeRatio * originalVolume) {
            printf("\nFRAGMENT IGNORED: too small");
            // do nothing (fragment too small)
            ++ignoreCount;
        } else if (segmentVolume <= m_params->smallFragmentVolumeRatio * originalVolume) { //parentVolume ){
            printf("\nSMALL FRAGMENT");
            std::stringstream name;
            name << parentName; // part of teh name is same as parent
            // make a small fragment
            name << "f" << ++fragmentCounter;
            printf("\n%% name %s\n", name.str().c_str());
            // use segmentVolume for mass computation
            // we want all masses to be VDB based so we have consistent measurements
            btRigidBody* fragmentRB;

            // some small fragments can be so small that they cause problems for bullet when handling collisions.
            // (objects passing through walls etc.). So we treat them as convex shapes which are easier to handle.
            #if 0
            if (segmentVolume <= m_params->smallConvexFragmentVolumeRatio * originalVolume) {
                printf("\n%% create as tiny convex fragment\n");
                // NOTE: "createConvexRigidBody" also saves the rendering mesh of tiny fragment
                fragmentRB = createConvexRigidBody(
                    fragmentLevelSetGridPtr, //segments[i],
                    *this,
                    pystring::os::path::join({ outDir, name.str() }),
                    segmentVolume,
                    simParams.renderMeshRemeshTarget);

				g_rigidBodyVolumes[name.str()] = segmentVolume;
            } else 
            #endif
            {
                printf("\n%% create as regular small meshed fragment\n");

                // bullet data for handling collisions
                ColliderData collMesh;
                //segments[i] = openvdb::tools::levelSetRebuild(*segments[i], 0.0, fractSim->getLevelSet().getBandwidth()); // make sure the narrow-band is properly resolved (also removes most interior cracks)
                fragmentLevelSetGridPtr->setName(name.str() + "Grid");
                fragmentLevelSetGridPtr->setGridClass(openvdb::GRID_LEVEL_SET);

                // save rendering mesh of small fragment
                trimesh::TriMesh* newRenderingMesh = nullptr;
                Eigen::Vector3f com;
                getCenteredMeshFromVBD(newRenderingMesh, com, fragmentLevelSetGridPtr, vdbVolToMeshAdaptivityForRenderMesh);
                ASSERT(newRenderingMesh != nullptr);
                #if 0
                decimateMesh(newRenderingMesh, simParams.renderMeshRemeshTarget * std::max(segmentVolume / originalVolume , 0.25));
                //renderingMesh->write(pystring::os::path::join({ outDir, name.str() + ".obj" }));
                #endif
                std::string fpath = pystring::os::path::join({ outDir, name.str() + ".obj" });
                printf("save rendering mesh: %s\n", fpath.c_str());
                newRenderingMesh->write(fpath);

                // here we now get the collision mesh from OpenVDB grid
                Eigen::Vector3f tmp;
                setCollMeshFromRenderingMesh(collMesh, newRenderingMesh, fragmentLevelSetGridPtr, pystring::os::path::join({ outDir, name.str() + "_colMesh.obj" })/*simParams.colMeshRemeshTarget*/); // limit the number of triangles to 20* the surface elements in the parent's BEM mesh
                // tmp = com - vdbCOM;
                tmp = com - this->centerOfMass; // "*this" is the parent
                fragmentRB = createFragmentRB(this->rigidBodyPtr, collMesh.shape, btVector3(tmp[0], tmp[1], tmp[2]), segmentVolume * m_params->density /*mat->getDensity()*/);
                storeData(collMesh);
                printf("%% meshed ");
				g_rigidBodyVolumes[name.str()] = segmentVolume;
            }

            if (fragmentRB) {
                fragmentRB->setUserPointer(new std::string(name.str()));
                fragmentRB->getCollisionShape()->setMargin(this->rigidBodyPtr->getCollisionShape()->getMargin());
                smallFragments.push_back(fragmentRB);
                printf("\n%% !!! small fragment has rest.coeff. %.3lf (%s)\n", fragmentRB->getRestitution(), name.str().c_str());
                printf("%% small fragment margin is %.3lg\n", fragmentRB->getCollisionShape()->getMargin());
            }

        } else if (segmentVolume >= m_params->originalFragmentVolumeRatio * parentVolume) { // originalVolume){// replace the implicit surface of the parent
            printf("\nUPDATE PARENT");
            std::stringstream name;
            name << parentName; // part of teh name is same as parent
            name << "" << ++(this->updateCounter);
            printf("\n%% name %s\n", name.str().c_str());
            keepParent = true;
            // update to reflect level-set shape update
            this->levelSetGridPtr = openvdb::tools::levelSetRebuild(*fragmentLevelSetGridPtr, 0.0); // make sure the narrow-band is properly resolved (also removes most interior cracks)
            this->levelSetGridPtr->setName(name.str());
            this->levelSetGridPtr->setGridClass(openvdb::GRID_LEVEL_SET);
            //fractSim->getLevelSet().setObjectGrid(segments[i]);

            // save updated rendering mesh
            delete this->renderingMesh; // delete old copy
            this->renderingMesh = nullptr;
            //updateCOM = this->centerOfMass; // assume that the center of mass does not change so much to prevent jitter in rendered mesh in blender
            //Eigen::Vector3f stub;
            getCenteredMeshFromVBD(this->renderingMesh, updateCOM, levelSetGridPtr, vdbVolToMeshAdaptivityForRenderMesh);
            //updateCOM = stub;
            ASSERT(this->renderingMesh != nullptr);
            //decimateMesh(renderingMesh, simParams.renderMeshRemeshTarget * std::max(segmentVolume / originalVolume , 0.25));
            //this->renderingMesh->write();

            std::string fpath = pystring::os::path::join({ outDir, name.str() + ".obj" });
            printf("save rendering mesh: %s\n", fpath.c_str());
            this->renderingMesh->write(fpath);

            if(this->iglRenderingMesh != nullptr)
            {
                delete this->iglRenderingMesh;
            }
            this->iglRenderingMesh = new iglTriMesh(fpath);

            printf("\n%% collision shape update on parent");

            //queue the update of the collision shape to a mesh of the remaining part ...
            if (pendingCollisionUpdate()) {
                this->collUpdate.deleteAll(); // we may have skipped an older update
            }

            setCollMeshFromRenderingMesh(this->collUpdate, renderingMesh, fragmentLevelSetGridPtr, pystring::os::path::join({ outDir, name.str() + "_colMesh.obj" }) /*simParams.colMeshRemeshTarget*/);
            updateVolume = segmentVolume;
            this->collUpdate.shape->setMargin(this->rigidBodyPtr->getCollisionShape()->getMargin());
            printf("\n%% collision shape update pending ... ");

            this->levelSetVolume = segmentVolume;
            //this->levelSetVolumeOrig = this->levelSetVolumeOrig;
            //this->centerOfMass = centerOfMass;
            //this->m_params = this->m_params;
            this->mass = segmentVolume * m_params->density;
            //this->fragmentCounter = this->fragmentCounter;
            //this->updateCounter = 0;

            this->initVolumeMeshData();
			g_rigidBodyVolumes[name.str()] = segmentVolume;
        } else { // make a large fragment
            printf("\nLARGE FRAGMENT");
            // ... get collision mesh
            ColliderData collMesh;
            std::stringstream name;
            name << parentName; // part of teh name is same as parent
            name << "F" << ++fragmentCounter;
            printf("\n%% name %s\n", name.str().c_str());

            //segments[i] = openvdb::tools::levelSetRebuild(*segments[i], 0.0, fractSim->getLevelSet().getBandwidth()); // make sure the narrow-band is properly resolved (also removes most interior cracks)
            fragmentLevelSetGridPtr->setName(name.str());
            fragmentLevelSetGridPtr->setGridClass(openvdb::GRID_LEVEL_SET);

            // update & save rendering mesh
            Eigen::Vector3f com;
            trimesh::TriMesh* newRenderingMesh = nullptr;
            getCenteredMeshFromVBD(newRenderingMesh, com, fragmentLevelSetGridPtr, vdbVolToMeshAdaptivityForRenderMesh);
            ASSERT(newRenderingMesh != nullptr);
            //decimateMesh(newRenderingMesh, simParams.renderMeshRemeshTarget * std::max(segmentVolume / originalVolume , 0.25));
            std::string fpath = pystring::os::path::join({ outDir, name.str() + ".obj" });
            printf("save rendering mesh: %s\n", fpath.c_str());
            newRenderingMesh->write(fpath);

            Eigen::Vector3f tmp;
            setCollMeshFromRenderingMesh(collMesh, newRenderingMesh, fragmentLevelSetGridPtr, pystring::os::path::join({ outDir, name.str() + "_colMesh.obj" }));

            // ... make rigid body
            //tmp = com - vdbCOM;
            tmp = com - this->centerOfMass; // "*this" is the parent
            btRigidBody* fragmentRB = createFragmentRB(this->rigidBodyPtr, collMesh.shape, btVector3(tmp[0], tmp[1], tmp[2]), segmentVolume * m_params->density /*mat->getDensity()*/);

            ASSERT(fragmentRB);

            fragmentRB->setUserPointer(new std::string(name.str())); // need to set this before calling initFractureSim on the FractureRB object
            fragmentRB->getCollisionShape()->setMargin(this->rigidBodyPtr->getCollisionShape()->getMargin());
            // ... create breakable rigid body from VDB
            Object* breakableFragment = new Object;
            breakableFragment->rigidBodyPtr = fragmentRB;
            breakableFragment->collInUse = collMesh;
            breakableFragment->levelSetGridPtr = fragmentLevelSetGridPtr;
            breakableFragment->m_params = this->m_params;

            breakableFragment->renderingMesh = newRenderingMesh;

            //breakableFragment->centerOfMass = renderingMeshCOM;
            breakableFragment->levelSetVolume = segmentVolume;
            breakableFragment->levelSetVolumeOrig = this->levelSetVolumeOrig;
            breakableFragment->iglRenderingMesh = new iglTriMesh(fpath);
            breakableFragment->centerOfMass = com;
            breakableFragment->m_params = this->m_params;
            breakableFragment->mass = segmentVolume * m_params->density;
            breakableFragment->fragmentCounter = this->fragmentCounter;
            breakableFragment->updateCounter = 0;

            breakableFragment->initVolumeMeshData();

            largeFragments.push_back(breakableFragment);
            printf("%% large fragment margin is %.3lg\n", fragmentRB->getCollisionShape()->getMargin());

			g_rigidBodyVolumes[name.str()] = segmentVolume;
        }
    }
    printf("\n%% ignored %u tiny fragments\n", ignoreCount);
    return keepParent ? (this->updateCounter) : 0;
}

#endif


#ifdef USE_MCUT 





// set collision mesh
void Object::setCollMesh_toColliderData(
    ColliderData &coll,
    trimesh::TriMesh *collisionMesh_Tri,
    std::string writeMeshFile)
{

    coll.trimeshLowRes = new trimesh::TriMesh; // NOTE: thsi memory must be freed
    coll.trimeshLowRes->vertices = collisionMesh_Tri->vertices;
    coll.trimeshLowRes->faces = collisionMesh_Tri->faces;

    coll.trimeshLowRes->write(writeMeshFile.c_str());

    coll.iface = new btTriangleIndexVertexArray(
        (int)coll.trimeshLowRes->faces.size(),
        (int *)&coll.trimeshLowRes->faces[0],
        3 * sizeof(int),
        (int)coll.trimeshLowRes->vertices.size(),
        (btScalar *)&coll.trimeshLowRes->vertices[0],
        3 * sizeof(float));

    btGImpactMeshShape *gimesh = new btGImpactMeshShape(coll.iface);
    gimesh->updateBound();
    coll.shape = gimesh;

    //return com;
}

// btRigidBody* Object::createConvexRigidBody_MeshCutting(
//     trimesh::TriMesh* mesh,
//     Object& parent,
//     std::string writeMeshFile,
//     Eigen::Vector3f fragmentCOM,
//     double fragmentVolume)
// {
//     std::string fpath = writeMeshFile.append(".obj").c_str();
//     printf("save rendering mesh: %s\n", fpath.c_str());
//     mesh->write(fpath);

//     // this method is supposed to be used for small fragments
//     // a convex hull collision shape should be accurate enough ...
//     btConvexHullShape* shape = new btConvexHullShape();
//     for (int i = 0; i < mesh->vertices.size(); ++i) {
//         //shape->addPoint(btVector3(mesh.vert[i].P()[0], mesh.vert[i].P()[1], mesh.vert[i].P()[2]), false); // don't update AABB of shape when inserting points
//         shape->addPoint(toBullet(mesh->vertices[i]), false); // don't update AABB of shape when inserting points
//     }
//     mesh->vertices.clear();
//     mesh->faces.clear();
//     shape->recalcLocalAabb();

//     fragmentCOM -= parent.centerOfMass; // parent.getCOMshift();

//     btRigidBody* rb = createFragmentRB(
//         parent.rigidBodyPtr, // parent.getRB(),
//         shape,
//         btVector3(fragmentCOM[0], fragmentCOM[1], fragmentCOM[2]),
//         (btScalar)(fragmentVolume * m_params->density /*parent.getMaterial()->getDensity()*/));

//     return rb;

// }

// split fragments using fTetwild
int Object::splitFragments_MeshCutting(
    const Params &simParams,
    const std::string &outDir,
    std::vector<meshObjFormat> *fragmentsCollision,
    std::vector<meshObjFormat> *fragmentsCollision_NoTri,
    std::vector<meshObjFormat> *fragmentsRendering,
    std::vector<std::vector<Eigen::Vector3d>> perturbVerticesMagVec,
    std::vector<Object *> &largeFragments,
    std::vector<btRigidBody *> &smallFragments)
{

    SCOPED_TIMER(__FUNCTION__);

    largeFragments.clear();
    smallFragments.clear();

    // convert fracture points to grid
    const std::string &parentName = *((std::string *)this->rigidBodyPtr->getUserPointer());
    printf("parentName: %s (the object we are breaking)\n", parentName.c_str());

    printf("\n%% have %d fragments", (int)fragmentsCollision->size()); // a bit of debug output

    printf("\n%% compute volume");
    bool keepParent = false;
    ////////////////////////////////////
    ////////////////////////////////////
    ////////////////////////////////////
    ////////////////////////////////////
    ////////////////////////////////////
    //double parentVolume = this->fragmentVolume;
    double parentVolume = this->levelSetVolume;

    printf("\n%% parent volume: %f", parentVolume);
    float originalVolume = this->levelSetVolumeOrig;
    unsigned int ignoreCount = 0;

    std::cout << "parentVolume = " << parentVolume << std::endl;
    std::cout << "originalVolume = " << originalVolume << std::endl;

    // for each fragment level
    for (unsigned int i = 0; i < fragmentsCollision->size(); ++i)
    {
        printf("\n%% ------ process fragment %d", i);
        // get pointer/iterator to the current fragment
        //openvdb::FloatGrid::Ptr fragmentLevelSetGridPtr = fragmentLevelSetGridPtrList[i];

        std::vector<Eigen::Vector3d> perturbVerticesMag = perturbVerticesMagVec[i];

        trimesh::TriMesh* collisionMesh = new trimesh::TriMesh;
        meshObjToTriMesh((*fragmentsCollision)[i], collisionMesh);

        meshObjFormat *highResCuttingMesh = new meshObjFormat;
        highResCuttingMesh->vertices = (*fragmentsCollision_NoTri)[i].vertices;
        highResCuttingMesh->faces = (*fragmentsCollision_NoTri)[i].faces;
 

        trimesh::TriMesh* highResRenderingMesh = new trimesh::TriMesh;
        meshObjToTriMesh((*fragmentsRendering)[i], highResRenderingMesh);

        // calculate volume of fragment
        float segmentVolume;
        btVector3 ComFragment = Object::getMeshVolumeAndCentreOfMass(highResRenderingMesh, segmentVolume);

        highResRenderingMesh->write("highResRenderingMesh_b4_recenter" + std::to_string(i) + ".obj");
        collisionMesh->write("collisionMesh_b4_recenter" + std::to_string(i) + ".obj");

        getCenteredMesh_Meshcutting(highResRenderingMesh, ComFragment);
        getCenteredMesh_Meshcutting(collisionMesh, ComFragment);

        highResRenderingMesh->write("highResRenderingMesh_after_recenter" + std::to_string(i) + ".obj");
        collisionMesh->write("collisionMesh_after_recenter" + std::to_string(i) + ".obj");

        std::string renderingMeshFilePath;
        // init the to-be name of the fragment level set
        printf(" vol ratio %f\n", segmentVolume / originalVolume);

        if (segmentVolume <= m_params->ignoredFragmentVolumeRatio * originalVolume)
        {
            printf("\nFRAGMENT IGNORED: too small");
            // do nothing (fragment too small)
            ++ignoreCount;
        }
        else if (segmentVolume <= m_params->smallFragmentVolumeRatio * originalVolume)
        { //parentVolume ){
            printf("\nSMALL FRAGMENT");
            std::stringstream name;
            name << parentName; // part of teh name is same as parent
            // make a small fragment
            name << "f" << ++fragmentCounter;
            printf("\n%% name %s\n", name.str().c_str());
            // use segmentVolume for mass computation
            // we want all masses to be VDB based so we have consistent measurements
            btRigidBody *fragmentRB;

            // some small fragments can be so small that they cause problems for bullet when handling collisions.
            // (objects passing through walls etc.). So we treat them as convex shapes which are easier to handle.
            if (segmentVolume <= m_params->smallConvexFragmentVolumeRatio * originalVolume)
            {
                ///////////////
                // bullet data for handling collisions
                ColliderData collMesh;
                // save rendering mesh of small fragment
                Eigen::Vector3f com = {ComFragment[0], ComFragment[1], ComFragment[2]};
                
                //getCenteredMesh_Meshcutting(highResRenderingMesh, ComFragment);
                //getCenteredMesh_Meshcutting(collisionMesh, ComFragment);

                std::string fpath = pystring::os::path::join({outDir, name.str() + ".obj"});
                renderingMeshFilePath = fpath;
                printf("save rendering mesh: %s\n", fpath.c_str());
                highResRenderingMesh->write(fpath);

                float ratio = segmentVolume / originalVolume;
                //trimesh::TriMesh *newCollisionMesh = collisionMesh;
                decimateMesh(collisionMesh, simParams.renderMeshRemeshTarget * std::max((double)ratio, 0.25));


                // here we now get the collision mesh from OpenVDB grid
                setCollMesh_toColliderData(collMesh, collisionMesh, pystring::os::path::join({outDir, name.str() + "_colMesh.obj"}));
                Eigen::Vector3f tmp = com - this->centerOfMass; // "*this" is the parent
                fragmentRB = createFragmentRB(this->rigidBodyPtr, collMesh.shape, btVector3(tmp[0], tmp[1], tmp[2]), segmentVolume * m_params->density /*mat->getDensity()*/);
                storeData(collMesh);
                printf("%% meshed ");
                g_rigidBodyVolumes[name.str()] = segmentVolume;
            }
            else
            {
                printf("\n%% create as regular small meshed fragment\n");

                ///////////////
                // bullet data for handling collisions
                ColliderData collMesh;
                // save rendering mesh of small fragment
                Eigen::Vector3f com = {ComFragment[0], ComFragment[1], ComFragment[2]};
                //getCenteredMesh_Meshcutting(highResRenderingMesh, ComFragment);
                //getCenteredMesh_Meshcutting(collisionMesh, ComFragment);

                std::string fpath = pystring::os::path::join({outDir, name.str() + ".obj"});
                renderingMeshFilePath = fpath;
                printf("save rendering mesh: %s\n", fpath.c_str());
                highResRenderingMesh->write(fpath);

                float ratio = segmentVolume / originalVolume;
                //trimesh::TriMesh *newCollisionMesh = collisionMesh;
                decimateMesh(collisionMesh, simParams.renderMeshRemeshTarget * std::max((double)ratio, 0.25));

                // here we now get the collision mesh from OpenVDB grid
                setCollMesh_toColliderData(collMesh, collisionMesh, pystring::os::path::join({outDir, name.str() + "_colMesh.obj"}));
                Eigen::Vector3f tmp = com - this->centerOfMass; // "*this" is the parent
                fragmentRB = createFragmentRB(this->rigidBodyPtr, collMesh.shape, btVector3(tmp[0], tmp[1], tmp[2]), segmentVolume * m_params->density /*mat->getDensity()*/);
                storeData(collMesh);
                printf("%% meshed ");
                g_rigidBodyVolumes[name.str()] = segmentVolume;
              
            }

            if (fragmentRB)
            {
                fragmentRB->setUserPointer(new std::string(name.str()));
                fragmentRB->getCollisionShape()->setMargin(this->rigidBodyPtr->getCollisionShape()->getMargin());
                smallFragments.push_back(fragmentRB);
                printf("\n%% !!! small fragment has rest.coeff. %.3lf (%s)\n", fragmentRB->getRestitution(), name.str().c_str());
                printf("%% small fragment margin is %.3lg\n", fragmentRB->getCollisionShape()->getMargin());
            }
        }
        else if (segmentVolume >= m_params->originalFragmentVolumeRatio * parentVolume)
        { // originalVolume){// replace the implicit surface of the parent
            printf("\nUPDATE PARENT");
            std::stringstream name;
            name << parentName; // part of teh name is same as parent
            name << "" << ++(this->updateCounter);
            printf("\n%% name %s\n", name.str().c_str());
            keepParent = true;
            // save updated rendering mesh
            delete this->renderingMesh; // delete old copy
            this->renderingMesh = nullptr;
            this->renderingMesh = highResRenderingMesh;


            // save updated cutting mesh
            delete this->cuttingMesh; // delete old copy
            this->cuttingMesh = nullptr;
            this->cuttingMesh = collisionMesh;


            // save rendering mesh of small fragment
            Eigen::Vector3f com = {ComFragment[0], ComFragment[1], ComFragment[2]};
            //getCenteredMesh_Meshcutting(this->renderingMesh, ComFragment);
            //getCenteredMesh_Meshcutting(collisionMesh, ComFragment);
            std::string fpath = pystring::os::path::join({outDir, name.str() + ".obj"});
            renderingMeshFilePath = fpath;
            printf("save rendering mesh: %s\n", fpath.c_str());
            this->renderingMesh->write(fpath);

            if(this->iglRenderingMesh != nullptr)
            {
                delete this->iglRenderingMesh;
            }
            this->iglRenderingMesh = new iglTriMesh(fpath);

            printf("\n%% collision shape update on parent");

            //queue the update of the collision shape to a mesh of the remaining part ...
            if (pendingCollisionUpdate())
            {
                this->collUpdate.deleteAll(); // we may have skipped an older update
            }



            float ratio = segmentVolume / originalVolume;
            //trimesh::TriMesh *newCollisionMesh = collisionMesh;
            //decimateMesh(newCollisionMesh, simParams.renderMeshRemeshTarget * std::max((double)ratio, 0.25));
            setCollMeshFromRenderingMesh(this->collUpdate, collisionMesh, pystring::os::path::join({outDir, name.str() + "_colMesh.obj"}), simParams.colMeshRemeshTarget);


            // here we now get the collision mesh from OpenVDB grid
            //setCollMesh_toColliderData(this->collUpdate, collisionMesh, pystring::os::path::join({outDir, name.str() + "_colMesh.obj"}));
            updateVolume = segmentVolume;
            this->collUpdate.shape->setMargin(this->rigidBodyPtr->getCollisionShape()->getMargin());
            printf("\n%% collision shape update pending ... ");

            this->perturbVerticesMag = perturbVerticesMag;

            this->levelSetVolume = segmentVolume;
            this->mass = segmentVolume * m_params->density;
            this->initVolumeMeshData();
            g_rigidBodyVolumes[name.str()] = segmentVolume;





            //meshObjFormat *newcollisionMesh_NoTri = highResCuttingMesh;
            getCenteredMesh_Obj_Meshcutting(highResCuttingMesh, ComFragment);
            this->cuttingMesh_NoTri = highResCuttingMesh;



        }
        else
        { 
            
            
            // make a large fragment
            printf("\nLARGE FRAGMENT");
            // ... get collision mesh
            ColliderData collMesh;
            std::stringstream name;
            name << parentName; // part of teh name is same as parent
            name << "F" << ++fragmentCounter;
            printf("\n%% name %s\n", name.str().c_str());

            // update & save rendering mesh
            Eigen::Vector3f com = {ComFragment[0], ComFragment[1], ComFragment[2]};
            //trimesh::TriMesh *newRenderingMesh = highResRenderingMesh;
            //getCenteredMesh_Meshcutting(newRenderingMesh, ComFragment);

            trimesh::TriMesh* collisionMesh_NoDecimate = collisionMesh;

            std::string fpath = pystring::os::path::join({outDir, name.str() + ".obj"});
            renderingMeshFilePath = fpath;
            printf("save rendering mesh: %s\n", fpath.c_str());
            highResRenderingMesh->write(fpath);



            Eigen::Vector3f tmp;


            float ratio = segmentVolume / originalVolume;
            //trimesh::TriMesh *newCollisionMesh = collisionMesh;
            //getCenteredMesh_Meshcutting(newCollisionMesh, ComFragment);
            decimateMesh(collisionMesh, simParams.renderMeshRemeshTarget * std::max((double)ratio, 0.25));
            setCollMeshFromRenderingMesh(collMesh, collisionMesh, pystring::os::path::join({outDir, name.str() + "_colMesh.obj"}), simParams.colMeshRemeshTarget);


            //meshObjFormat *newcollisionMesh_NoTri = highResCuttingMesh;
            //getCenteredMesh_Obj_Meshcutting(highResCuttingMesh, ComFragment);





            // ... make rigid body
            //tmp = com - vdbCOM;
            tmp = com - this->centerOfMass; // "*this" is the parent
            btRigidBody *fragmentRB = createFragmentRB(this->rigidBodyPtr, collMesh.shape, btVector3(tmp[0], tmp[1], tmp[2]), segmentVolume * m_params->density /*mat->getDensity()*/);

            ASSERT(fragmentRB);

            fragmentRB->setUserPointer(new std::string(name.str())); // need to set this before calling initFractureSim on the FractureRB object
            fragmentRB->getCollisionShape()->setMargin(this->rigidBodyPtr->getCollisionShape()->getMargin());
            // ... create breakable rigid body from VDB
            Object *breakableFragment = new Object;
            breakableFragment->rigidBodyPtr = fragmentRB;
            breakableFragment->collInUse = collMesh;
            breakableFragment->m_params = this->m_params;
            breakableFragment->renderingMesh = highResRenderingMesh;

            breakableFragment->cuttingMesh = collisionMesh_NoDecimate; 
            breakableFragment->cuttingMesh_NoTri = highResCuttingMesh; 
            breakableFragment->perturbVerticesMag = perturbVerticesMag;

            breakableFragment->iglRenderingMesh = new iglTriMesh(renderingMeshFilePath);
            //breakableFragment->centerOfMass = renderingMeshCOM;
            breakableFragment->levelSetVolume = segmentVolume;
            breakableFragment->levelSetVolumeOrig = this->levelSetVolumeOrig;
            breakableFragment->centerOfMass = com;
            breakableFragment->m_params = this->m_params;
            breakableFragment->mass = segmentVolume * m_params->density;
            breakableFragment->fragmentCounter = this->fragmentCounter;
            breakableFragment->updateCounter = 0;

            breakableFragment->initVolumeMeshData();

            largeFragments.push_back(breakableFragment);
            printf("%% large fragment margin is %.3lg\n", fragmentRB->getCollisionShape()->getMargin());

            g_rigidBodyVolumes[name.str()] = segmentVolume;




            openvdb::FloatGrid::Ptr levelSetGridPtr = meshToVDBLevelSetGrid(collisionMesh, simParams.vdbVoxelSize, name.str(), levelSetVolume);
            breakableFragment->levelSetGridPtr = levelSetGridPtr;

        }
    }
    printf("\n%% ignored %u tiny fragments\n", ignoreCount);
    return keepParent ? (this->updateCounter) : 0;
}

#endif

#ifdef USE_LEVEL_SET_FULLYCUT


// 1. segment the level-set, get grids
// 2. check size - can be almost original, large or small
// 3. if we need to update original, do that
// 4. mesh the grids, make RBs
// 5. make breakable RBs for large ones
int Object::splitFragments_OpenVdbFullyCut(
    const parametersSim& simParams,
    const std::vector<GenericMesh>& crackSurfaceMesh,
    const std::string& outDir,
    std::vector<Object*>& largeFragments,
    std::vector<btRigidBody*>& smallFragments)
{
    SCOPED_TIMER(__FUNCTION__);
    
    printf("\n%% splitFragments ------------------------------------------------------------ \n");

    largeFragments.clear();
    smallFragments.clear();

    // convert fracture points to grid
    const std::string& parentName = *((std::string*)this->rigidBodyPtr->getUserPointer());
    printf("parentName: %s (the object we are breaking)\n", parentName.c_str());

    printf("convert crack-mesh to VDB level set\n");


    // list of fragment pieces after cutting (as level set grids)
    std::vector<openvdb::FloatGrid::Ptr> fragmentLevelSetGridPtrList;
    for(int i = 0; i < crackSurfaceMesh.size(); i++)
    {
        openvdb::FloatGrid::Ptr fractureSurfaceLevelSetGridPtr = crackSurfaceMeshToLevelSetGrid(
            crackSurfaceMesh[i],
            simParams.vdbVoxelSize,
            parentName);

        openvdb::FloatGrid::Ptr copyOfCrackGrid = fractureSurfaceLevelSetGridPtr->deepCopy();
        openvdb::FloatGrid::Ptr copyOfObjGrid = this->levelSetGridPtr->deepCopy();

        // Compute the difference (A / B) of the two level sets.
        openvdb::tools::csgIntersection(*copyOfObjGrid, *copyOfCrackGrid);

        openvdb::FloatGrid::Ptr csgSubtractedObjGrid = copyOfObjGrid; // cutted piece
        csgSubtractedObjGrid->setName(csgSubtractedObjGrid->getName() + "Cutted");

        openvdb::tools::segmentSDF(*csgSubtractedObjGrid, fragmentLevelSetGridPtrList);
    }

    printf("\n%% have %d fragments", (int)fragmentLevelSetGridPtrList.size()); // a bit of debug output

    // if there's only one segment we don't need to do anything, if there are
    // no segments, something went wrong and we can't do anything ... should
    // probably print a warning or something
    if (fragmentLevelSetGridPtrList.empty()) {
        return -1;
    }

    printf("\n%% compute volume");
    bool keepParent = false;
    //double parentVolume = openvdb::tools::levelSetVolume(*fractSim->getLevelSet().getObjectGrid());
    double parentVolume = openvdb::tools::levelSetVolume(*this->levelSetGridPtr);
    printf("\n%% parent volume: %f", parentVolume);
    float originalVolume = this->levelSetVolumeOrig;
    unsigned int ignoreCount = 0;

    if (parentVolume < 0.0)
    {
        printf("\n%% !!! negative parent volume !!!");
        parentVolume *= -1.0;
    }

    // for each fragment level
    for (unsigned int i = 0; i < fragmentLevelSetGridPtrList.size(); ++i)
    {
        printf("\n%% ------ process fragment %d", i);


        // get pointer/iterator to the current fragment
        openvdb::FloatGrid::Ptr fragmentLevelSetGridPtr = fragmentLevelSetGridPtrList[i];

        // calculate volume of fragment
        double segmentVolume = openvdb::tools::levelSetVolume(*fragmentLevelSetGridPtr);

        if (segmentVolume < 0.0)
        {
            printf("\n%% !!! negative segment volume (%.2lg) !!!", segmentVolume);
            segmentVolume *= -1.0; // something odd happened
        }

        printf("\n%% fragment volume %f", segmentVolume);

        // init the to-be name of the fragment level set

        printf(" vol ratio %f\n", segmentVolume / originalVolume);

        std::string renderingMeshFilePath;

        if (segmentVolume <= m_params->ignoredFragmentVolumeRatio * originalVolume)
        {
            printf("\nFRAGMENT IGNORED: too small");
            // do nothing (fragment too small)
            ++ignoreCount;
        }
        else if (segmentVolume <= m_params->smallFragmentVolumeRatio * originalVolume)
        { //parentVolume ){
            printf("\nSMALL FRAGMENT");
            std::stringstream name;

            
#ifndef TIMEASNAME
            name << parentName; // part of teh name is same as parent
            // make a small fragment
            name << "f" << ++fragmentCounter;
            printf("\n%% name %s\n", name.str().c_str());
#else
            for(std::string::size_type i = 0; i < parentName.size(); ++i)
            {
                name << parentName[i];
                if(parentName[i] == '_')
                {
                    break;
                }   
            }

            time_t rawtime;
            struct tm * timeinfo;
            char buffer[80];
            time (&rawtime);
            timeinfo = localtime(&rawtime);
            strftime(buffer,sizeof(buffer),"%m%d%H%M%S",timeinfo);
            std::string str(buffer);


            name << str; // part of teh name is same as parent
            // make a small fragment
            name << "f" << ++fragmentCounter;
            printf("\n%% name %s\n", name.str().c_str());

#endif

            // use segmentVolume for mass computation
            // we want all masses to be VDB based so we have consistent measurements
            btRigidBody *fragmentRB;

#if 0
            // some small fragments can be so small that they cause problems for bullet when handling collisions.
            // (objects passing through walls etc.). So we treat them as convex shapes which are easier to handle.
            if (segmentVolume <= m_params->smallConvexFragmentVolumeRatio * originalVolume)
            {
                printf("\n%% create as tiny fragment\n");
                ///////////////
                // bullet data for handling collisions
                ColliderData collMesh;
                // save rendering mesh of small fragment
                trimesh::TriMesh *newRenderingMesh = nullptr;
                Eigen::Vector3f com;
                getCenteredMeshFromVBD(newRenderingMesh, com, fragmentLevelSetGridPtr, vdbVolToMeshAdaptivityForRenderMesh);
                ASSERT(newRenderingMesh != nullptr);
                //decimateMesh(newRenderingMesh, simParams.renderMeshRemeshTarget * std::max(segmentVolume / originalVolume, 0.25));

                std::string fpath = pystring::os::path::join({outDir, name.str() + ".obj"});
                renderingMeshFilePath = fpath;
                printf("save rendering mesh: %s\n", fpath.c_str());

                newRenderingMesh->write(fpath);
 

                // here we now get the collision mesh from OpenVDB grid
                Eigen::Vector3f tmp;
                setCollMeshFromRenderingMesh(collMesh, newRenderingMesh, fragmentLevelSetGridPtr, pystring::os::path::join({outDir, name.str() + "_colMesh.obj"})); // limit the number of triangles to 20* the surface elements in the parent's BEM mesh    
                tmp = com - this->centerOfMass; // "*this" is the parent
                //tmp = com + comTemp - this->centerOfMass;
                fragmentRB = createFragmentRB(this->rigidBodyPtr, collMesh.shape, btVector3(tmp[0], tmp[1], tmp[2]), segmentVolume * m_params->density /*mat->getDensity()*/);
                storeData(collMesh);
                printf("%% meshed ");
                g_rigidBodyVolumes[name.str()] = segmentVolume;
            }
            else
#endif
            {
                printf("\n%% create as small fragment\n");

                // bullet data for handling collisions
                ColliderData collMesh;
                // save rendering mesh of small fragment
                trimesh::TriMesh *newRenderingMesh = nullptr;
                Eigen::Vector3f com;
                getCenteredMeshFromVBD(newRenderingMesh, com, fragmentLevelSetGridPtr, vdbVolToMeshAdaptivityForRenderMesh);
                ASSERT(newRenderingMesh != nullptr);

                std::string fpath = pystring::os::path::join({outDir, name.str() + ".obj"});
                renderingMeshFilePath = fpath;
                printf("save rendering mesh: %s\n", fpath.c_str());
                newRenderingMesh->write(fpath);

                // here we now get the collision mesh from OpenVDB grid
                Eigen::Vector3f tmp;
                setCollMeshFromRenderingMesh(collMesh, newRenderingMesh, fragmentLevelSetGridPtr, pystring::os::path::join({outDir, name.str() + "_colMesh.obj"})); // limit the number of triangles to 20* the surface elements in the parent's BEM mesh
                tmp = com - this->centerOfMass; // "*this" is the parent
                //tmp = com + comTemp - this->centerOfMass;
                fragmentRB = createFragmentRB(this->rigidBodyPtr, collMesh.shape, btVector3(tmp[0], tmp[1], tmp[2]), segmentVolume * m_params->density /*mat->getDensity()*/);
                storeData(collMesh);
                printf("%% meshed ");
                g_rigidBodyVolumes[name.str()] = segmentVolume;
            } 

            if (fragmentRB)
            {
                fragmentRB->setUserPointer(new std::string(name.str()));
                fragmentRB->getCollisionShape()->setMargin(this->rigidBodyPtr->getCollisionShape()->getMargin());
                smallFragments.push_back(fragmentRB);
                printf("\n%% !!! small fragment has rest.coeff. %.3lf (%s)\n", fragmentRB->getRestitution(), name.str().c_str());
                printf("%% small fragment margin is %.3lg\n", fragmentRB->getCollisionShape()->getMargin());
            }
        }
        else if (segmentVolume >= m_params->originalFragmentVolumeRatio * parentVolume)
        { // originalVolume){// replace the implicit surface of the parent
            printf("\nUPDATE PARENT");
            std::stringstream name;
            name << parentName; // part of teh name is same as parent
            name << "" << ++(this->updateCounter);
            printf("\n%% name %s\n", name.str().c_str());


            keepParent = true;
            // update to reflect level-set shape update
            this->levelSetGridPtr = openvdb::tools::levelSetRebuild(*fragmentLevelSetGridPtr, 0.0); // make sure the narrow-band is properly resolved (also removes most interior cracks)
            this->levelSetGridPtr->setName(name.str());
            this->levelSetGridPtr->setGridClass(openvdb::GRID_LEVEL_SET);

            // save updated rendering mesh
            delete this->renderingMesh; // delete old copy
            this->renderingMesh = nullptr;
            getCenteredMeshFromVBD(this->renderingMesh, updateCOM, levelSetGridPtr, vdbVolToMeshAdaptivityForRenderMesh);
            ASSERT(this->renderingMesh != nullptr);

            std::string fpath = pystring::os::path::join({outDir, name.str() + ".obj"});
            renderingMeshFilePath = fpath;
            printf("save rendering mesh: %s\n", fpath.c_str());
            this->renderingMesh->write(fpath);

            if(this->iglRenderingMesh != nullptr)
            {
                delete this->iglRenderingMesh;
            }
            this->iglRenderingMesh = new iglTriMesh(fpath);


            printf("\n%% collision shape update on parent");

            //queue the update of the collision shape to a mesh of the remaining part ...
            if (pendingCollisionUpdate())
            {
                this->collUpdate.deleteAll(); // we may have skipped an older update
            }

            setCollMeshFromRenderingMesh(this->collUpdate, renderingMesh, levelSetGridPtr, pystring::os::path::join({outDir, name.str() + "_colMesh.obj"}));
            updateVolume = segmentVolume;
            this->collUpdate.shape->setMargin(this->rigidBodyPtr->getCollisionShape()->getMargin());
            printf("\n%% collision shape update pending ... ");

            this->levelSetVolume = segmentVolume;
            this->mass = segmentVolume * m_params->density;
            this->initVolumeMeshData();
            g_rigidBodyVolumes[name.str()] = segmentVolume;
        }
        else
        { // make a large fragment
            printf("\nLARGE FRAGMENT");
            // ... get collision mesh
            ColliderData collMesh;
            std::stringstream name;


#ifndef TIMEASNAME
            name << parentName; // part of teh name is same as parent
            name << "F" << ++fragmentCounter;
            printf("\n%% name %s\n", name.str().c_str());
#else
            for(std::string::size_type i = 0; i < parentName.size(); ++i)
            {
                name << parentName[i];
                if(parentName[i] == '_')
                {
                    break;
                }  
            }


            time_t rawtime;
            struct tm * timeinfo;
            char buffer[80];
            time (&rawtime);
            timeinfo = localtime(&rawtime);
            strftime(buffer,sizeof(buffer),"%m%d%H%M%S",timeinfo);
            std::string str(buffer);


            name << str; // part of teh name is same as parent
            // make a small fragment
            name << "F" << ++fragmentCounter;
            printf("\n%% name %s\n", name.str().c_str());
#endif

            //segments[i] = openvdb::tools::levelSetRebuild(*segments[i], 0.0, fractSim->getLevelSet().getBandwidth()); // make sure the narrow-band is properly resolved (also removes most interior cracks)
            fragmentLevelSetGridPtr->setName(name.str());
            fragmentLevelSetGridPtr->setGridClass(openvdb::GRID_LEVEL_SET);

            // update & save rendering mesh
            Eigen::Vector3f com;
            trimesh::TriMesh *newRenderingMesh = nullptr;
            getCenteredMeshFromVBD(newRenderingMesh, com, fragmentLevelSetGridPtr, vdbVolToMeshAdaptivityForRenderMesh);
            ASSERT(newRenderingMesh != nullptr);
            //decimateMesh(newRenderingMesh, simParams.renderMeshRemeshTarget);
            std::string fpath = pystring::os::path::join({outDir, name.str() + ".obj"});
            renderingMeshFilePath = fpath;
            printf("save rendering mesh: %s\n", fpath.c_str());
            newRenderingMesh->write(fpath);

            Eigen::Vector3f tmp;
            setCollMeshFromRenderingMesh(collMesh, newRenderingMesh, fragmentLevelSetGridPtr, pystring::os::path::join({outDir, name.str() + "_colMesh.obj"}));

            // ... make rigid body
            //tmp = com - this->centerOfMass; // "*this" is the parent
            tmp = com  - this->centerOfMass;
            btRigidBody *fragmentRB = createFragmentRB(this->rigidBodyPtr, collMesh.shape, btVector3(tmp[0], tmp[1], tmp[2]), segmentVolume * m_params->density /*mat->getDensity()*/);

            ASSERT(fragmentRB);

            fragmentRB->setUserPointer(new std::string(name.str())); // need to set this before calling initFractureSim on the FractureRB object
            fragmentRB->getCollisionShape()->setMargin(this->rigidBodyPtr->getCollisionShape()->getMargin());
            // ... create breakable rigid body from VDB
            Object *breakableFragment = new Object;
            breakableFragment->rigidBodyPtr = fragmentRB;
            breakableFragment->collInUse = collMesh;
            breakableFragment->levelSetGridPtr = fragmentLevelSetGridPtr;
            breakableFragment->m_params = this->m_params;
            breakableFragment->renderingMesh = newRenderingMesh;
            breakableFragment->iglRenderingMesh = new iglTriMesh(renderingMeshFilePath);
            breakableFragment->levelSetVolume = segmentVolume;
            breakableFragment->levelSetVolumeOrig = this->levelSetVolumeOrig;
            breakableFragment->centerOfMass = com;
            breakableFragment->m_params = this->m_params;
            breakableFragment->mass = segmentVolume * m_params->density;
            breakableFragment->fragmentCounter = this->fragmentCounter;
            breakableFragment->updateCounter = 0;

            breakableFragment->initVolumeMeshData();

            largeFragments.push_back(breakableFragment);
            printf("%% large fragment margin is %.3lg\n", fragmentRB->getCollisionShape()->getMargin());

            g_rigidBodyVolumes[name.str()] = segmentVolume;
        }

    }
    printf("\n%% ignored %u tiny fragments\n", ignoreCount);
    return keepParent ? (this->updateCounter) : 0;
}


#endif

