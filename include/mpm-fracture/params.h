#pragma once

#include "mpm-fracture/utils.h"
#include <unordered_map>

#define STALE_FLOAT_PARAM (-1.0f)
#define STALE_FLOAT3_PARAM btVector3(-1e10, -1e10, -1e10)
#define STALE_INT_PARAM (-1)

struct BreakableObjectParams {
    float youngsModulus = STALE_FLOAT_PARAM;
    float poissonsRatio = STALE_FLOAT_PARAM;
    float density = STALE_FLOAT_PARAM;
    float referenceStress = STALE_FLOAT_PARAM; // a particle will be damaged when its stress is larger than this threshold
    float modeIFractureEnergy = STALE_FLOAT_PARAM; // mode 1 fracture energy Gf, refer to damage partitioning paper

    float forceScaleFactor = STALE_FLOAT_PARAM; // mode 1 fracture energy Gf, refer to damage partitioning paper

    float meshRefineSpacing = STALE_FLOAT_PARAM; // float value(0, 1].It is ratio between refined triangle edge length to the particle spacing
    int chevronMarkMaxLength = STALE_INT_PARAM; // the maximum length of chevron marks
    float influenceRadius = STALE_FLOAT_PARAM; // float value(1, ).It is ratio of the radius to meshRefineSpacing; The radius is an area where inside no vertices should be visited.
    float smoothRadius = STALE_FLOAT_PARAM; // float value(1, ).It is ratio of the radius to meshRefineSpacing; The larger it is, the wider the mount is
    float pertubationMagitude = STALE_FLOAT_PARAM; //float value(1, ).It is ratio of the radius to meshRefineSpacing; pertubation magnitude for each vertex
    float ratioOfDeepenMarks = STALE_FLOAT_PARAM; //float value(0, 1).ratio of chevron marks that should be deepened
    float deepenMagnitude = STALE_FLOAT_PARAM; //float value(0, ).magnitude of deepened chevron marks

    // if a fragment has more volume than "originalFragmentVolumeRatio" times the original volume we'll replace
    // the parent's implicit surface with the fragment's
    // while keeping the high-res mesh etc. intact
    float originalFragmentVolumeRatio = STALE_FLOAT_PARAM;
    // ignore fragments that have less than "ignoredFragmentVolumeRatio" times the original volume
    float ignoredFragmentVolumeRatio = STALE_FLOAT_PARAM;
    // fragments that have less volume than FRAGMENT_SMALL_THR times the original volume are treated as small fragments
    // small fragments will get a standard rigid body and will not be breakable anymore
    float smallFragmentVolumeRatio = STALE_FLOAT_PARAM;
    // small fragments get convex hull collision shapes if smaller than "smallConvexFragmentVolumeRatio", otherwise meshes are used
    // meshes are not working well for very small objects
    float smallConvexFragmentVolumeRatio = STALE_FLOAT_PARAM;

    float particleDamageUpdateThreshold = STALE_FLOAT_PARAM; // real number(0...X)
    // # max ratio between damaged an undamaged particles before stopping MPM fracture sim
    float fractureSimDamageRatioThreshold = STALE_FLOAT_PARAM; // real number(0...1)
    //  # number of material points per object
    int materialPointCount = STALE_INT_PARAM;
    //  # mpm simulation timestep
    float mpmTimestep = STALE_FLOAT_PARAM;
    //  # maximum number of simulation timestep
    int mpmTimeStepsMax = STALE_INT_PARAM;
};

struct KinematicObjectParams {
    btVector3 initialLinearVelocity = STALE_FLOAT3_PARAM;
    btVector3 initialAngularVelocity = STALE_FLOAT3_PARAM;
};

struct Params {

    Params(const std::string& demoSettingsFilePath)
    {
        std::printf("\nload simulation parameters\n");

        loadSimulationParams(demoSettingsFilePath);
    }

    float rigidBodyImpulseThreshold;
    float rigidBodyForceThreshold;
    int rigidBodyTimeSteps;
    float rigidBodyTimeStepSize;
    float vdbVoxelSize;
    int colMeshRemeshTarget;
    int renderMeshRemeshTarget;

    std::map<std::string, BreakableObjectParams> breakableObjectParams;
    std::map<std::string, KinematicObjectParams> kinematicObjectParams;
    std::vector<std::pair<std::string, std::string>> collisionFilters;

private:
    Eigen::Vector3f parseVector3(const std::string& s) const
    {
        ASSERT(!s.empty());
        Eigen::Vector3f vector(0, 0, 0);
        std::vector<std::string> l;
        pystring::split(s, l, ",");
        ASSERT_(l.size() == 3, "expected 3 components, got %zu\n", l.size());

        for (int j = 0; j < 3; ++j) {
            vector[j] = std::atof(l[j].c_str());
        }
        return vector;
    }

    int parseInt(const std::string& s)
    {
        ASSERT(!s.empty());
        return std::atoi(s.c_str());
    }

    float parseReal(const std::string& s)
    {
        ASSERT(!s.empty());
        return std::atof(s.c_str());
    }

    bool parseBool(const std::string& s)
    {
        ASSERT(!s.empty());
        ASSERT(s == "true" || s == "false");
        return s == "true" ? true : false;
    }

    std::vector<std::string> parseStringList(const std::string& s)
    {
        ASSERT(!s.empty());
        std::vector<std::string> result;
        pystring::split(s, result, ";");

        if (result.size() == 1) {
            result.clear();
            pystring::split(s, result, ","); // try to split on comma
        }

        return result;
    }

    template <typename T, typename U>
    void copyGenericParamToObject(T& t, const U u, const U staleValue)
    {
        if (t == staleValue) {
            t = u;
            ASSERT(t != staleValue); // i.e. generic parameter was not defined in config file
        }
    }

    void loadSimulationParams(
        const std::string& demoSettingsFilePath)
    {
        std::string settings = readTextFile(demoSettingsFilePath);
        std::istringstream stream;
        stream.str(settings);
        std::string line;

        // using std::vector to maintain order as given in text file, which is necessary due to order-dependent definitions
        std::vector<std::pair<std::string, std::string>> paramVals;

        while (std::getline(stream, line)) {
            bool isComment = pystring::startswith(line, "#");

            if (isComment || line.empty()) {
                continue;
            }

            std::vector<std::string> chunks;
            pystring::split(line, chunks, "=", 1); // split on '='
            std::string paramName = chunks[0];
            std::string paramValue = chunks[1];

            if (paramValue.find('=') != std::string::npos) {
                printf("error: '=' on RHS\nline: %s\n", line.c_str());
                std::exit(1);
            }

            if (paramValue.find('#') != std::string::npos) {
                paramValue = paramValue.substr(0, paramValue.find('#')); // remove trailing comment
            }

            paramValue = pystring::rstrip(paramValue); // remove trailing whitespace

            ASSERT(!paramValue.empty());

            if (paramValue == "undefined") {
                continue; // skip undefined param
            }

            paramVals.push_back(std::make_pair(paramName, paramValue));
            //paramVals[paramName] = paramValue;
        }

        ASSERT(!paramVals.empty());

        for (std::vector<std::pair<std::string, std::string>>::const_iterator iter = paramVals.cbegin();
             iter != paramVals.cend(); ++iter) {

            const std::string& paramName = iter->first;
            const std::string& paramValue = iter->second;

            printf("param=`%s`, value=`%s`\n", paramName.c_str(), paramValue.c_str());

            if (paramName == "rigidBodyTimeSteps") {
                rigidBodyTimeSteps = parseInt(paramValue);
                ASSERT(rigidBodyTimeSteps > 0);
            } else if (paramName == "rigidBodyTimeStepSize") {
                rigidBodyTimeStepSize = parseReal(paramValue);
                ASSERT(rigidBodyTimeStepSize > 0 && rigidBodyTimeStepSize < 0.016);
            } else if (paramName == "rigidBodyImpulseThreshold") {
                rigidBodyImpulseThreshold = parseReal(paramValue);
                ASSERT(rigidBodyImpulseThreshold > 0);
            } else if (paramName == "rigidBodyForceThreshold") {
                rigidBodyForceThreshold = parseReal(paramValue);
                ASSERT(rigidBodyForceThreshold > 0);
            } else if (paramName == "vdbVoxelSize") {
                vdbVoxelSize = parseReal(paramValue);
                ASSERT(vdbVoxelSize > 0 && vdbVoxelSize <= 1.0);
            } else if (paramName == "colMeshRemeshTarget") {
                colMeshRemeshTarget = parseInt(paramValue);
                ASSERT(colMeshRemeshTarget > 0);
            } else if (paramName == "renderMeshRemeshTarget") {
                renderMeshRemeshTarget = parseInt(paramValue);
                ASSERT(renderMeshRemeshTarget > (2 ^ 15) && renderMeshRemeshTarget < (2 << 20));
            } else if (paramName == "breakableObjects") {
                std::vector<std::string> breakableObjects = parseStringList(paramValue);
                for (auto& i : breakableObjects) {
                    const std::string& breakableObjectName = i;
                    if (breakableObjectName == "undefined") {
                        continue;
                    }
                    breakableObjectParams.emplace(breakableObjectName, BreakableObjectParams()).second;
                }
            } else if (paramName == "kinematicObjects") {
                std::vector<std::string> kinematicObjects = parseStringList(paramValue);
                for (auto& i : kinematicObjects) {
                    const std::string& kinematicObjectName = i;

                    if (kinematicObjectName == "undefined") {
                        continue;
                    }
                    ASSERT(kinematicObjectName.find("_k") != std::string::npos); // kinematic objs are suffixed with "_k"
                    kinematicObjectParams.emplace(kinematicObjectName, KinematicObjectParams()).second;
                }
            }
            else if (paramName == "collisionFilters")
            {
              std::vector<std::string> collisionFiltersVec = parseStringList(paramValue);
              for (auto& i : collisionFiltersVec) {
                const std::string& filter = i;

                if (filter == "undefined") {
                  continue;
                }
               
                std::vector<std::string> l;
                pystring::split(filter, l, "&");
                ASSERT_(l.size() == 2, "expected 2 values, got %zu\n", l.size());
                collisionFilters.push_back(std::pair<std::string, std::string>(l.front(), l.back()));
              }
              
            }

            //
            // object (kinematic | breakable) specific parameters
            //

            int dotPos = pystring::index(paramName, ".");

            if (dotPos != -1) {

                // get object name as prefix

                std::string objectName = paramName.substr(0, dotPos);
                std::map<std::string, BreakableObjectParams>::iterator biter = breakableObjectParams.find(objectName);
                std::map<std::string, KinematicObjectParams>::iterator kiter = kinematicObjectParams.find(objectName);

                bool objIsBreakable = (biter != breakableObjectParams.end());
                bool objIsKinematic = (kiter != kinematicObjectParams.end());

                //ASSERT( !(objIsBreakable && objIsKinematic)); // cannot be both!

                if (!objIsBreakable && !objIsKinematic) { // wildcard

                    if (objectName == "*") {
                        biter = breakableObjectParams.emplace(objectName, BreakableObjectParams()).first;
                    }

                    if (objectName == "*_k") {
                        kiter = kinematicObjectParams.emplace(objectName, KinematicObjectParams()).first;
                    }

                    if (kiter == kinematicObjectParams.end() && biter == breakableObjectParams.end()) {
                        continue;
                    }
                }

                BreakableObjectParams& breakableObjParams = biter->second;
                KinematicObjectParams& kinematicObjParams = kiter->second;

                std::string paramNameSuffix = paramName.substr(dotPos + 1);

                if (paramNameSuffix == "youngsModulus") {
                    breakableObjParams.youngsModulus = parseReal(paramValue);
                    ASSERT(breakableObjParams.youngsModulus > 0);
                } else if (paramNameSuffix == "poissonsRatio") {
                    breakableObjParams.poissonsRatio = parseReal(paramValue);
                    ASSERT(breakableObjParams.poissonsRatio > 0 && breakableObjParams.poissonsRatio < 0.5);
                } else if (paramNameSuffix == "density") {
                    breakableObjParams.density = parseReal(paramValue);
                    ASSERT(breakableObjParams.density > 0);
                } else if (paramNameSuffix == "referenceStress") {
                    breakableObjParams.referenceStress = parseReal(paramValue);
                    ASSERT(breakableObjParams.referenceStress > 0);
                } else if (paramNameSuffix == "modeIFractureEnergy") {
                    breakableObjParams.modeIFractureEnergy = parseReal(paramValue);
                    ASSERT(breakableObjParams.modeIFractureEnergy > 0);
                } else if (paramNameSuffix == "forceScaleFactor") {
                    breakableObjParams.forceScaleFactor = parseReal(paramValue);
                    ASSERT(breakableObjParams.forceScaleFactor > 0);
                } else if (paramNameSuffix == "meshRefineSpacing") {
                    breakableObjParams.meshRefineSpacing = parseReal(paramValue);
                    ASSERT(breakableObjParams.meshRefineSpacing > 0);
                } else if (paramNameSuffix == "chevronMarkMaxLength") {
                    breakableObjParams.chevronMarkMaxLength = parseReal(paramValue);
                    ASSERT(breakableObjParams.chevronMarkMaxLength > 0);
                } else if (paramNameSuffix == "influenceRadius") {
                    breakableObjParams.influenceRadius = parseReal(paramValue);
                    ASSERT(breakableObjParams.influenceRadius > 0);
                } else if (paramNameSuffix == "smoothRadius") {
                    breakableObjParams.smoothRadius = parseReal(paramValue);
                    ASSERT(breakableObjParams.smoothRadius > 0);
                } else if (paramNameSuffix == "pertubationMagitude") {
                    breakableObjParams.pertubationMagitude = parseReal(paramValue);
                    ASSERT(breakableObjParams.pertubationMagitude > 0);
                } else if (paramNameSuffix == "ratioOfDeepenMarks") {
                    breakableObjParams.ratioOfDeepenMarks = parseReal(paramValue);
                    ASSERT(breakableObjParams.ratioOfDeepenMarks > 0);
                } else if (paramNameSuffix == "deepenMagnitude") {
                    breakableObjParams.deepenMagnitude = parseReal(paramValue);
                    ASSERT(breakableObjParams.deepenMagnitude > 0);
                } else if (paramNameSuffix == "originalFragmentVolumeRatio") {
                    breakableObjParams.originalFragmentVolumeRatio = parseReal(paramValue);
                    ASSERT(breakableObjParams.originalFragmentVolumeRatio > 0);
                } else if (paramNameSuffix == "ignoredFragmentVolumeRatio") {
                    breakableObjParams.ignoredFragmentVolumeRatio = parseReal(paramValue);
                    ASSERT(breakableObjParams.ignoredFragmentVolumeRatio > 0);
                } else if (paramNameSuffix == "smallFragmentVolumeRatio") {
                    breakableObjParams.smallFragmentVolumeRatio = parseReal(paramValue);
                    ASSERT(breakableObjParams.smallFragmentVolumeRatio > 0);
                } else if (paramNameSuffix == "smallConvexFragmentVolumeRatio") {
                    breakableObjParams.smallConvexFragmentVolumeRatio = parseReal(paramValue);
                    ASSERT(breakableObjParams.smallConvexFragmentVolumeRatio > 0);
                } else if (paramNameSuffix == "particleDamageUpdateThreshold") {
                    breakableObjParams.particleDamageUpdateThreshold = parseReal(paramValue);
                    ASSERT(breakableObjParams.particleDamageUpdateThreshold > 0);
                } else if (paramNameSuffix == "fractureSimDamageRatioThreshold") {
                    breakableObjParams.fractureSimDamageRatioThreshold = parseReal(paramValue);
                    ASSERT(breakableObjParams.fractureSimDamageRatioThreshold > 0 && breakableObjParams.fractureSimDamageRatioThreshold < 1.0f);
                } else if (paramNameSuffix == "materialPointCount") {
                    breakableObjParams.materialPointCount = parseInt(paramValue);
                    ASSERT(breakableObjParams.materialPointCount > 0);
                } else if (paramNameSuffix == "mpmTimestep") {
                    breakableObjParams.mpmTimestep = parseReal(paramValue);
                } else if (paramNameSuffix == "mpmTimeStepsMax") {
                    breakableObjParams.mpmTimeStepsMax = parseInt(paramValue);
                    ASSERT(breakableObjParams.mpmTimeStepsMax > 0);
                }

                // kinematic objects
                else if (paramNameSuffix == "initialLinearVelocity") {
                    kinematicObjParams.initialLinearVelocity = toBullet(parseVector3(paramValue));
                } else if (paramNameSuffix == "initialAngularVelocity") {
                    kinematicObjParams.initialAngularVelocity = toBullet(parseVector3(paramValue));
                }
            }
        }

        //ASSERT(breakableObjectParams.size() >= 2); // ... including wildcard "*" parameters

        ASSERT(breakableObjectParams.size() >= 1); // ... wildcard "*" parameters

        // set any unspecified object-specific parameters from the wildcard parameters

        std::map<std::string, BreakableObjectParams>::const_iterator genericBreakableParamsIter = breakableObjectParams.find("*");

        for (std::map<std::string, BreakableObjectParams>::iterator iter = breakableObjectParams.begin();
             iter != breakableObjectParams.end(); ++iter) {

            if (iter->first == genericBreakableParamsIter->first) {
                continue; // we care only about copying from generic-to-object-specific
            }

            std::cout << "copy generic params for " << iter->first << std::endl;

            copyGenericParamToObject(iter->second.youngsModulus, genericBreakableParamsIter->second.youngsModulus, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.poissonsRatio, genericBreakableParamsIter->second.poissonsRatio, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.density, genericBreakableParamsIter->second.density, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.referenceStress, genericBreakableParamsIter->second.referenceStress, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.modeIFractureEnergy, genericBreakableParamsIter->second.modeIFractureEnergy, STALE_FLOAT_PARAM);

            copyGenericParamToObject(iter->second.forceScaleFactor, genericBreakableParamsIter->second.forceScaleFactor, STALE_FLOAT_PARAM);

            copyGenericParamToObject(iter->second.meshRefineSpacing, genericBreakableParamsIter->second.meshRefineSpacing, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.chevronMarkMaxLength, genericBreakableParamsIter->second.chevronMarkMaxLength, STALE_INT_PARAM);
            copyGenericParamToObject(iter->second.influenceRadius, genericBreakableParamsIter->second.influenceRadius, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.smoothRadius, genericBreakableParamsIter->second.smoothRadius, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.pertubationMagitude, genericBreakableParamsIter->second.pertubationMagitude, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.ratioOfDeepenMarks, genericBreakableParamsIter->second.ratioOfDeepenMarks, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.deepenMagnitude, genericBreakableParamsIter->second.deepenMagnitude, STALE_FLOAT_PARAM);

            copyGenericParamToObject(iter->second.originalFragmentVolumeRatio, genericBreakableParamsIter->second.originalFragmentVolumeRatio, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.ignoredFragmentVolumeRatio, genericBreakableParamsIter->second.ignoredFragmentVolumeRatio, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.smallFragmentVolumeRatio, genericBreakableParamsIter->second.smallFragmentVolumeRatio, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.smallConvexFragmentVolumeRatio, genericBreakableParamsIter->second.smallConvexFragmentVolumeRatio, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.particleDamageUpdateThreshold, genericBreakableParamsIter->second.particleDamageUpdateThreshold, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.fractureSimDamageRatioThreshold, genericBreakableParamsIter->second.fractureSimDamageRatioThreshold, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.materialPointCount, genericBreakableParamsIter->second.materialPointCount, STALE_INT_PARAM);
            copyGenericParamToObject(iter->second.mpmTimestep, genericBreakableParamsIter->second.mpmTimestep, STALE_FLOAT_PARAM);
            copyGenericParamToObject(iter->second.mpmTimeStepsMax, genericBreakableParamsIter->second.mpmTimeStepsMax, STALE_INT_PARAM);
        }

        std::map<std::string, KinematicObjectParams>::const_iterator genericKinematicParamsIter = kinematicObjectParams.find("*_k");

        for (std::map<std::string, KinematicObjectParams>::iterator iter = kinematicObjectParams.begin();
             iter != kinematicObjectParams.end(); ++iter) {

            if (iter->first == genericKinematicParamsIter->first) {
                continue; // we care only about copying from generic-to-object-specific
            }

            std::cout << "copy generic params for " << iter->first << std::endl;

            // kinematic parameters
            copyGenericParamToObject(iter->second.initialLinearVelocity, genericKinematicParamsIter->second.initialLinearVelocity, STALE_FLOAT3_PARAM);
            copyGenericParamToObject(iter->second.initialAngularVelocity, genericKinematicParamsIter->second.initialAngularVelocity, STALE_FLOAT3_PARAM);
        }
    }
};
