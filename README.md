
# Setup

* `git clone --recursive`
* `cd extern/openvdb`
* `git checkout master`
* install dependencies(eigen3) on system
* 
* `mkdir build && cd build`
* call cmake
    - if windows
        * `cmake -DCMAKE_TOOLCHAIN_FILE=C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows -A x64 -DOPENVDB_BUILD_VDB_VIEW=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_NO_SYSTEM_FROM_IMPORTED:BOOL=TRUE -DCMAKE_BUILD_TYPE=Release ..`
          - NOTE: the cmake flags are needed because of OpenVDB
    - else linux
        * `cmake -DBULLET2_USE_TBB_MULTITHREADING=ON -DCMAKE_BUILD_TYPE=Release ..`
    
* build the code
    - If you are on Ubuntu 
        * run `make -j2 mpm-fracture` 
    - else (Windows)
        * use visual studio

# Explicit Dependancies (most are git submodules)

* See the "extern" directory for all dependencies
    - Some of these projects have "sub-dependencies" (e.g. OpenVDB)

# How to setup a demo in Blender

1. Start new project
2. Create scene objects (ground, objects etc.)
  - To add an object
    1) create/load it
    2) set position
    3) scale it
    4) rotate it
    5) set center-of-mass to volume center-of-mass
    6) apply scaling
    7) GOTO step 3. below
3. Set physic of each object in the scene 
  - collision object must be a "mesh"
  - non-moving objects must be set as "static" 
  - the breakable object must be set as "dynamic" 
  - Make sure to apply ONLY SCALING for each object.
4. Save the blender project into folder `demo-config/<DEMO_NAME>/` 
5. Export JSON data 
  - open `scripts/blenderSceneExport.py` via the built-in text-editer
  - run the script (to generate meta-data about the scene for our simulator)

Once you have complete the above steps, the directory `demo-config/<DEMO_NAME>/`
will have the following files/directories:
- `demo-config/<DEMO_NAME>/meshes/` : meshes exported from blender
- `demo-config/<DEMO_NAME>/output/` : stores output from the simulator
- `demo-config/<DEMO_NAME>/<DEMO_NAME>.blend` : blender project file
- `demo-config/<DEMO_NAME>/<DEMO_NAME>.json` : metadata exported from blender

# Running the simulator

* `build/mpm-fracture <DEMO_NAME>`

The application will use the `<DEMO_NAME>` argument to infer where to 
look for the inputs. These inputs are from the demo's directory and are
listed as follows:

* The metadata from the JSON file
* The meshes exported from blender

There two outputs of the program which are saved in the folder `demo-config/<DEMO_NAME>/output`. 
The first is a script called `<DEMO_NAME>-blender.py`, which is used to script blender keyframes.
The second is a set of meshes (a lot) which will be loaded into blender by the script.

# Version control

## The demo files you should commit to the git-repo

* `demo-config/<DEMO_NAME>/<DEMO_NAME>.blend`
* `demo-config/<DEMO_NAME>/<DEMO_NAME>.conf`

Do not commit anything in the `build/` directory!

## Storing the simulation output files

* The simulation data are stored in /build/simulation-data folder. It contains meshes of rigid objects and a python file which imports meshes into blender. 

## Understanding Bullet engine

* Forums: https://pybullet.org/Bullet/phpBB3/index.php
* https://andysomogyi.github.io/mechanica/bullet.html

## How to remove a git submodule

Source: https://coderwall.com/p/csriig/remove-a-git-submodule

1) Delete the relevant line from the .gitmodules file.
2) Delete the relevant section from .git/config.
3) Run git rm --cached pathtosubmodule (no trailing slash).
4) Commit the superproject.
5) Delete the now untracked submodule files.

## Change to Bullet3

I made a minor change to the Bullet3 in the file `BulletCollision/CollisionDispatch/btCollisionDispatcherMt.cpp`

There is a bug is the source code, where btPersistentManifold objects are not cleared before 
searching for new contacts inside stepSimulation. This only happens when running Bullet3
in multi-threaded mode.

BEFORE
m_batchUpdating = true;
btParallelFor(0, pairCount, m_grainSize, updater);
m_batchUpdating = false;

AFTER

m_batchUpdating = true;
// FLOYD UPDATE: manuallY clear "m_manifoldsPtr" because the author of bullet forgot to
for (int i = m_manifoldsPtr.size()-1; i >= 0; --i)
{
  m_manifoldsPtr[i]->clearManifold();
  m_manifoldsPtr.pop_back();
}
btParallelFor(0, pairCount, m_grainSize, updater);
m_batchUpdating = false;

# TODO list

* Collision bug [done]
* Segfault 
* Measure average time to compute MPM timestep 
* Measure average number of MPM timesteps per fracture sim (after collision event)
* Record total number of MPM timesteps per fracture sim (for full demo)