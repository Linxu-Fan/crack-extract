
# What is this?

* This respository contains the code for the paper "Simulating brittle fracture with material points". The code uses voronoi diagram to approximate a medial surface from damaged points cloud. Although the code is used to extract crack for brittle fracture, it could nevertheless be used to approximate the medial surface of any shape with a particle representation. The input points have their coordinates(x,y,z), and phase filed value(0-1) where 0 represents a healthy point and 1 is a fully damaged point. The output is a surface approximating the medial surface which could be: 1) the exact approximation, or 2) faces that only define a complete detached shape. Our code enables cutting with the surface mesh. The domain can be cut with level set(OPENVDB) or direct mesh cutting(MCUT). For more information, pelease refer to our paper.

# How to build?

* `git clone --recursive`
* `git checkout master`
* `mkdir build && cd build`
* call cmake
    - if windows
        * `cmake -DCMAKE_TOOLCHAIN_FILE=C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows -A x64 -DOPENVDB_BUILD_VDB_VIEW=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_NO_SYSTEM_FROM_IMPORTED:BOOL=TRUE -DCMAKE_BUILD_TYPE=Release ..`
          - NOTE: the cmake flags are needed because of OpenVDB
    - else linux
        * `cmake -DCMAKE_BUILD_TYPE=Release ..`
    
* build the code
    - If you are on Ubuntu 
        * run `make -j6 crackExtraction` 
    - else (Windows)
        * use visual studio

# How to use?
* Write a text file and put the file into the build directory. Then pass the text file as the argument to the executable. The text file has 6 lines. See the demo file in the folder "example".
    - 1) the folder contains the object and point clouds file.
    - 2) the point cloud file name. The file has four columns which are separated by space. The first three columns are the (x,y,z) coordinate and the last one is the phase value.
    - 3) the object file name. It should be in .obj format.
    - 4) the cutting method.
        * "OPENVDB_FULL": full cut with openVDB. No partial crack is enabled.
        * "OPENVDB_PARTIAL": partial cut with openVDB. Partial crack is enabled.
        * "MCUT": direct mesh cutting with MCUT.        
    - 5) the medial surface resolution which should be roughly twice of the average point space.
    - 6) the openVDB voxel size. The cutting gap is thinner with a small voxel size.

# Examples: 1) Mode 1 test; 2) Glass

* Mode 1 test 
!(./example/mode1/mode1.png)

* Glass
!(./example/glass/glass.png)

