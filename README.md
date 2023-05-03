# plasmatic

Multiphysics finite element simulation code mostly written with C++ using PETSc and Eigen for linear algebra. Plasmatic is capable of thermal simulation in both 2D and 3D. It reads input mesh geometry in the Gmsh file format and outputs results in the VTK file format for visualization in ParaView.

## Building

```bash
git clone git@github.com:jam4375/plasmatic-cpp.git
cd plasmatic-cpp
mkdir build
cd build
cmake ..
make -j <NUM_CORES>
```

## Building Documentation

```bash
make docs
```

## Running Tests

```bash
cd <REPO_ROOT>/build
ctest -j <NUM_CORES>
```
