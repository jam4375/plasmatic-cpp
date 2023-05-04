# plasmatic

Multiphysics finite element simulation code mostly written with C++ using PETSc and Eigen for linear algebra. Plasmatic is capable of thermal simulation in both 2D and 3D. It reads input mesh geometry in the Gmsh file format and outputs results in the VTK file format for visualization in ParaView.

## Building

Use the following commands to build the code: (replace `<NUM_CORES>` with the number of CPU cores on your machine)

```bash
git clone git@github.com:jam4375/plasmatic-cpp.git
cd plasmatic-cpp
mkdir build
cd build
cmake ..
make -j <NUM_CORES>
```

## Building & Viewing Documentation

```bash
make plasmatic_docs
open docs/index.html
```

## Running Tests

```bash
cd <REPO_ROOT>/build
ctest -j <NUM_CORES>
```
