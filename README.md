# plasmatic

Multiphysics finite element simulation code mostly written with C++ using PETSc and Eigen for linear algebra. Plasmatic is capable of thermal and mechanical simulation in 2D and 3D. It reads input mesh geometry in the Gmsh file format and outputs results in the VTK file format for visualization in ParaView.

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

## Running a Mechanical Simulation

Run the following command to start a mechanical simulation. The output will be written to the file `mechanical.vtk` in the same directory and can be viewed with ParaView.
```bash
cd build/bin
./plasmatic -i mechanical.json
```
where `mechanical.json` contains input configuration similar to below:
```json
{
  "command": "run_mechanical_sim",
  "mesh_filepath": "assets/ProblemTypes/mesh3d_quadratic.msh",
  "youngs_modulus": 69.0e9,
  "poisson_ratio": 0.32,
  "displacement_bcs": [{ "surface_name": "fixed", "value": [0.0, 0.0, 0.0] }],
  "traction_bcs": [{ "surface_name": "load", "value": [0.0, -100.0, 0.0] }],
  "output_file": "mechanical"
}
```
