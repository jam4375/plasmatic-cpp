# templo
Template for new C++ projects

## Building
```bash
git clone git@github.com:jam4375/templo.git
cd templo
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