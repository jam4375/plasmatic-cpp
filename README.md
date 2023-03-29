# plasmatic
Finite element simulation code

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