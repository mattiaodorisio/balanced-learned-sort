# Balanced Learned Sort

This is the implementation of the Balanced Learned Sort, presented in the paper [Balanced Learned Sort: a new learned model for fast and balanced item bucketing](link).

## Dependencies
`cmake` to build the project, `python` and `numpy` to download/prepare the real-world datasets to reproduce the experiments.

## Usage

To use BLS in your project: add it as submodule,
copy the `include` directory to your system's or project's include path.
And use it (have a look also at `src/example.cpp` for a working example):

```C++
#include "bls.h"

std::vector<int> v = {5, 1, 2, 4, 3};
ls_framework::bls(v.begin(), v.end());
```

The reproduce the experiments, clone this repository and its dependencies:

```bash
git clone https://github.com/mattiaodorisio/balanced-learned-sort
cd balanced-learned-sort
git submodule update --init
mkdir build && cd build
cmake ..
make -j 8
cd ..
./scripts/run_synth.sh
```

The reproduce the experiments on real datasets you need also to download them:
```bash
./scripts/download.sh
./scripts/run_reals.sh
```
