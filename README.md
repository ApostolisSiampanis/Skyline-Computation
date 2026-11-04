# Skyline Computation

## [University of Piraeus](https://www.unipi.gr/en/home/) | [Department of Informatics](https://cs.unipi.gr/en/)
**MSc course**: Algorithmic Techniques and Applications

**Semester**: 1

**Project Completion Year**: 2025

## Project Overview

This project computes the **skyline** of multi-dimensional datasets using the **BSkyTree** algorithm in C++.

Given a set of points, the skyline is the subset of points that are not dominated by any other point. A point *p* dominates a point *q* when *p* is at least as good as *q* in every dimension and strictly better in at least one. Skyline queries are used in multi-criteria decision making, for example to find the best trade-offs between price, quality and distance.

The implementation is based on [SkyBench](https://github.com/sean-chester/SkyBench) by Darius Šidlauskas, Sean Chester and Kenneth S. Bøgh. This repository keeps only the sequential BSkyTree algorithm: the other skyline algorithms and the OpenMP parallelisation were removed, and the benchmark driver was simplified accordingly.

## Key Features

- **BSkyTree algorithm:** Pivot selection, lattice-based partitioning and dominance tests to compute the skyline efficiently.
- **AVX dominance tests:** Dominance checks between points use AVX/SSE vector instructions.
- **Benchmark driver:** Reads a dataset, runs the algorithm and reports input reading time, initialisation time, runtime and skyline size.
- **Real-world workloads:** Three datasets are included for testing.
- **Debug and profiling options:** Verbose output, dominance test counters and a debug build that works with Valgrind.

## Datasets

The `workloads` folder contains three datasets. File names follow the pattern `<name>-U-<dimensions>-<points>.csv`, and all values are normalised to the range [0, 1].

| Dataset | Dimensions | Points |
|---|---|---|
| `nba-U-8-17264.csv` | 8 | 17,264 |
| `house-U-6-127931.csv` | 6 | 127,931 |
| `elv_weather-U-15-566268.csv` | 15 | 566,268 |

## Technologies Used

- C++ (C++11)
- GCC with AVX support (`-mavx`)
- GNU Make

## Project Structure

```
/src
    testdriver.cpp           # Benchmark driver (entry point)
    /bskytree                # BSkyTree algorithm, tree nodes and pivot selection
    /common                  # Shared definitions, skyline interface, AVX dominance tests
    /util                    # Timing and utility functions
/workloads                   # Input datasets (CSV)
makefile
```

## Setup Instructions

### Requirements

- Linux (the build links against `librt`)
- `g++` and `make`
- A CPU with AVX support

### Build

The number of dimensions is fixed at compile time, so build with the `DIMS` value of the dataset you want to run:

```bash
git clone https://github.com/ApostolisSiampanis/Skyline-Computation.git
cd Skyline-Computation
make DIMS=6
```

This creates the executable `bin/SkyBench`. Run `make clean` before building again with a different `DIMS` value.

Other build options:

- `make dbg`: Debug build without optimisations.
- `make DT=1`: Count the dominance tests per point.

### Run

```bash
./bin/SkyBench -f workloads/house-U-6-127931.csv -s "bskytree"
```

| Option | Description |
|---|---|
| `-f` | Input dataset (CSV) |
| `-s` | Algorithms to run (only `bskytree` is supported) |
| `-v` | Verbose mode (not for performance measurements) |

## Acknowledgments
This project was developed as part of the "Algorithmic Techniques and Applications" MSc course at the University of Piraeus. It is based on [SkyBench](https://github.com/sean-chester/SkyBench) by Darius Šidlauskas, Sean Chester and Kenneth S. Bøgh, released under the MIT License. Contributions and feedback are always welcome!

## License

This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.
