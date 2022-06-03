# GeneratingFunctionForCones

This program produces the generating function for cones in polyhedra. Four cases are taken into account: cone(A), v + cone(A), Ax <= 0, Ax <= b. 
To build it, one must invoke 
```clang++ -O3 -std=c++20 -march=native main.cpp```
and have the prerequisites installed. This hasn't been tested on MSVC, but should work provided you configure the prerequisites by yourself.
**Prerequisites**
1. tl::expected (https://github.com/TartanLlama/expected)
2. Eigen (https://eigen.tuxfamily.org/index.php?title=Main_Page)
3. tsl::robin_map (https://github.com/Tessil/robin-map)
4. xoroshiro256 (included in the source code)
