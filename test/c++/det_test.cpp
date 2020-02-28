#include <triqs/arrays/matrix.hpp>
#include <iostream>
using namespace triqs::arrays;

int main() { auto A = matrix<double>{{-0.0000000e+00, 1.9887142e-310}, {-1.5483636e-176, 1.8928403e+00}}; 
  std::cout << A<< "\n\n";
  std::cout << determinant(A) << "\n\n" ;
}
