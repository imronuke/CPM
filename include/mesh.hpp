#ifndef MESH_HPP
#define MESH_HPP

#include <vector>

class Mesh {
 public:
  static int n_ring;
  static double* radius;
  static double* volume;
  double* mesh;

  Mesh();   // Constructor declaration
  ~Mesh();  // Destructor declaration

  static void calculate_volume(std::vector<double>);
};

#endif