#include "mesh.hpp"

#include <cmath>
#include <vector>

int Mesh::n_ring = 0;
double* Mesh::radius = nullptr;
double* Mesh::volume = nullptr;

// Constructor implementation
Mesh::Mesh() {
  mesh = new double[n_ring]{0.0};
}

// Destructor implementation
Mesh::~Mesh() {
  delete[] mesh;
}

void Mesh::calculate_volume(const std::vector<double> inp) {

  if (radius == nullptr) {
    radius = new double[n_ring]{0.0};
  }

  for (int i = 0; i < n_ring; ++i) {
    radius[i] = inp[i];
  }

  if (volume == nullptr) {
    volume = new double[n_ring]{0.0};
  }
  double tmp = 0.0;
  for (int i = 0; i < n_ring; ++i) {
    volume[i] = M_PI * radius[i] * radius[i] - tmp;
    tmp += volume[i];
  }
}