#ifndef CPM_HPP
#define CPM_HPP

#include "mesh.hpp"  // Include Mesh header

class Group {
 public:
  static int ng;
  static int n_gauss;

  unsigned int g;
  double chi;
  double jplus;
  double jminus;
  double jnet;
  double albedo;
  double jext;

  Mesh sigt;
  Mesh* sigs_in;  // Incoming scattering from other group scattering
  Mesh nusigf;
  Mesh src;
  Mesh fluq;
  Mesh fluj;
  Mesh flux;
  Mesh flux_prev;
  Mesh fiss_src;

  Group();  // Constructor
  ~Group();

  void get_response_matrix();
  void get_flux(double* source);
  void solve_collision_probability();

 private:
  double** LU;
  double** pij;
  double* gamma;

  static void LU_decompose(double** A, int n);
  static double* LU_solve(double** LU, const double* b, int n);
  static double KI3(double xx);
  static double yi(double y, double gam, double alpha);
  void get_pij_annular(double** p, double* gam);
};

extern Group* group;

constexpr int ngauss_max = 6;
#ifdef _GG
// Gauss Quadrature
constexpr double gauss_quadrature[ngauss_max][ngauss_max] = {
    {0.0, 0.00000000, -0.57735027, -0.77459667, -0.86113631, -0.90617985},
    {2.00000000, 0.0, 0.57735027, 0.0000000, -0.33998104, -0.53846931},
    {1.00000000, 1.00000000, 0.0, 0.77459667, 0.33998104, 0.00000000},
    {0.55555555, 0.88888889, 0.55555555, 0.0, 0.86113631, 0.53846931},
    {0.34785485, 0.65214515, 0.65214515, 0.34785485, 0.0, 0.90617985},
    {0.23692688, 0.47862867, 0.56888889, 0.47862867, 0.23692688, 0.0}};
#else
// Gauss-Jacobi Quadrature
constexpr double gauss_quadrature[ngauss_max][ngauss_max] = {
    {0.0, 0.55555556, 0.87393877, 0.95491150, 0.98046718, 0.99029084},
    {2.00000000, 0.0, 0.28606124, 0.65127016, 0.82660307, 0.90725799},
    {0.72783448, 1.27216552, 0.0, 0.16932809, 0.47704397, 0.68412769},
    {0.27930792, 0.91696442, 0.80372766, 0.0, 0.11094751, 0.35681753},
    {0.12472388, 0.51939018, 0.81385828, 0.54202764, 0.0, 0.07803490},
    {0.06299166, 0.29563548, 0.58554794, 0.66869856, 0.38712636, 0.0}};
#endif

#endif
