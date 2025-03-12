#ifndef OUTER_HPP
#define OUTER_HPP

#include "cpm.hpp"
#include "mesh.hpp"

enum Problem { Eigenvalue, Fixed_source };

class Solver {
 public:
  Solver();
  ~Solver();

  Problem problem = Fixed_source;
  double Ke = 1.0;
  Group* group;

  static const int outer_max = 100;

  void solve(bool print_iter);
  void display_results();

 private:
  Mesh fission_src;
  Mesh fission_src_prev;
  Mesh total_src;
  void calculate_fission_source();
  double max_rel_error_fission();
  double max_rel_error_flux();
  void get_total_source(int g);
  void get_removal(int g, int i, double& sigr, double& removal);
  void solve_fixed_source(bool print_iter);
  void solve_eigenvalue(bool print_iter);
};

extern Solver outer;

#endif