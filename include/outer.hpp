#ifndef OUTER_HPP
#define OUTER_HPP

#include "cpm.hpp"
#include "mesh.hpp"

enum Problem { eigenvalue, fixed_source };

class Solver {
 public:
  Solver();
  ~Solver();

  Problem problem = fixed_source;
  double Ke = 1.0;
  Group* group;

  static const int outer_max = 100;

  void fixed_source(bool print_iter);
  void eigenvalue(bool print_outer);
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
};

extern Solver outer;

#endif