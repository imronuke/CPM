#include "outer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iterator>

#include "cpm.hpp"
#include "mesh.hpp"

Solver::Solver() {
  group = new Group[Group::ng];
}
Solver::~Solver() {
  delete[] group;
}

void Solver::calculate_fission_source() {
  // Loop over the spatial cells (i) and energy groups (g)
  for (int i = 0; i < Mesh::n_ring; i++) {
    fission_src.mesh[i] = 0.0;
    for (int g = 0; g < Group::ng; g++) {
      fission_src.mesh[i] += group[g].nusigf.mesh[i] * group[g].flux.mesh[i];
      // printf("%.7f  %.7f  %.7f \n", fission_src.mesh[i],
      //        group[g].nusigf.mesh[i], group[g].flux.mesh[i]);
    }
  }
}

double Solver::max_rel_error_fission() {
  // Purpose: To calculate Max Relative error
  double rel = 0.0;
  double error;

  for (int i = 0; i < Mesh::n_ring; i++) {
    if (std::abs(fission_src.mesh[i]) > 1.e-20) {
      error = std::abs(fission_src.mesh[i] - fission_src_prev.mesh[i]) /
              std::abs(fission_src.mesh[i]);
      if (error > rel) {
        rel = error;
      }
    }
  }

  return rel;
}

// Get total source for group g
void Solver::get_total_source(const int g) {
  double scat_src[Mesh::n_ring]{0.0};

  for (int i = 0; i < Mesh::n_ring; i++) {
    scat_src[i] = 0.0;
    for (int h = 0; h < Group::ng; h++) {
      if (g != h) {
        scat_src[i] += group[h].sigs_in[g].mesh[i] * group[h].flux.mesh[i];
      }
    }
  }

  for (int i = 0; i < Mesh::n_ring; i++) {
    total_src.mesh[i] = group[g].chi * fission_src.mesh[i] / Ke + scat_src[i] +
                        group[g].src.mesh[i];
  }
}

double Solver::max_rel_error_flux() {
  // Purpose: To calculate Max Relative error for flux
  double rel = 0.0;
  double error;

  for (int i = 0; i < Mesh::n_ring; i++) {
    for (int g = 0; g < Group::ng; g++) {
      if (std::abs(group[g].flux.mesh[i]) > 1.e-20) {
        error = std::abs(group[g].flux.mesh[i] - group[g].flux_prev.mesh[i]) /
                std::abs(group[g].flux.mesh[i]);
        if (error > rel) {
          rel = error;
        }
      }
    }
  }
  return rel;
}

double integrate(const double* s) {
  double result = 0.0;
  for (int i = 0; i < Mesh::n_ring; i++) {
    result += Mesh::volume[i] * s[i];
  }
  return result;
}

void Solver::display_results() {
  printf("\n");
  printf("  Number of group = %3d\n", Group::ng);

  printf("  Albedo = %-11.4e    Jext = %-11.4e    Chi = %-11.4e\n",
         group[0].albedo, group[0].jext, group[0].chi);
  for (int g = 1; g < Group::ng; ++g) {
    printf("           %-11.4e           %-11.4e          %-11.4e\n",
           group[g].albedo, group[g].jext, group[g].chi);
  }

  printf("    Jnet = %-11.4e    Jout = %-11.4e    Jin = %-11.4e\n",
         group[0].jnet, group[0].jplus, group[0].jminus);
  for (int g = 1; g < Group::ng; ++g) {
    printf("           %-10.4e            %-10.4e           %-10.4e\n",
           group[g].jnet, group[g].jplus, group[g].jminus);
  }

  double tot_removal[Group::ng]{0.0};

  printf(
      "  Region      Radius      SigRem      SigScat       Source    SRC "
      "Flux     J Flux    Total Flux    Removal\n");
  for (int i = 0; i < Mesh::n_ring; ++i) {
    double sigr, removal;

    get_removal(0, i, sigr, removal);
    tot_removal[0] += removal;

    printf(
        "  %4d      %-11.4e %-11.4e  %-11.4e %-11.4e %-11.4e %-11.4e "
        "%-11.4e "
        "%-11.4e\n",
        i + 1, Mesh::radius[i], sigr, group[0].sigs_in[0].mesh[i],
        group[0].src.mesh[i], group[0].fluq.mesh[i], group[0].fluj.mesh[i],
        group[0].flux.mesh[i], removal);

    for (int g = 1; g < Group::ng; ++g) {
      get_removal(g, i, sigr, removal);
      tot_removal[g] += removal;

      printf(
          "                        %-11.4e  %-11.4e %-11.4e %-11.4e "
          "%-11.4e %-11.4e %-11.4e\n",
          sigr, group[g].sigs_in[g].mesh[i], group[g].src.mesh[i],
          group[g].fluq.mesh[i], group[g].fluj.mesh[i], group[g].flux.mesh[i],
          removal);
    }
  }

  printf("    Rtot = %-11.4e    qtot = %-11.4e\n", tot_removal[0],
         tot_removal[0] + group[0].jnet);
  for (int g = 1; g < Group::ng; ++g) {
    printf("           %-11.4e           %-11.4e\n", tot_removal[g],
           tot_removal[g] + group[g].jnet);
  }
}

void Solver::get_removal(int g, int i, double& sigr, double& removal) {
  sigr = group[g].sigt.mesh[i] - group[g].sigs_in[g].mesh[i];
  removal = sigr * group[g].flux.mesh[i] * Mesh::volume[i];
}

void Solver::fixed_source(bool print_outer) {
  calculate_fission_source();

  for (int g = 0; g < Group::ng; g++) {
    group[g].get_response_matrix();
  }

  if (print_outer) {
    printf("     ==== OUTER ITERATION ====\n");
    printf("#iter    Fission error   Flux Error\n");
  }
  for (int k = 1; k <= outer_max; k++) {
    std::copy(fission_src.mesh, fission_src.mesh + Mesh::n_ring,
              fission_src_prev.mesh);

    for (int g = 0; g < Group::ng; g++) {
      std::copy(group[g].flux.mesh, group[g].flux.mesh + Mesh::n_ring,
                group[g].flux_prev.mesh);
    }

    for (int g = 0; g < Group::ng; g++) {
      get_total_source(g);
      group[g].get_flux(total_src.mesh);
    }

    calculate_fission_source();
    double flux_error = max_rel_error_flux();
    double fsrc_error = max_rel_error_fission();

    if (print_outer) {
      printf("%3d      %-11.5e     %-11.5e\n", k, fsrc_error, flux_error);
    }

    if (flux_error < 1.e-5 && fsrc_error < 1.e-5)
      break;
  }
}

void Solver::eigenvalue(bool print_outer) {
  calculate_fission_source();

  double fsrc = integrate(fission_src.mesh);

  for (int g = 0; g < Group::ng; g++) {
    group[g].get_response_matrix();
  }

  if (print_outer) {
    printf("     ==== OUTER ITERATION ====\n");
    printf("#iter    Eigenvalue   Fission error   Flux Error\n");
  }
  for (int k = 1; k <= outer_max; k++) {
    double Ke_prev = Ke;
    double fsrc_prev = fsrc;

    std::copy(fission_src.mesh, fission_src.mesh + Mesh::n_ring,
              fission_src_prev.mesh);

    for (int g = 0; g < Group::ng; g++) {
      std::copy(group[g].flux.mesh, group[g].flux.mesh + Mesh::n_ring,
                group[g].flux_prev.mesh);
    }

    for (int g = 0; g < Group::ng; g++) {
      get_total_source(g);
      group[g].get_flux(total_src.mesh);
    }

    calculate_fission_source();
    fsrc = integrate(fission_src.mesh);
    Ke = Ke_prev * fsrc / fsrc_prev;
    double flux_error = max_rel_error_flux();
    double fsrc_error = max_rel_error_fission();

    if (print_outer) {
      printf("%3d      %-10.6f  %-11.5e     %-11.5e\n", k, Ke, fsrc_error,
             flux_error);
    }

    if (flux_error < 1.e-5 && fsrc_error < 1.e-5)
      break;
  }
}
