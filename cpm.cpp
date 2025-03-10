#include "cpm.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>

#include "mesh.hpp"

int Group::ng = 0;
int Group::n_gauss = 0;
const int& n_ring = Mesh::n_ring;

// Constructor Implementation
Group::Group()
    : sigt(), nusigf(), src(), fluq(), fluj(), flux(), flux_prev(), fiss_src() {
  chi = 0.0;
  jplus = 0.0;
  jminus = 0.0;
  jnet = 0.0;
  albedo = 0.0;
  jext = 0.0;

  sigs_in = new Mesh[ng];

  for (int g = 0; g < ng; ++g) {
    gamma = new double[n_ring]{0.0};
    LU = new double*[n_ring];
    pij = new double*[n_ring];
    for (int i = 0; i < n_ring; ++i) {
      LU[i] = new double[n_ring]{0.0};
      pij[i] = new double[n_ring]{0.0};
    }
  }

  for (int i = 0; i < n_ring; ++i) {
    flux.mesh[i] = 1.0;
  }
}

Group::~Group() {
  for (int i = 0; i < n_ring; ++i) {
    delete[] LU[i];
    delete[] pij[i];
  }

  delete[] LU;
  delete[] pij;
  delete[] gamma;
  delete[] sigs_in;
}

// Function to perform LU decomposition in-place
void Group::LU_decompose(double** A, int n) {
  for (int i = 0; i < n; ++i) {
    // Check for zero pivot
    if (std::abs(A[i][i]) < 1.0e-10) {
      throw std::runtime_error(
          "LU decomposition failed: Zero pivot encountered.");
    }

    // Compute the multipliers and update the lower and upper parts of the
    // matrix
    for (int j = i + 1; j < n; ++j) {
      A[j][i] /= A[i][i];  // Store the multiplier in the lower part (L)
      for (int k = i + 1; k < n; ++k) {
        A[j][k] -= A[j][i] * A[i][k];  // Update the upper part (U)
      }
    }
  }
}

// Function to solve LUx = b using forward and backward substitution
// double* LU_solve(const double* const LU[], const double* b, int n) {
double* Group::LU_solve(double** LU, const double* b, int n) {
  // Forward substitution: Solve Ly = b (L is unit lower triangular)
  double* x = new double[n];
  double y[n];
  for (int i = 0; i < n; ++i) {
    y[i] = b[i];
    for (int j = 0; j < i; ++j) {
      y[i] -= LU[i][j] * y[j];
    }
  }

  // Backward substitution: Solve Ux = y (U is upper triangular)
  for (int i = n - 1; i >= 0; --i) {
    x[i] = y[i];
    for (int j = i + 1; j < n; ++j) {
      x[i] -= LU[i][j] * x[j];
    }
    x[i] /= LU[i][i];
  }

  return x;
}

// Function to calculate the 3rd order Bickely function
double Group::KI3(double xx) {
  double f, x;

  // Take absolute value of xx
  x = std::abs(xx);

  if (x < 1.e-20) {
    f = 0.25 * M_PI;
  } else if (x < 0.25) {
    if (x < 0.05) {
      f = (0.7266088 * x - 0.9990226) * x + 0.7853961;
    } else if (x < 0.10) {
      f = (0.6466375 * x - 0.9912340) * x + 0.7852024;
    } else if (x < 0.15) {
      f = (0.5856605 * x - 0.9791293) * x + 0.7845986;
    } else if (x < 0.20) {
      f = (0.5346648 * x - 0.9638914) * x + 0.7834577;
    } else {
      f = (0.4907827 * x - 0.9463843) * x + 0.7817094;
    }
  } else if (x < 0.50) {
    if (x < 0.30) {
      f = (0.4521752 * x - 0.9271152) * x + 0.7793031;
    } else if (x < 0.35) {
      f = (0.4177388 * x - 0.9064822) * x + 0.7762107;
    } else if (x < 0.40) {
      f = (0.3869945 * x - 0.8849865) * x + 0.7724519;
    } else if (x < 0.45) {
      f = (0.3590753 * x - 0.8626685) * x + 0.7679903;
    } else {
      f = (0.3338676 * x - 0.8400133) * x + 0.7628988;
    }
  } else if (x < 1.0) {
    if (x < 0.60) {
      f = (0.2998569 * x - 0.8054172) * x + 0.7540982;
    } else if (x < 0.70) {
      f = (0.2609154 * x - 0.7587821) * x + 0.7401279;
    } else if (x < 0.80) {
      f = (0.2278226 * x - 0.7125290) * x + 0.7239594;
    } else if (x < 0.90) {
      f = (0.1994999 * x - 0.6672761) * x + 0.7058777;
    } else {
      f = (0.1751248 * x - 0.6234536) * x + 0.6861762;
    }
  } else if (x < 3.0) {
    if (x < 1.4) {
      f = ((-0.05337485 * x + 0.3203223) * x - 0.7538355) * x + 0.7247294;
    } else if (x < 1.8) {
      f = ((-0.03146833 * x + 0.2295280) * x - 0.6279752) * x + 0.6663720;
    } else if (x < 2.2) {
      f = ((-0.01906198 * x + 0.1631667) * x - 0.5094124) * x + 0.5956163;
    } else if (x < 2.6) {
      f = ((-0.01174752 * x + 0.1152418) * x - 0.4046007) * x + 0.5191031;
    } else {
      f = ((-0.007328415 * x + 0.08097913) * x - 0.3159648) * x + 0.4425954;
    }
  } else if (x < 5.0) {
    if (x < 3.4) {
      f = ((-0.004617254 * x + 0.05669960) * x - 0.2434341) * x + 0.3703178;
    } else if (x < 3.8) {
      f = (0.007923547 * x - 0.07158569) * x + 0.1684022;
    } else if (x < 4.2) {
      f = (0.005095111 * x - 0.05016344) * x + 0.1278307;
    } else if (x < 4.6) {
      f = (0.003286040 * x - 0.03501524) * x + 0.09611422;
    } else {
      f = (0.002126242 * x - 0.02437465) * x + 0.07170491;
    }
  } else if (x < 9.0) {
    if (x < 5.8) {
      f = (0.001123687 * x - 0.01425519) * x + 0.04616317;
    } else if (x < 6.6) {
      f = (4.762937E-4 * x - 6.810124E-3) * x + 0.02475115;
    } else if (x < 7.4) {
      f = (2.031843E-4 * x - 3.232035E-3) * x + 0.01302864;
    } else if (x < 8.2) {
      f = (8.701440E-5 * x - 1.524126E-3) * x + 6.749972E-3;
    } else {
      f = (3.742673E-5 * x - 7.157367E-4) * x + 3.454768E-3;
    }
  } else {
    std::cerr << "No value for x = " << x << " in the KI3 function"
              << std::endl;
    throw std::out_of_range("Value out of range in KI3 function");
  }

  return f;
}

void Group::get_pij_annular(double** p, double* gam) {
  double s[n_ring + 1][n_ring + 1] = {};  // initialize to 0
  double tau[n_ring]{0.0};

  double y_start = 0.0;
  constexpr auto& quad = gauss_quadrature;

  // Sweep upward over y-axis levels
  const double* radius = Mesh::radius;
  for (int i = 0; i < n_ring; ++i) {
    double dely = radius[i] - y_start;

    // Loop over Gauss quadrature points
    for (int k = 0; k <= n_gauss; ++k) {
      // Calculate Gauss points and weights
#ifdef _GG
      double y_point = y_start + 0.5 * (quad[n_gauss - 1][k] + 1.0) * dely;
      double y_weight = quad[k][n_gauss - 1] * dely;
#else
      double y_point = y_start + quad[k][n_gauss] * dely;
      double y_weight = quad[n_gauss][k] * dely;
#endif

      double x_start = 0.0;
      double tau_prev = 0.0;

      // Determine optical distance
      for (int j = i; j < n_ring; ++j) {
        // std::cout << i << " " << j << "\n";
        double xj = sqrt(radius[j] * radius[j] - y_point * y_point);
        double delxj = xj - x_start;
        tau[j] = tau_prev + sigt.mesh[j] * delxj;
        x_start = xj;
        tau_prev = tau[j];

        // Accumulate sij
        for (int l = i; l <= j; ++l) {
          double taup = tau[j] + tau[l];
          double taum = tau[j] - tau[l];
          s[l + 1][j + 1] += y_weight * (KI3(taup) - KI3(taum));
        }  // loop over origins
      }  // loop over destinations
    }  // loop over gauss points
    y_start = radius[i];
  }  // loop y-axis levels over levels

#ifdef _CHECK
  std::cout << "SIJ" << std::endl;
  for (int i = 0; i <= n_ring; ++i) {
    for (int j = 0; j <= n_ring; ++j) {
      printf("%12.4E", s[i][j]);
    }
    printf("\n");
  }
#endif

  const double* volume = Mesh::volume;
  for (int i = n_ring; i >= 1; --i) {
    s[i][i - 1] = s[i - 1][i];
    for (int j = i; j >= 1; --j) {
      p[j - 1][i - 1] = s[j][i] + s[j - 1][i - 1] - (s[j][i - 1] + s[j - 1][i]);
      p[i - 1][j - 1] = p[j - 1][i - 1];
    }
    double sigvol = volume[i - 1] * sigt.mesh[i - 1];
    gam[i - 1] = sigvol;
    p[i - 1][i - 1] = sigvol + p[i - 1][i - 1];
  }

  double rsb = 2.0 / (M_PI * radius[n_ring - 1]);
  for (int i = 0; i < n_ring; ++i) {
    double sum = 0.0;
    for (int j = 0; j < n_ring; ++j) {
      sum += p[i][j];
    }
    gam[i] = (gam[i] - sum) * rsb;
  }
#ifdef _CHECK
  std::cout << "SIJ" << std::endl;
  for (int i = 0; i <= n_ring; ++i) {
    for (int j = 0; j <= n_ring; ++j) {
      printf("%12.4E", s[i][j]);
    }
    printf("\n");
  }
#endif
#ifdef _CHECK
  std::cout << "PIJ" << std::endl;
  for (int i = 0; i < n_ring; ++i) {
    for (int j = 0; j < n_ring; ++j) {
      printf("%12.4E", p[i][j]);
    }
    printf("\n");
  }
  std::cout << "GAMMA" << std::endl;
  for (int i = 0; i < n_ring; ++i) {
    printf("%12.4E\n", gam[i]);
  }
  exit(0);
#endif
}

void Group::get_response_matrix() {
  get_pij_annular(pij, gamma);

  const double* volume = Mesh::volume;
  for (int i = 0; i < n_ring; i++) {
    double c = sigs_in[g].mesh[i] / sigt.mesh[i];
    double sigvol = volume[i] * sigt.mesh[i];

    for (int j = 0; j < n_ring; j++) {
      if (i == j) {
        LU[j][i] = sigvol - pij[i][j] * c;
      } else {
        LU[j][i] = -pij[i][j] * c;
      }
    }
  }

#ifdef _CHECK
  std::cout << "MATRIX A" << std::endl;
  for (int i = 0; i < n_ring; ++i) {
    for (int j = 0; j < n_ring; ++j) {
      std::cout << LU[i][j] << " ";
    }
    std::cout << std::endl;
  }
#endif

  LU_decompose(LU, n_ring);

#ifdef _CHECK
  std::cout << "MATRIX LU" << std::endl;
  for (int i = 0; i < n_ring; ++i) {
    for (int j = 0; j < n_ring; ++j) {
      std::cout << LU[i][j] << " ";
    }
    std::cout << std::endl;
  }
#endif
}

// Function to calculate yi(alpha)
double Group::yi(double y, double gam, double alpha) {
  // Calculate the result
  double f = y / (1.0 - (1.0 - gam) * alpha);
  return f;
}

void Group::get_flux(double* source) {
  double sb = 2.0 * M_PI * Mesh::radius[n_ring - 1];

  // create RHS and solve for X
  double* x[n_ring];
  double b[n_ring];
  for (int i = 0; i < n_ring; ++i) {
    for (int j = 0; j < n_ring; ++j) {
      b[j] = pij[i][j] / sigt.mesh[i];
    }
    x[i] = new double[n_ring]{0.0};
    x[i] = LU_solve(LU, b, n_ring);
  }

  // solve for Y
  double* y = new double[n_ring];
  y = LU_solve(LU, gamma, n_ring);

  double sum_gam = 0.0;
  double xk[n_ring];
  for (int i = 0; i < n_ring; ++i) {
    // Calculate capital gamma
    sum_gam += (sigt.mesh[i] - sigs_in[g].mesh[i]) * y[i] * Mesh::volume[i];
    xk[i] = 0.25 * sb * Mesh::volume[i] * y[i];
  }

  for (int i = 0; i < n_ring; ++i) {
    // Get flux due to sources
    fluq.mesh[i] = 0.0;
    for (int k = 0; k < n_ring; ++k) {
      fluq.mesh[i] +=
          (x[k][i] + albedo * xk[k] * yi(y[i], sum_gam, albedo)) * source[k];
    }

    // Get flux due to currents
    fluj.mesh[i] = yi(y[i], sum_gam, albedo) * jext;

    // Get total flux
    flux.mesh[i] = fluq.mesh[i] + fluj.mesh[i];
  }

  double sumxk = 0.0;
  for (int i = 0; i < n_ring; ++i) {
    sumxk += xk[i] * source[i];
  }

  jplus = (sumxk + (1.0 - sum_gam) * jext) / (1.0 - albedo * (1.0 - sum_gam));
  jminus = (albedo * sumxk + jext) / (1.0 - albedo * (1.0 - sum_gam));
  jnet = jplus - jminus;

  for (int i = 0; i < n_ring; ++i) {
    delete[] x[i];
  }
  delete[] y;

#ifdef _CHECK
  std::cout << "CAPITAL GAMMA : " << sum_gam << std::endl;
  std::cout << "SMALL GAMMA  SRC Flux J Flux Total Flux : " << std::endl;
  for (int i = 0; i < n_ring; ++i) {
    printf("%12.5E%12.5E%12.5E%12.5E\n", gamma[i], fluq.mesh[i], fluj.mesh[i],
           flux.mesh[i]);
  }
#endif
}

void Group::solve_collision_probability() {
  for (int g = 0; g < Group::ng; g++) {
    get_response_matrix();
    get_flux(src.mesh);
  }
}