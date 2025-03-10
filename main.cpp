#include <chrono>
#include <cstdio>
#include <cstdio>  // For tmpfile()
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "cpm.hpp"
#include "mesh.hpp"
#include "outer.hpp"
#include "read_input.hpp"

// Get random number in range(a, b)
// double get_random_number(const double a, const double b) {
//     // Seed the random number generator with the current time
//     std::random_device rd;
//     std::mt19937 gen(rd());

//     // Define a uniform real distribution between 0 and 1
//     std::uniform_real_distribution<> dist(0.0, 1.0);

//     // Generate and print 10 random floating-point numbers between 0 and 1
//     double random_number = dist(gen);

//     return (b - a) * random_number + a;
// }

void remove_comments(const std::string& filename, char* tempFilename) {
  std::ifstream input(filename);
  if (!input) {
    std::cerr << "Error: " << filename << " not found!" << std::endl;
    throw std::runtime_error("Input file not found!");
  }

  std::string buff;
  std::vector<std::string> initstr;

  // Read lines and remove comments (!)
  while (getline(input, buff)) {
    size_t pos = buff.find('!');
    if (pos != std::string::npos) {
      buff = buff.substr(0, pos);  // Remove everything from '!' onwards
    }
    initstr.push_back(buff);
  }
  input.close();

  std::ofstream tempFile(tempFilename);
  if (!tempFile) {
    std::cerr << "Error: Unable to create temporary file!" << std::endl;
    throw std::runtime_error("Unable to create temporary file!");
  }

  // Write modified content to the temporary file
  for (const auto& line : initstr) {
    tempFile << line << std::endl;  // Write cleaned lines
  }
  tempFile.close();
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Error: You must provide input.xml path" << std::endl;
    std::cerr << "Usage: ./cpm path/to/input.xml" << std::endl;
    return 1;  // Return an error code
  }

  // Create a temporary file where ! mark is removed
  char tempFilename[L_tmpnam];  // Buffer to store temp filename

  std::string filename = argv[1];
  remove_comments(filename, tempFilename);

  // I/O directory
  const std::string input_xml = std::string(tempFilename);
  open_xml(input_xml);

  read_problem_definition(Mesh::n_ring, Group::n_gauss, Group::ng);

  Outer outer;

  read_data(outer);

  if (outer.problem == fixed_source) {
    outer.solve_fixed_source(true);
  } else {
    outer.solve_eigenvalue(true);
  }

  outer.display_results();

  // bool performance_test = false;

  // if (performance_test) {
  //     constexpr int n_calc = 10000;
  //     constexpr int n_print = 1000;
  //     auto sigt = outer.group[0].sigt;

  //     printf("\n");
  //     printf("  === PERFORMANCE TEST === \n");

  //     auto start = std::chrono::high_resolution_clock::now();

  //     int p = 0;
  //     for (int k = 1; k <= n_calc; ++k) {
  //         auto rx = get_random_number(-0.01, 0.01);

  //         for (int i = 0; i < Mesh::n_ring; ++i) {
  //             outer.group[0].sigt.mesh[i] = sigt.mesh[i] + rx *
  //             sigt.mesh[i];
  //         }

  //         for (int g = 0; g < Group::ng; ++g) {
  //             outer.group[g].solve_collision_probability();
  //         }

  //         if (k % n_print == 0) {
  //             p += n_print;
  //             printf("  Now is reaching %7d calculations\n", p);
  //         }
  //     }

  //     auto end = std::chrono::high_resolution_clock::now();
  //     std::chrono::duration<double> duration = end - start;
  //     std::cout << "Elapsed time: " << duration.count() << " seconds\n";
  // }

  remove(tempFilename);
}