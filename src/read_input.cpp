#include "read_input.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

#include "cpm.hpp"
#include "mesh.hpp"
#include "outer.hpp"
#include "pugixml.hpp"

// Function to trim whitespace from both ends of the string
std::string trim(const std::string& str) {
  // Find the first non-whitespace character
  auto start = str.find_first_not_of(" \t\n\r\f\v");

  // If the string is empty or only contains whitespace, return an empty
  // string
  if (start == std::string::npos)
    return "";

  // Find the last non-whitespace character
  auto end = str.find_last_not_of(" \t\n\r\f\v");

  // Extract the substring without the leading and trailing whitespace
  return str.substr(start, end - start + 1);
}

// Function to convert a string to lowercase
std::string LowerCase(const std::string& str) {
  std::string result = str;  // Make a copy of the string to modify
  std::transform(result.begin(), result.end(), result.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return result;
}

// Function to trim and convert a string to lowercase
std::string std_str(const std::string& str) {
  std::string result = trim(LowerCase(str));
  return result;
}

std::tuple<std::vector<double>, int> string_to_vec(const std::string& input) {
  // Vector to store the result in one-dimensional form
  std::vector<double> result;
  int n_data = 0;

  // Use stringstream to parse the entire string
  std::stringstream ss(input);
  std::string line;

  // Process each line (representing each row)
  while (std::getline(ss, line)) {
    std::stringstream lineStream(line);
    double value;

    // Process each number in the line and convert it to double
    while (lineStream >> value) {
      result.push_back(value);
      n_data++;  // Increment column count for this row
    }

    // If this is the first row, set numCols (column count)
  }

  return {result, n_data};
}

// XML input file
pugi::xml_document doc;
pugi::xml_parse_result result;

void open_xml(std ::string input_xml) {
  result = doc.load_file(input_xml.c_str());

  if (!result) {
    std::cerr << "XML [" << input_xml << "] parsed with errors, attr value: ["
              << doc.child("node").attribute("attr").value() << "]\n";
    std::cerr << "Error description: " << result.description() << "\n";
    std::cerr << "Error offset: " << result.offset << " (error at [..."
              << (input_xml.c_str() + result.offset) << "]\n\n";
    throw std::runtime_error("Invalid input detected!");
  }
}

void read_problem_definition(int& n_ring, int& n_gauss, int& ng) {
  const auto& head = doc.child("cpm");

  n_ring = std::stoi(head.child("data").attribute("n_ring").value());
  n_gauss = std::stoi(head.child("data").attribute("n_gauss").value());
  ng = std::stoi(head.child("data").attribute("ng").value());
}

void read_data(Solver& outer) {
  const auto& head = doc.child("cpm");

  const int& n_ring = Mesh::n_ring;
  const int& ng = Group::ng;
  const auto& group = outer.group;

  if (head.child("mode")) {
    const std::string mode = std_str(head.child_value("mode"));
    if (mode == "eigenvalue") {
      outer.problem = Eigenvalue;
    } else if (mode == "fixed_source") {
      outer.problem = Fixed_source;
    } else {
      std::cerr << "Invalid input <mode>" << std::endl;
      throw std::runtime_error("Invalid input detected!");
    }
  }

  {
    auto [tmp, n_data] = string_to_vec(head.child_value("rings"));
    if (n_data != n_ring) {
      std::cerr << "Invalid number of data in <rings>. Expected : " << n_ring
                << ". Got: " << n_data << std::endl;
      throw std::runtime_error("Invalid input detected!");
    } else {
      for (int i = 0; i < n_ring; ++i) {}
    }

    Mesh::calculate_volume(tmp);
  }

  {
    auto [tmp, n_data] = string_to_vec(head.child_value("xsec"));

    int n_expect = 2 * n_ring * ng + ng * ng * n_ring;

    if (n_data != n_expect) {
      std::cerr << "Invalid number of data in <xsec>. Expected : " << n_expect
                << ". Got: " << n_data << std::endl;
      throw std::runtime_error("Invalid input detected!");
    } else {
      int n1 = 0;
      int n2 = n_ring;
      for (int g = 0; g < ng; ++g) {
        std::vector<double> sliced(tmp.begin() + n1, tmp.begin() + n2);
        for (int i = 0; i < n_ring; ++i) {
          group[g].sigt.mesh[i] = sliced[i];
        }
        n1 = n2;
        n2 = n2 + n_ring;
      }

      for (int g = 0; g < ng; ++g) {
        std::vector<double> sliced(tmp.begin() + n1, tmp.begin() + n2);
        for (int i = 0; i < n_ring; ++i) {
          group[g].nusigf.mesh[i] = sliced[i];
        }
        n1 = n2;
        n2 = n2 + n_ring;
      }

      for (int g = 0; g < ng; ++g) {
        group[g].g = g;
        for (int h = 0; h < ng; ++h) {
          std::vector<double> sliced(tmp.begin() + n1, tmp.begin() + n2);
          for (int i = 0; i < n_ring; ++i) {
            group[g].sigs_in[h].mesh[i] = sliced[i];
          }
          n1 = n2;
          n2 = n2 + n_ring;
        }
      }
    }
  }

  if (outer.problem == Fixed_source) {
    auto [tmp, n_data] = string_to_vec(head.child_value("src"));
    if (n_data != ng * n_ring) {
      std::cerr << "Invalid number of data in <src>. Expected : " << ng * n_ring
                << ". Got: " << n_data << std::endl;
      throw std::runtime_error("Invalid input detected!");
    } else {
      int n1 = 0;
      int n2 = n_ring;
      for (int g = 0; g < ng; ++g) {
        std::vector<double> sliced(tmp.begin() + n1, tmp.begin() + n2);
        for (int i = 0; i < n_ring; ++i) {
          group[g].src.mesh[i] = sliced[i];
        }
        n1 = n2;
        n2 = n2 + n_ring;
      }
    }
  }

  {
    auto [tmp, n_data] = string_to_vec(head.child_value("chi"));
    if (n_data != ng) {
      std::cerr << "Invalid number of data in <chi>. Expected : " << ng
                << ". Got: " << n_data << std::endl;
      throw std::runtime_error("Invalid input detected!");
    } else {
      double sum = 0.0;
      for (int g = 0; g < ng; ++g) {
        sum += tmp[g];
        group[g].chi = tmp[g];
      }

      if (std::abs(sum - 1.0) > 1.0e-5 && outer.problem == Eigenvalue) {
        std::cerr << "Sum of fission spectrum is equal to : " << sum
                  << ". It should be equal to one\n";
        throw std::runtime_error("Invalid input detected!");
      }
    }
  }

  {
    auto [tmp, n_data] = string_to_vec(head.child_value("albedo"));
    if (n_data != ng) {
      std::cerr << "Invalid number of data in <albedo>. Expected : " << ng
                << ". Got: " << n_data << std::endl;
      throw std::runtime_error("Invalid input detected!");
    } else {
      for (int g = 0; g < ng; ++g) {
        group[g].albedo = tmp[g];
      }
    }
  }

  {
    auto [tmp, n_data] = string_to_vec(head.child_value("j_ext"));
    if (outer.problem == Eigenvalue) {
      if (n_data > 0) {
        std::cerr << "External current is not necessary for eigenvalue "
                     "problem";
        throw std::runtime_error("Invalid input detected!");
      }
    } else {
      if (n_data != ng) {
        std::cerr << "Invalid number of data in <j_ext>. Expected : " << ng
                  << ". Got: " << n_data << std::endl;
        throw std::runtime_error("Invalid input detected!");
      } else {
        for (int g = 0; g < ng; ++g) {
          group[g].jext = tmp[g];
        }
      }
    }
  }
}