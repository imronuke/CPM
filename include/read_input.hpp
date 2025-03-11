#ifndef READ_INPUT_HPP
#define READ_INPUT_HPP

#include <string>

#include "cpm.hpp"
#include "outer.hpp"

void open_xml(std ::string input_xml);
void read_problem_definition(int& n_ring, int& n_gauss, int& ng);
void read_data(Solver& outer);

#endif