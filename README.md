## Overview

Welcome to the **Collision Probability Method (CPM)** repository! This project features an integral neutron transport solver for circular geometry using the CPM technique. The method implemented in this solver is based on the COLPROB code, as detailed in the book by Prof. Stammler and Prof. Abbate.

In this repo, the implementation is extended to multi-group and to solve the eigenvalue problem.

## Installation

To get started with the CPM solver, clone the repository and follow the installation instructions:

```bash
git clone https://github.com/imronuke/CPM.git
cd CPM
make
```

## Usage

Here's a quick example of how to use the CPM solver:

```bash
# Example command to run the solver
./cpm example/test_1.xml
```

## Dependencies

You just need g++ or other c++ compilers.


## License

This project is licensed under the [MIT License](https://github.com/imronuke/CPM/blob/main/LICENSE).


## Acknowledgements

This CPM code was converted from the Fortran version that I submitted for homework. Functions ```get_pij_annular``` (to calculate the PIJ matrix) and ```KI3``` (to calculate the 3rd order Bickley-Naylor function) were provided in FORTRAN-77 as homework supplements. Pugixml to read the XML input file is from this [repository](https://github.com/zeux/pugixml).

---

Feel free to explore, use, and contribute to the **Collision Probability Method (CPM)** project. Happy coding!