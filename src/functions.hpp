#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP
#include "molecule.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>

#define MAX_LENGTH 100.0

std::vector<Molecule> molecule_set_generate(int number_of_molecules, int seed) {
  std::srand(seed);
  std::vector<Molecule> molecules_set;
  // normalization to unit volume
  for (int i = 0; i < number_of_molecules; i++) {
    molecules_set.push_back(Molecule());
    std::array<double, 3> coords = {
        static_cast<double>(std::rand()) / RAND_MAX * MAX_LENGTH,
        static_cast<double>(std::rand()) / RAND_MAX * MAX_LENGTH,
        static_cast<double>(std::rand()) / RAND_MAX * MAX_LENGTH};
    molecules_set[i].set_coordinate(coords);
  }
  return molecules_set;
}

#endif
