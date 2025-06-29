#include "src/functions.hpp"
#include "src/molecule.hpp"
#include "src/volume.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>

int main() {
  int seed, number_of_molecules = 0;
  std::cin >> number_of_molecules >> seed;

  std::vector<Molecule> molecules_set;
  molecules_set = molecule_set_generate(number_of_molecules, seed);

  Volume volume(300.0, 25.0, 1.0, molecules_set);
  for (int i = 0; i < volume.get_space().get_amount_of_molecules(); i++) {
    auto &mol = volume.get_space().get_molecule(i);
    std::cout << "Molecule " << mol.get_id() << ": ";
    mol.print_coordinate();
  }

  return 0;
}
