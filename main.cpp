#include "src/molecule.hpp"
#include "src/space.hpp"
#include <iostream>
#include <vector>

int main() {
  Molecule water;
  water.print_full_information();

  unsigned int id = 1;
  std::array<float, 3> coord = {1.0, 2.0, 3.0};
  std::array<float, 3> veloc = {0.1, 0.2, 0.3};
  std::array<float, 3> accel = {0.01, 0.01, 0.01};

  Molecule ethanol(id, coord, veloc, accel);
  ethanol.print_full_information();
  std::cout << "Molecule test complete!" << "\n\n";

  std::vector<Molecule> molecules_set;
  for (int i = 0; i < 5; i++) {
    std::array<float, 3> test = {0, 0, 0};
    test[0] += i;
    test[1] += i;
    test[2] += i;
    molecules_set.push_back(Molecule(i, test, test, test));
  }

  Volume flaska = Volume(1, 2, 3, molecules_set);
  flaska.print_full_information();
  flaska.change_space(4, 5, 6);
  flaska.print_space();
  std::cout << "Space test complete!" << "\n\n";
  return 0;
}
