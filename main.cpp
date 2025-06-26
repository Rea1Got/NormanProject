#include "src/molecule.hpp"

int main() {
  Molecule water;
  water.get_full_information();

  unsigned int id = 1;
  std::array<float, 3> coord = {1.0, 2.0, 3.0};
  std::array<float, 3> veloc = {0.1, 0.2, 0.3};
  std::array<float, 3> accel = {0.01, 0.01, 0.01};

  Molecule ethanol(id, coord, veloc, accel);
  ethanol.get_full_information();
  return 0;
}
