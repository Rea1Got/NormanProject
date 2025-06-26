#include "space.hpp"
#include "molecule.hpp"
#include <iostream>

void Space::print_space() {
  std::cout << "Current space { x: " << lenght_x << " y: " << lenght_y
            << " z: " << lenght_z << " }" << std::endl;
}

void Space::change_space(int x, int y, int z) {
  lenght_x = x;
  lenght_y = y;
  lenght_z = z;
}

void Volume::print_full_information() {
  std::cout << "Space parameters: " << std::endl;
  print_space();
  std::cout << "Molecules parameters: " << std::endl;
  for (int i = 0; i < molecules.size(); i++) {
    molecules[i].print_full_information();
  }
}
