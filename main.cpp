// #include "src/functions.hpp"
#include "src/molecule.hpp"
#include "src/space.hpp"
#include "src/volume.hpp"
#include <cstdlib>
#include <iostream>

int main() {
  int seed, number_of_molecules = 0;
  std::cin >> number_of_molecules >> seed;

  Volume volume(number_of_molecules, seed);

  for (int i = 0; i < volume.get_space().get_amount_of_molecules(); i++) {
    auto &mol = volume.get_space().get_molecule(i);
    std::cout << "Molecule " << mol.get_id() << ":" << std::endl;
    mol.print_coordinate();
    mol.print_velocity();
  }

  double v_sum, v2_sum = 0.0;
  std::array<double, 3> cur_vel;
  for (int i = 0; i < volume.get_space().get_amount_of_molecules(); i++) {
    cur_vel = volume.get_space().get_molecule(i).get_velocity();
    for (int j = 0; j < 3; j++) {
      v_sum += cur_vel[j];
      v2_sum += cur_vel[j] * cur_vel[j];
    }
  }
  std::cout << "V_center_mass: " << v_sum << " V^2: " << v2_sum << std::endl;

  return 0;
}
