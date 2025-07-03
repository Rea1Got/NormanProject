// #include "src/functions.hpp"
#include "src/molecule.hpp"
#include "src/space.hpp"
#include "src/volume.hpp"
#include <cstdlib>
#include <iostream>

int main() {
  int seed, number_of_molecules = 0;
  std::cin >> number_of_molecules >> seed;
  double temperature = 0;
  std::cin >> temperature;
  Volume volume(number_of_molecules, seed, temperature, 2.5, 0.34,
                120 * 1.34e-23);

  std::vector<std::array<double, 3>> force;
  force = volume.calculate_force();
  for (int i = 0; i < volume.get_space().get_amount_of_molecules(); i++) {
    Molecule &mol = volume.get_space().get_molecule(i);
    mol.print_full_information();
    std::cout << "Current force: " << force[i][0] << " " << force[i][1] << " "
              << force[i][2] << "\n\n";
  }
  std::array<double, 3> pot_energy = force.back();
  std::cout << "Pot_energy: " << pot_energy[0] << "\n";

  double v_sum = 0.0, v2_sum = 0.0;
  for (int i = 0; i < volume.get_space().get_amount_of_molecules(); i++) {
    auto velocity = volume.get_space().get_molecule(i).get_velocity();
    for (int j = 0; j < 3; j++) {
      v_sum += velocity[j];
      v2_sum += velocity[j] * velocity[j];
    }
  }
  std::cout << "Sum V: " << v_sum << " Sum V^2: " << v2_sum << "\n\n";

  return 0;
}
