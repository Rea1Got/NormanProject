#include "include/json.hpp" //nlohmann
#include "src/molecule.hpp"
#include "src/space.hpp"
#include "src/volume.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

using json = nlohmann::json;

int main() {
  std::ifstream cfg_file("cfg/cfg.json");
  json cfg;
  cfg_file >> cfg;

  std::string velocity_file = cfg["velocity_file"];
  std::string coord_abs_file = cfg["coord_abs_file"];
  int seed = cfg["seed"];
  int num_molecules = cfg["num_molecules"];
  double temperature = cfg["temperature"];
  double dt = cfg["dt"];
  int total_steps = cfg["total_steps"];
  double r_cut = cfg["r_cut"];
  double sigma = cfg["sigma"];
  double epsilon = cfg["epsilon"];
  int snapshot = cfg["snapshot"]; // how often you want to see system state
  double x = cfg["length_x"];
  double y = cfg["length_y"];
  double z = cfg["length_z"];
  double density = cfg["density"];
  double mass_real = cfg["mass_real"];
  double sigma_real = cfg["sigma_real"];

  Volume volume(num_molecules, seed, density, mass_real, sigma_real,
                temperature, r_cut, sigma, epsilon);
  Space &space = volume.get_space();

  std::cout << "===== Начало симуляции =====\n";
  std::cout << "Молекул: " << num_molecules << "\n";
  std::cout << "Шаг по времени: " << dt << "\n";
  std::cout << "Температура: " << temperature << "\n";
  std::cout << "Длина объема: " << volume.get_length_x() << "\n\n";

  for (int step = 0; step <= total_steps; step++) {
    auto force = volume.calculate_force();
    auto [avg_kinetic, total_energy] = volume.integrate_verle(force, dt);

    if (step % snapshot == 0) {
      std::cout << "Шаг " << step << ":\n";
      std::cout << "  Средняя кинетическая энергия: " << avg_kinetic << "\n";
      std::cout << "  Полная энергия: " << total_energy << "\n";
      for (int i = 0; i < 2; i++) {
        space.get_molecule(i).print_full_information();
      }
      // space.impulse_print();

      space.write_velocity(velocity_file);
      space.write_coord_abs(coord_abs_file);
    }
  }
  std::cout << "===== Симуляция завершена =====\n";

  return 0;
}
