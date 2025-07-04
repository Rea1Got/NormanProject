// #include "src/functions.hpp"
#include "src/molecule.hpp"
#include "src/space.hpp"
#include "src/volume.hpp"
#include <cstdlib>
#include <iostream>

int main() {
  const int seed = 2;
  const int num_molecules = 500;
  const double temperature = 1.0;
  const double dt = 0.001;
  const int total_steps = 100;
  const double r_cut = 2.5;
  const double sigma = 1;
  const double epsilon = 1;

  Volume volume(num_molecules, seed, temperature, r_cut, sigma, epsilon, 100,
                100, 100);

  std::cout << "===== Начало симуляции =====\n";
  std::cout << "Молекул: " << num_molecules << "\n";
  std::cout << "Шаг по времени: " << dt << "\n";
  std::cout << "Температура: " << temperature << "\n\n";

  std::cout << "Скорость на 0 шаге: "
            << volume.get_space().get_molecule(0).get_velocity()[2] << "\n";
  for (int step = 0; step <= total_steps; step++) {
    auto force = volume.calculate_force();

    auto [avg_kinetic, total_energy] = volume.integrate_verle(force, dt);

    if (step % 10 == 0) {
      std::cout << "Шаг " << step << ":\n";
      std::cout << "  Средняя кинетическая энергия: " << avg_kinetic << "\n";
      std::cout << "  Полная энергия: " << total_energy << "\n";

      for (int i = 0; i < 2; i++) {
        volume.get_space().get_molecule(i).print_full_information();
      }
      std::cout << "----------------------\n";
    }
  }

  std::cout << "===== Симуляция завершена =====\n";
  return 0;
}
