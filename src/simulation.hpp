#ifndef SIMULATION_HPP
#define SIMULATION_HPP 1

#include "volume.hpp"
#include <iostream>

void simulation(Volume &volume, int answer, int num_molecules,
                double temperature, double dt, double total_steps, int snapshot,
                std::string coord_abs_file, std::string coord_file,
                std::string velocity_file) {
  Space &space = volume.get_space();

  std::cerr << "===== Начало симуляции =====\n";
  std::cerr << "Молекул: " << num_molecules << "\n";
  std::cerr << "Шаг по времени: " << dt << "\n";
  std::cerr << "Температура: " << temperature << "\n";
  std::cerr << "Длина объема: " << volume.get_length_x() << "\n\n";

  for (int step = 0; step < total_steps; step++) {
    if (!answer) {                  // init_modeling
      if (step > total_steps - 3) { // last 2 coords
        space.write_coord(coord_file);
      }
      if (step == total_steps - 1) { // last velocity
        space.write_velocity(velocity_file);
      }
    }
    if (answer) { // simulation starting from init modeling
      space.write_coord_abs(coord_abs_file);
      if (step % snapshot == 0) {
        space.write_velocity(velocity_file);
        space.write_coord(coord_file);
      }
    }

    auto force = volume.calculate_force();
    auto [avg_kinetic, total_energy] = volume.integrate_verle(force, dt);
    if (step % 1000 == 0) {
      std::cerr << "Шаг " << step << ":\n";
      std::cerr << "  Средняя кинетическая энергия: " << avg_kinetic << "\n";
      std::cerr << "  Полная энергия: " << total_energy << "\n\n";
    }
  }
  std::cerr << "===== Симуляция завершена =====\n";
}

#endif
