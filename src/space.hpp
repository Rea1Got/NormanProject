#ifndef SPACE_HPP
#define SPACE_HPP

#include "molecule.hpp"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

class Space {
private:
  std::vector<Molecule> molecules;

public:
  Space() = default;

  Space(int number_of_molecules) {
    molecules.reserve(number_of_molecules);
    for (int i = 0; i < number_of_molecules; i++) {
      molecules.push_back(Molecule());
      molecules.back().set_id(i);
    }
  }

  Space(std::vector<Molecule> init_molecules) : molecules(init_molecules) {};

  Molecule &get_molecule(int index) { return molecules[index]; }
  const Molecule &get_molecule(int index) const { return molecules[index]; }

  int get_amount_of_molecules() { return molecules.size(); }
  void change_molecules(const std::vector<Molecule> &new_molecules) {
    molecules = new_molecules;
  };

  void remove_total_momentum(double temperature, double dt = 0.001) {
    if (molecules.empty())
      return;

    std::array<double, 3> total_velocity = {0, 0, 0};
    std::array<double, 3> total_velocity_2 = {0, 0, 0};
    int num_mol = molecules.size();

    for (int i = 0; i < num_mol; i++) {
      Molecule &mol = molecules[i];
      auto vel = mol.get_velocity();
      for (int j = 0; j < 3; j++) {
        total_velocity[j] += vel[j];
        total_velocity_2[j] += vel[j] * vel[j];
      }
    }

    for (int i = 0; i < 3; i++) {
      total_velocity[i] /= num_mol;
      total_velocity_2[i] /= num_mol;
    }

    std::array<double, 3> rescale_velocity{
        std::sqrt(3 * temperature / total_velocity_2[0]),
        std::sqrt(3 * temperature / total_velocity_2[1]),
        std::sqrt(3 * temperature / total_velocity_2[2])};

    for (int i = 0; i < num_mol; i++) {
      Molecule &mol = molecules[i];
      std::array<double, 3> coordinates_prev;
      std::array<double, 3> coordinates = mol.get_coordinate();
      std::array<double, 3> vel = mol.get_velocity();
      for (int j = 0; j < 3; j++) {
        vel[j] = (vel[j] - total_velocity[j]) * rescale_velocity[j];
        coordinates_prev[j] = coordinates[j] - vel[j] * dt;
      }

      mol.set_coordinate_prev(coordinates_prev);
      mol.set_velocity(vel);
    }
  }

  void write_velocity(std::string file_name) {
    std::ofstream f_mol_vel(file_name, std::ios::app);
    for (int i = 0; i < get_amount_of_molecules(); i++) {
      for (int j = 0; j < 3; j++) {
        f_mol_vel << get_molecule(i).get_velocity()[j] << " ";
      }
    }
    f_mol_vel << "\n";
    f_mol_vel.close();
  }

  void impulse_print() {
    std::array<double, 3> v = {0, 0, 0};
    for (int i = 0; i < get_amount_of_molecules(); i++) {
      std::array<double, 3> vel = get_molecule(i).get_velocity();
      for (int j = 0; j < 3; j++) {
        v[j] += vel[j];
      }
    }
    std::cout << "Импульс: "
              << std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) << "\n ";
    std::cout << "----------------------\n";
  }
};

#endif
