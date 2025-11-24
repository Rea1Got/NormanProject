#ifndef SPACE_HPP
#define SPACE_HPP 1

#include "../include/json.hpp"
#include "functions.hpp"
#include "molecule.hpp"
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

class Space {
private:
  std::vector<Molecule> molecules;
  double dt;

public:
  Space() = default;

  Space(int number_of_molecules, double mass = 1.0) {
    molecules.reserve(number_of_molecules);
    for (int i = 0; i < number_of_molecules; i++) {
      molecules.push_back(Molecule());
      molecules.back().set_id(i);
      molecules.back().set_mass(mass);
    }
    std::ifstream cfg_file("cfg/cfg.json");
    nlohmann::json cfg;
    cfg_file >> cfg;
    dt = cfg["dt"];
  }

  Space(std::vector<Molecule> init_molecules) : molecules(init_molecules) {};

  Molecule &get_molecule(int index) { return molecules[index]; }
  const Molecule &get_molecule(int index) const { return molecules[index]; }
  const double get_total_mass() const {
    double mass = 0.0;
    for (int i = 0; i < molecules.size(); i++) {
      mass += molecules[i].get_mass();
    }
    return mass;
  }

  int get_amount_of_molecules() { return molecules.size(); }
  void change_molecules(const std::vector<Molecule> &new_molecules) {
    molecules = new_molecules;
  };

  void rescale_velocity(double target_temperature) {
    int num_molecules = molecules.size();

    double current_kinetic_energy = 0.0;
    std::array<double, 3> total_momentum = {0, 0, 0};

    for (int i = 0; i < num_molecules; i++) {
      Molecule &mol = molecules[i];
      std::array<double, 3> vel = mol.get_velocity();

      for (int j = 0; j < 3; j++) {
        current_kinetic_energy +=
            0.5 * vel[j] * vel[j]; // E_kin = 1/2 * m * v^2 (m=1)
        total_momentum[j] += vel[j];
      }
    }

    for (int j = 0; j < 3; j++) {
      total_momentum[j] /= num_molecules;
    }

    // T = (2/3) * E_kin / N
    double current_temperature =
        (2.0 / 3.0) * current_kinetic_energy / num_molecules;

    double rescale_factor = 1.0;
    if (current_temperature > 0) {
      rescale_factor = std::sqrt(target_temperature / current_temperature);
    }

    for (int i = 0; i < num_molecules; i++) {
      Molecule &mol = molecules[i];
      std::array<double, 3> coordinates = mol.get_coordinate();
      std::array<double, 3> vel = mol.get_velocity();
      std::array<double, 3> coordinates_prev;

      for (int j = 0; j < 3; j++) {
        vel[j] = (vel[j] - total_momentum[j]) * rescale_factor;

        coordinates_prev[j] = coordinates[j] - vel[j] * dt;
      }

      mol.set_coordinate_prev(coordinates_prev);
      mol.set_velocity(vel);
    }
  }

  void remove_total_momentum(double temperature) {
    if (molecules.empty())
      return;

    std::array<double, 3> total_velocity = {0, 0, 0};
    int num_mol = molecules.size();

    for (int i = 0; i < num_mol; i++) {
      Molecule &mol = molecules[i];
      std::array<double, 3> vel = mol.get_velocity();
      for (int j = 0; j < 3; j++) {
        total_velocity[j] += vel[j];
      }
    }

    rescale_velocity(temperature);
  }

  void set_velocity(const std::string file) {
    std::vector<double> vel_data = read_input(file);
    int num_molecules = get_amount_of_molecules();

    for (int i = 0; i < num_molecules; i++) {
      std::array<double, 3> tmp = {vel_data[3 * i], vel_data[3 * i + 1],
                                   vel_data[3 * i + 2]};
      get_molecule(i).set_velocity(tmp);
    }
  }

  void set_coordinate(const std::string file) {
    std::vector<double> coord_data = read_input(file);
    std::array<double, 3> zeros = {0.0, 0.0, 0.0};
    int num_mol = get_amount_of_molecules();

    for (int i = 0; i < num_mol; i++) {
      Molecule &mol = get_molecule(i);
      std::array<double, 3> coord_prev = {
          coord_data[3 * i], coord_data[3 * i + 1], coord_data[3 * i + 2]};
      mol.set_coordinate_prev(coord_prev);
      mol.set_coordinate_abs(zeros);
    }

    for (int i = 0; i < num_mol; i++) {
      Molecule &mol = get_molecule(i);
      std::array<double, 3> coord = {coord_data[3 * num_mol + 3 * i],
                                     coord_data[3 * num_mol + 3 * i + 1],
                                     coord_data[3 * num_mol + 3 * i + 2]};
      mol.set_coordinate(coord);
    }
  }

  void write_velocity(std::string file) {
    // create folder if needed
    std::filesystem::create_directories(
        std::filesystem::path(file).parent_path());

    std::ofstream f_mol_vel(file, std::ios::app);
    for (int i = 0; i < get_amount_of_molecules(); i++) {
      for (int j = 0; j < 3; j++) {
        f_mol_vel << get_molecule(i).get_velocity()[j] << " ";
      }
    }
    f_mol_vel << "\n";
    f_mol_vel.close();
  }

  void write_coord(std::string file) {
    // create folder if needed
    std::filesystem::create_directories(
        std::filesystem::path(file).parent_path());

    std::ofstream f_coord(file, std::ios::app);
    for (int i = 0; i < get_amount_of_molecules(); i++) {
      for (int j = 0; j < 3; j++) {
        f_coord << get_molecule(i).get_coordinate()[j] << " ";
      }
    }
    f_coord << "\n";
    f_coord.close();
  }

  void write_coord_abs(std::string file) {
    // create folder if needed
    std::filesystem::create_directories(
        std::filesystem::path(file).parent_path());

    std::ofstream f_coord_abs(file, std::ios::app);
    for (int i = 0; i < get_amount_of_molecules(); i++) {
      for (int j = 0; j < 3; j++) {
        f_coord_abs << get_molecule(i).get_coordinate_abs()[j] << " ";
      }
    }
    f_coord_abs << "\n";
    f_coord_abs.close();
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
