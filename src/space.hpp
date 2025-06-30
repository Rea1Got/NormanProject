#ifndef SPACE_HPP
#define SPACE_HPP

#include "molecule.hpp"
#include <cmath>
#include <cstdlib>
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

  // Space(int number_of_molecules, int seed, double MAX_VELOCITY = 1.0) {
  //   std::srand(seed);
  //   std::vector<Molecule> molecules_set;
  //   // normalization to unit volume
  //   for (int i = 0; i < number_of_molecules; i++) {
  //     molecules.push_back(Molecule());
  //     std::array<double, 3> coords = {
  //         static_cast<double>(std::rand()) / RAND_MAX,
  //         static_cast<double>(std::rand()) / RAND_MAX,
  //         static_cast<double>(std::rand()) / RAND_MAX};
  //     std::array<double, 3> velocity = {
  //         static_cast<double>(std::rand()) / RAND_MAX * MAX_VELOCITY -
  //             0.5 * MAX_VELOCITY,
  //         static_cast<double>(std::rand()) / RAND_MAX * MAX_VELOCITY -
  //             0.5 * MAX_VELOCITY,
  //         static_cast<double>(std::rand()) / RAND_MAX * MAX_VELOCITY -
  //             0.5 * MAX_VELOCITY};
  //     molecules[i].set_coordinate(coords);
  //     molecules[i].set_velocity(velocity);
  //     molecules[i].set_id(i);
  //   }
  // };

  Molecule &get_molecule(int index) { return molecules[index]; }
  const Molecule &get_molecule(int index) const { return molecules[index]; }

  int get_amount_of_molecules() { return molecules.size(); }
  void change_molecules(const std::vector<Molecule> &new_molecules) {
    molecules = new_molecules;
  };

  void remove_total_momentum(double temperature) {
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
        std::sqrt(3 * temperature * total_velocity_2[1]),
        std::sqrt(3 * temperature * total_velocity_2[2])};

    for (int i = 0; i < num_mol; i++) {
      Molecule &mol = molecules[i];
      auto vel = mol.get_velocity();
      for (int j = 0; j < 3; j++) {
        vel[j] = (vel[j] - total_velocity[j]) * rescale_velocity[j];
      }
      mol.set_velocity(vel);
    }
  }
};

#endif
