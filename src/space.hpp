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
  Space(int number_of_molecules, int seed, double MAX_VELOCITY = 1.0) {
    std::srand(seed);
    std::vector<Molecule> molecules_set;
    // normalization to unit volume; if needed bigger
    for (int i = 0; i < number_of_molecules; i++) {
      molecules.push_back(Molecule());
      std::array<double, 3> coords = {
          static_cast<double>(std::rand()) / RAND_MAX,
          static_cast<double>(std::rand()) / RAND_MAX,
          static_cast<double>(std::rand()) / RAND_MAX};
      std::array<double, 3> velocity = {
          static_cast<double>(std::rand()) / RAND_MAX * MAX_VELOCITY -
              0.5 * MAX_VELOCITY,
          static_cast<double>(std::rand()) / RAND_MAX * MAX_VELOCITY -
              0.5 * MAX_VELOCITY,
          static_cast<double>(std::rand()) / RAND_MAX * MAX_VELOCITY -
              0.5 * MAX_VELOCITY};
      molecules[i].set_coordinate(coords);
      molecules[i].set_velocity(velocity);
      molecules[i].set_id(i);
    }
  };

  Space(std::vector<Molecule> init_molecules) : molecules(init_molecules) {};

  Molecule &get_molecule(int index) { return molecules[index]; }
  const Molecule &get_molecule(int index) const { return molecules[index]; }

  int get_amount_of_molecules() { return molecules.size(); }
  void change_molecules(const std::vector<Molecule> &new_molecules) {
    molecules = new_molecules;
  };
};

#endif
