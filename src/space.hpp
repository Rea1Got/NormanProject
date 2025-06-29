#ifndef SPACE_HPP
#define SPACE_HPP

#include "molecule.hpp"
#include <cmath>
#include <vector>

class Space {
private:
  std::vector<Molecule> molecules;

public:
  Space(std::vector<Molecule> init_molecules) : molecules(init_molecules) {};

  Molecule &get_molecule(int index) { return molecules[index]; }
  const Molecule &get_molecule(int index) const { return molecules[index]; }

  int get_amount_of_molecules() { return molecules.size(); }
  void change_molecules(const std::vector<Molecule> &new_molecules) {
    molecules = new_molecules;
  };
};

#endif
