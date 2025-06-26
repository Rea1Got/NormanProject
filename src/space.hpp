#ifndef SPACE_HPP
#define SPACE_HPP

#include "molecule.hpp"
#include <vector>

class Space {
private:
  int lenght_x = 0;
  int lenght_y = 0;
  int lenght_z = 0;

public:
  Space(int x, int y, int z) : lenght_x(x), lenght_y(y), lenght_z(z) {};
  void print_space();
  void change_space(int, int, int);
};

class Volume : public Space {
private:
  std::vector<Molecule> molecules;

public:
  Volume(int x, int y, int z, std::vector<Molecule> init_molecules)
      : Space(x, y, z), molecules(init_molecules) {}
  void print_full_information();
};
#endif
