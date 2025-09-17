#include "molecule.hpp"
#include <iostream>

void Molecule::print_id() const { std::cout << "Id: " << id << std::endl; }
void Molecule::print_mass() const {
  std::cout << "Mass: " << mass << std::endl;
}
void Molecule::print_coordinate() const {
  std::cout << "Coordinates { " << coordinate[0] << "; " << coordinate[1]
            << "; " << coordinate[2] << " }" << std::endl;
}
void Molecule::print_velocity() const {
  std::cout << "Velocity { " << velocity[0] << "; " << velocity[1] << "; "
            << velocity[2] << " }" << std::endl;
}

void Molecule::print_full_information() const {
  print_id();
  print_coordinate();
  print_velocity();
  std::cout << "Coord_abs { " << coordinate_abs[0] << "; " << coordinate_abs[1]
            << "; " << coordinate_abs[2] << " }\n";
  std::cout << std::endl;
}
