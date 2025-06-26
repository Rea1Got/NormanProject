#include "molecule.hpp"
#include <iostream>

void Molecule::print_id() { std::cout << "Id: " << id << std::endl; }

void Molecule::print_coordinate() {
  std::cout << "Coordinates { " << coordinate[0] << "; " << coordinate[1]
            << "; " << coordinate[2] << " }" << std::endl;
}
void Molecule::print_velocity() {
  std::cout << "Velocity { " << velocity[0] << "; " << velocity[1] << "; "
            << velocity[2] << " }" << std::endl;
}
void Molecule::print_acceleration() {
  std::cout << "Acceleration { " << acceleration[0] << "; " << acceleration[1]
            << "; " << acceleration[2] << " }" << std::endl;
}

void Molecule::print_full_information() {
  print_id();
  print_coordinate();
  print_velocity();
  print_acceleration();
  std::cout << std::endl;
}
