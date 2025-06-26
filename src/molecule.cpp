#include "molecule.hpp"
#include <iostream>

void Molecule::get_id() { std::cout << "Id: " << id << std::endl; }

void Molecule::get_coordinate() {
  std::cout << "Coordinates { x: " << coordinate[0] << "; y: " << coordinate[1]
            << "; z: " << coordinate[2] << " }" << std::endl;
}
void Molecule::get_velocity() {
  std::cout << "Velocity { v_x: " << velocity[0] << "; v_y: " << velocity[1]
            << "; v_z: " << velocity[2] << " }" << std::endl;
}
void Molecule::get_acceleration() {
  std::cout << "Acceleration { a_x: " << acceleration[0]
            << "; a_y: " << acceleration[1] << "; a_z: " << acceleration[2]
            << " }" << std::endl;
}

void Molecule::get_full_information() {
  get_id();
  get_coordinate();
  get_velocity();
  get_acceleration();
  std::cout << std::endl;
}
