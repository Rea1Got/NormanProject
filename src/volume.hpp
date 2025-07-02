#ifndef VOLUME_HPP
#define VOLUME_HPP

#include "space.hpp"
#include <iostream>

class Volume {
private:
  double length_x = 1;
  double length_y = 1;
  double length_z = 1;
  Space space;

  void adjust_coordinate(double &coord, float length) {
    coord = std::fmod(coord, length);
    if (coord < 0) {
      coord += length;
    }
  }

  void border_periodic() {
    for (int i = 0; i < space.get_amount_of_molecules(); i++) {
      auto &current_molecule = space.get_molecule(i); // border-periodic
      std::array<double, 3> coord = current_molecule.get_coordinate();
      adjust_coordinate(coord[0], length_x);
      adjust_coordinate(coord[1], length_y);
      adjust_coordinate(coord[2], length_z);
      current_molecule.set_coordinate(coord);
    }
  }

  void init_cubic_grid() {
    int n =
        static_cast<int>(std::ceil(std::cbrt(space.get_amount_of_molecules())));

    if (n == 0)
      n = 1;

    double dx = length_x / n;
    double dy = length_y / n;
    double dz = length_z / n;

    int count = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
          if (count >= space.get_amount_of_molecules())
            break;

          auto &curr_mol = space.get_molecule(count);
          std::array<double, 3> coord = {i * dx, j * dy, k * dz};
          curr_mol.set_coordinate(coord);

          count++;
        }
      }
    }
  }

  void init_velocity(int seed) {
    std::srand(seed);
    for (int i = 0; i < space.get_amount_of_molecules(); i++) {
      std::array<double, 3> velocity = {
          static_cast<double>(std::rand()) / RAND_MAX - 0.5,
          static_cast<double>(std::rand()) / RAND_MAX - 0.5,
          static_cast<double>(std::rand()) / RAND_MAX - 0.5};
      space.get_molecule(i).set_velocity(velocity);
    }
  }

public:
  Volume(int number_of_molecules, int seed, double temperature = 1.0,
         double x = 1.0, double y = 1.0, double z = 1.0)
      : length_x(x), length_y(y), length_z(z), space(number_of_molecules) {
    init_cubic_grid();
    init_velocity(seed);
    border_periodic();
    space.remove_total_momentum(temperature);
  }

  Volume(Space &space, double temperature = 1.0, double x = 1.0, double y = 1.0,
         double z = 1.0)
      : length_x(x), length_y(y), length_z(z), space(space) {
    space.remove_total_momentum(temperature);
  }

  void set_volume(double x, float y, float z) {
    length_x = x;
    length_y = y;
    length_z = z;
    border_periodic();
  }

  void print_volume_specs() {
    std::cout << "Current space { x: " << length_x << " y: " << length_y
              << " z: " << length_z << " }" << std::endl;
    std::cout << "Num of molecules:  " << space.get_amount_of_molecules()
              << std::endl;
  }
  Space &get_space() { return space; }
  const Space &get_space() const { return space; }
};

#endif
