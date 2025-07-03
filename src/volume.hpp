#ifndef VOLUME_HPP
#define VOLUME_HPP

#include "space.hpp"
#include <iostream>

class Volume {
private:
  double length_x = 1;
  double length_y = 1;
  double length_z = 1;
  double radius_cut_abs = 1;
  double epsilon = 1;
  double sigma = 1;
  double sigma_6 = 1;
  double sigma_12 = 1;
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
         double radius_cut_in_sigma = 1.0, double sigma = 1.0,
         double epsilon = 1.0, double x = 1.0, double y = 1.0, double z = 1.0)
      : length_x(x), length_y(y), length_z(z), epsilon(epsilon), sigma(sigma),
        space(number_of_molecules) {
    sigma_6 = std::pow(sigma, 6);
    sigma_12 = sigma_6 * sigma_6;
    radius_cut_abs = radius_cut_in_sigma * sigma;

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

  double lennard_jones(double r2) {
    double force = 0;
    if (r2 <= radius_cut_abs * radius_cut_abs && r2 > 1e-15) {
      double r2inv = 1.0 / r2;              // 1/r^2
      double r6inv = r2inv * r2inv * r2inv; // 1/r^6
      double r8inv = r6inv * r2inv;         // 1/r^8

      force = 48 * epsilon * r8inv * (sigma_12 * r6inv - 0.5 * sigma_6);
    }
    return force;
  }

  std::vector<std::array<double, 3>> calculate_force() {
    double radius_2 = 0;
    double force_scalar = 0;
    auto &space = get_space();
    int number_of_molecules = space.get_amount_of_molecules();
    std::vector<std::array<double, 3>> force(number_of_molecules, {0, 0, 0});

    for (int i = 0; i < number_of_molecules - 1; i++) {
      std::array<double, 3> coord_i = space.get_molecule(i).get_coordinate();
      for (int j = i + 1; j < number_of_molecules; j++) {
        std::array<double, 3> coord_j = space.get_molecule(j).get_coordinate();
        std::array<double, 3> dr = {0, 0, 0};

        for (int k = 0; k < 3; k++) {
          dr[k] = coord_i[k] - coord_j[k];
          // periodic_border
          dr[k] -= get_length(k) * std::round(dr[k] / get_length(k));
        }
        radius_2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        force_scalar = lennard_jones(radius_2);

        for (int k = 0; k < 3; k++) {
          force[i][k] += force_scalar * dr[k];
          force[j][k] -= force_scalar * dr[k];
        }
      }
    }

    return force;
  }

  int get_radius_cut() const { return radius_cut_abs; }
  double get_length_x() const { return length_x; }
  double get_length_y() const { return length_y; }
  double get_length_z() const { return length_z; }
  double get_length(int axis) const {
    if (axis == 0)
      return length_x;
    else if (axis == 1)
      return length_y;
    else
      return length_z;
  }

  void set_radius_cut(double new_radius) { radius_cut_abs = new_radius; }
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
