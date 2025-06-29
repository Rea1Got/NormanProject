#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <array>

class Molecule {
private:
  unsigned int id = 0;
  double mass = 1;
  std::array<double, 3> coordinate = {0.0, 0.0, 0.0};
  std::array<double, 3> velocity = {0.0, 0.0, 0.0};
  std::array<double, 3> acceleration = {0.0, 0.0, 0.0};

public:
  Molecule() = default;
  Molecule(unsigned int id, double mass, const std::array<double, 3> &coord,
           const std::array<double, 3> &veloc,
           const std::array<double, 3> &accel)
      : id(id), mass(mass), coordinate(coord), velocity(veloc),
        acceleration(accel) {};

  void print_id() const;
  void print_mass() const;
  void print_coordinate() const;
  void print_velocity() const;
  void print_acceleration() const;
  void print_full_information() const;

  int get_id() const { return id; }
  double get_mass() const { return mass; }
  const std::array<double, 3> &get_coordinate() const { return coordinate; }
  const std::array<double, 3> &get_velocity() const { return velocity; }
  const std::array<double, 3> &get_acceleration() const { return acceleration; }

  void set_id(unsigned int new_id) { id = new_id; }
  void set_mass(double new_mass) { mass = new_mass; }
  void set_coordinate(std::array<double, 3> &new_coordinate) {
    coordinate = new_coordinate;
  };
  void set_velocity(std::array<double, 3> &new_velocity) {
    velocity = new_velocity;
  };
  void set_acceleration(std::array<double, 3> &new_acceleration) {
    acceleration = new_acceleration;
  };
};

#endif
