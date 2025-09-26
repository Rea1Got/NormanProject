#ifndef MOLECULE_HPP
#define MOLECULE_HPP 1

#include <array>

class Molecule {
private:
  unsigned int id = 0;
  double mass = 1;
  std::array<double, 3> coordinate = {0.0, 0.0, 0.0};
  std::array<double, 3> coordinate_prev = {0.0, 0.0, 0.0};
  std::array<double, 3> coordinate_abs = {0.0, 0.0, 0.0};
  std::array<double, 3> velocity = {0.0, 0.0, 0.0};

public:
  Molecule() = default;
  Molecule(unsigned int id, double mass, const std::array<double, 3> &coord,
           const std::array<double, 3> &coord_prev,
           const std::array<double, 3> &veloc)
      : id(id), mass(mass), coordinate(coord), coordinate_prev(coord_prev),
        coordinate_abs(coord), velocity(veloc) {};

  void print_id() const;
  void print_mass() const;
  void print_coordinate() const;
  void print_velocity() const;
  void print_full_information() const;

  int get_id() const { return id; }
  double get_mass() const { return mass; }
  const std::array<double, 3> &get_coordinate() const { return coordinate; }
  const std::array<double, 3> &get_coordinate_prev() const {
    return coordinate_prev;
  }
  const std::array<double, 3> &get_coordinate_abs() const {
    return coordinate_abs;
  }
  const std::array<double, 3> &get_velocity() const { return velocity; }

  void set_id(unsigned int new_id) { id = new_id; }
  void set_mass(double new_mass) { mass = new_mass; }
  void set_coordinate(std::array<double, 3> &new_coordinate) {
    coordinate = new_coordinate;
  };
  void set_coordinate_prev(std::array<double, 3> &new_coordinate) {
    coordinate_prev = new_coordinate;
  }
  void set_coordinate_abs(std::array<double, 3> &new_coordinate) {
    coordinate_abs = new_coordinate;
  }
  void set_velocity(std::array<double, 3> &new_velocity) {
    velocity = new_velocity;
  };
};

#endif
