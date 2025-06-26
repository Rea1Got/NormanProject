#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <array>

class Molecule {
private:
  unsigned int id = 0;
  std::array<float, 3> coordinate = {0.0, 0.0, 0.0};
  std::array<float, 3> velocity = {0.0, 0.0, 0.0};
  std::array<float, 3> acceleration = {0.0, 0.0, 0.0};

public:
  Molecule() = default;
  Molecule(unsigned int id, const std::array<float, 3> &coord,
           const std::array<float, 3> &veloc, const std::array<float, 3> &accel)
      : id(id), coordinate(coord), velocity(veloc), acceleration(accel) {};

  void print_id();
  void print_coordinate();
  void print_velocity();
  void print_acceleration();
  void print_full_information();

  std::array<float, 3> get_coordinate() { return coordinate; }
  std::array<float, 3> get_velocity() { return velocity; }
  std::array<float, 3> get_acceleration() { return acceleration; }

  void change_coordinate(std::array<float, 3> new_coordinate) {
    coordinate = new_coordinate;
  };
  void change_velocity(std::array<float, 3> new_velocity) {
    velocity = new_velocity;
  };
  void change_acceleration(std::array<float, 3> new_acceleration) {
    acceleration = new_acceleration;
  };
};

#endif
