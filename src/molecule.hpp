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

  void get_id();
  void get_coordinate();
  void get_velocity();
  void get_acceleration();
  void get_full_information();
};

#endif
