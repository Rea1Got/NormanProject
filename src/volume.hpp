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
  double energy_c = 1;
  double mass_real = 1E-21;
  double sigma_real = 1E-10;
  Space space;

  void adjust_coordinate(double &coord, double length) {
    double original = coord;
    coord = std::fmod(coord, length);
    if (coord < 0) {
      coord += length;
    }
  }

  void border_periodic() {
    for (int i = 0; i < space.get_amount_of_molecules(); i++) {
      Molecule &current_molecule = space.get_molecule(i);
      std::array<double, 3> coord = current_molecule.get_coordinate();
      std::array<double, 3> coord_prev = current_molecule.get_coordinate_prev();
      adjust_coordinate(coord[0], length_x);
      adjust_coordinate(coord[1], length_y);
      adjust_coordinate(coord[2], length_z);
      adjust_coordinate(coord_prev[0], length_x);
      adjust_coordinate(coord_prev[1], length_y);
      adjust_coordinate(coord_prev[2], length_z);

      current_molecule.set_coordinate(coord);
      current_molecule.set_coordinate_prev(coord_prev);
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
          std::array<double, 3> coord_abs_init = {0.0, 0.0, 0.0};
          curr_mol.set_coordinate(coord);
          curr_mol.set_coordinate_abs(coord_abs_init);
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

  void init_density(double density_SI) {
    double total_mass_SI = space.get_amount_of_molecules() * mass_real;
    double volume_SI = total_mass_SI / density_SI;
    double length_SI = std::pow(volume_SI, 1.0 / 3.0);
    double length = length_SI / sigma_real;

    length_x = length;
    length_y = length;
    length_z = length;

    validate_parameters();
  }

  void validate_parameters() {
    double V_star = length_x * length_y * length_z;
    double rho_star = space.get_amount_of_molecules() / V_star;
    if (rho_star < 0.7 || rho_star > 1.2) {
      std::cerr << "WARNING: Unusual density: " << rho_star
                << " (expected 0.8-0.9 for liquid Ar)\n";
    }
  }

public:
  // Volume(int num_molecules, int seed, double temperature = 1.0,
  //        double radius_cut_in_sigma = 1.0, double sigma = 1.0,
  //        double epsilon = 1.0, double x = 1.0, double y = 1.0, double z
  //        = 1.0)
  //     : epsilon(epsilon), sigma(sigma), space(num_molecules) {
  //   sigma_6 = std::pow(sigma, 6);
  //   sigma_12 = sigma_6 * sigma_6;
  //   radius_cut_abs = radius_cut_in_sigma * sigma;
  //   energy_c = potential_energy(std::pow(radius_cut_abs, -2));
  //   length_x = x * sigma;
  //   length_y = y * sigma;
  //   length_z = z * sigma;
  //
  //   init_cubic_grid();
  //   init_velocity(seed);
  //   border_periodic();
  //
  //   space.remove_total_momentum(temperature);
  // }
  Volume(int num_molecules, int seed, double density_SI, double mass_real,
         double sigma_real, double temperature = 1.0,
         double radius_cut_in_sigma = 2.5, double sigma = 1.0,
         double epsilon = 1.0)
      : epsilon(epsilon), sigma(sigma), space(num_molecules),
        sigma_real(sigma_real), mass_real(mass_real) {
    sigma_6 = std::pow(sigma, 6);
    sigma_12 = sigma_6 * sigma_6;
    radius_cut_abs = radius_cut_in_sigma * sigma;
    energy_c = potential_energy(std::pow(radius_cut_abs, -2));

    init_density(density_SI);
    validate_parameters();

    init_cubic_grid();
    init_velocity(seed);
    border_periodic();

    space.remove_total_momentum(temperature);
  }

  Volume(int num_molecules, int seed, double density_SI, double mass_real,
         double sigma_real, std::string f_vel, std::string f_coord,
         double temperature = 1.0, double radius_cut_in_sigma = 2.5,
         double sigma = 1.0, double epsilon = 1.0)
      : epsilon(epsilon), sigma(sigma), space(num_molecules),
        sigma_real(sigma_real), mass_real(mass_real) {
    sigma_6 = std::pow(sigma, 6);
    sigma_12 = sigma_6 * sigma_6;
    radius_cut_abs = radius_cut_in_sigma * sigma;
    energy_c = potential_energy(std::pow(radius_cut_abs, -2));
    //
    // init_density(density_SI);
    // validate_parameters();
    //
    // init_cubic_grid();
    // init_velocity(seed);
    // border_periodic();
    //
    // space.remove_total_momentum(temperature);

    space.set_velocity(f_vel);
    space.set_coordinate(f_coord);
    init_density(density_SI);
    space.remove_total_momentum(temperature);
  }

  double potential_energy(double r2inv) {
    double r6inv = r2inv * r2inv * r2inv;
    return 4 * epsilon * r6inv * (sigma_12 * r6inv - sigma_6);
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
    double radius_2 = 0.0;
    double force_scalar = 0.0;
    double pot_energy = 0.0;
    Space &space = get_space();
    int num_molecules = space.get_amount_of_molecules();
    std::vector<std::array<double, 3>> force(num_molecules, {0, 0, 0});

    for (int i = 0; i < num_molecules - 1; i++) {
      std::array<double, 3> coord_i = space.get_molecule(i).get_coordinate();
      for (int j = i + 1; j < num_molecules; j++) {
        std::array<double, 3> coord_j = space.get_molecule(j).get_coordinate();
        std::array<double, 3> dr = {0, 0, 0};

        for (int k = 0; k < 3; k++) {
          dr[k] = coord_i[k] - coord_j[k];
          // periodic_border rewrite to adjust_coordinate
          dr[k] -= get_length(k) * std::round(dr[k] / get_length(k));
        }
        radius_2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

        pot_energy += potential_energy(1.0 / radius_2) - energy_c;

        force_scalar = lennard_jones(radius_2);
        for (int k = 0; k < 3; k++) {
          force[i][k] += force_scalar * dr[k];
          force[j][k] -= force_scalar * dr[k];
        }
      }
    }
    //  bruh
    //  need to separate calculation of dr[], calculation of forces and energy
    force.push_back({pot_energy, 0.0, 0.0});

    return force;
  }

  std::pair<double, double>
  integrate_verle(std::vector<std::array<double, 3>> &force,
                  double dt = 0.001) {
    std::pair<double, double> result; // (moment temp and energy) on 1 molecule
    double potential_energy = force.back()[0];
    double kinetic_energy = 0.0;
    std::array<double, 3> coordinates_new;
    std::array<double, 3> coordinates;
    std::array<double, 3> coordinates_prev;
    std::array<double, 3> coordinates_abs;
    std::array<double, 3> coordinates_abs_new;
    std::array<double, 3> velocity;
    std::array<double, 3> velocity_curr;
    double mass = 1.0;
    double acceleration = 0.0;

    unsigned int num_molecules = space.get_amount_of_molecules();
    for (int i = 0; i < num_molecules; i++) {
      Molecule &molecule = space.get_molecule(i);
      mass = molecule.get_mass();
      coordinates = molecule.get_coordinate();
      coordinates_prev = molecule.get_coordinate_prev();
      coordinates_abs = molecule.get_coordinate_abs();
      velocity_curr = molecule.get_velocity();

      for (int j = 0; j < 3; j++) {
        acceleration = force[i][j] / mass;
        coordinates_new[j] =
            2 * coordinates[j] - coordinates_prev[j] + acceleration * dt * dt;

        coordinates_abs_new[j] =
            coordinates_abs[j] + coordinates_new[j] - coordinates_prev[j];

        double volume_len = get_length(j);
        adjust_coordinate(coordinates_new[j], volume_len);

        // Minimum Image Convention
        double delta_r = coordinates_new[j] - coordinates_prev[j];
        if (delta_r > 0.5 * volume_len) {
          delta_r -= volume_len;
        } else if (delta_r < -0.5 * volume_len) {
          delta_r += volume_len;
        }
        velocity[j] = delta_r / (2 * dt);

        kinetic_energy += 0.5 * mass * velocity[j] * velocity[j];
      }
      molecule.set_coordinate_prev(coordinates);
      molecule.set_coordinate(coordinates_new);
      molecule.set_coordinate_abs(coordinates_abs_new);
      molecule.set_velocity(velocity);
    }

    result.first = kinetic_energy / num_molecules;
    result.second = (potential_energy + kinetic_energy) / num_molecules;

    return result;
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
