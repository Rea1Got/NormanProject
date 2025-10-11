#include "include/json.hpp" //nlohmann
#include "src/simulation.hpp"
#include "src/volume.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

using json = nlohmann::json;

int main() {
  std::ifstream cfg_file("cfg/cfg.json");
  json cfg;
  cfg_file >> cfg;

  std::string velocity_file = cfg["velocity_file"];
  std::string coord_file = cfg["coord_file"];
  std::string coord_abs_file = cfg["coord_abs_file"];

  int seed = cfg["seed"];
  int num_molecules = cfg["num_molecules"];
  double temperature = cfg["temperature"];
  double dt = cfg["dt"];
  int total_steps = cfg["total_steps"];
  double r_cut = cfg["r_cut"];
  double sigma = cfg["sigma"];
  double epsilon = cfg["epsilon"];
  int snapshot = cfg["snapshot"];
  double density = cfg["density"];
  double mass_real = cfg["mass_real"];
  double sigma_real = cfg["sigma_real"];

  int answer;
  std::cout << "Init modeling or simulation from a snapshot? (0; 1) \n";
  std::cin >> answer;
  if (answer != 0 and answer != 1) {
    std::cout << "Invalid input.";
  }

  if (!answer) {
    Volume volume(num_molecules, seed, density, mass_real, sigma_real,
                  temperature, r_cut, sigma, epsilon);
    simulation(volume, answer, num_molecules, temperature, dt, total_steps,
               snapshot, coord_abs_file, coord_file, velocity_file);
  } else {
    std::string velocity_file_input = cfg["velocity_file_input"];
    std::string coord_file_input = cfg["coord_file_input"];
    Volume volume(num_molecules, seed, density, mass_real, sigma_real,
                  velocity_file_input, coord_file_input, temperature, r_cut,
                  sigma, epsilon);
    simulation(volume, answer, num_molecules, temperature, dt, total_steps,
               snapshot, coord_abs_file, coord_file, velocity_file);
  }

  return 0;
}
