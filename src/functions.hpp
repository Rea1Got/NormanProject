#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP
#include "molecule.hpp"
#include "volume.hpp"
#include <vector>

// std::vector<std::array<double, 3>> calculate_force(Volume volume) {
//   double radius_2 = 0;
//   double force_scalar = 0;
//   int number_of_molecules = volume.get_space().get_amount_of_molecules();
//   std::vector<std::array<double, 3>> force(number_of_molecules, {0, 0, 0});
//
//   for (int i = 0; i < number_of_molecules - 1; i++) {
//     std::array<double, 3> coord_i =
//         volume.get_space().get_molecule(i).get_coordinate();
//     for (int j = i + 1; j < number_of_molecules; j++) {
//       std::array<double, 3> coord_j =
//           volume.get_space().get_molecule(j).get_coordinate();
//       std::array<double, 3> dr = {0, 0, 0};
//
//       for (int k = 0; k < 3; k++) {
//         dr[k] = coord_i[k] - coord_j[k];
//         // periodic_border
//         dr[k] -=
//             volume.get_length(k) * std::round(dr[k] / volume.get_length(k));
//       }
//       radius_2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
//       force_scalar = volume.lennard_jones(radius_2);
//
//       for (int k = 0; k < 3; k++) {
//         force[i][k] += force_scalar * dr[k];
//         force[j][k] -= force_scalar * dr[k];
//       }
//     }
//   }
//
//   return force;
// }
#endif
