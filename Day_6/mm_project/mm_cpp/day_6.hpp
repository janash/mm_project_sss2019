#pragma once

#include <Eigen/Dense>


/*! \brief Round all elements of a vector to integers
 *
 * The return type is still contains doubles, however
 * it should contain integer values.
 */
Eigen::Vector3d round_vector(Eigen::Vector3d v);


/*! \brief Calculate the distance between two points using the minimum image convention.
 *
 * The minimum image convention is for calculating distances in boxes with periodic
 * boundary conditions.
 *
 * \param [in] r_i The coordinates of the first point
 * \param [in] r_j The coordinates of the second point
 * \param [in] box_length The length of each side of the box
 *
 * \return The square distance between the two particles
 */
double minimum_image_distance(Eigen::Vector3d r_i, Eigen::Vector3d r_j, double box_length);


/*! \brief Calculate the Lennard Jones pairwise potential between two particles based on
 *         a separation distance.
 *
 * \param [in] rij2 The squared distance between two particles.
 *
 * \return The Lennard Jones interaction energy for two particles.
 */
double lennard_jones_potential(double rij2);


/*! \brief Calculate the interaction energy of a particle with its environment
 *         (all other particles in the system)
 *
 * \param [in] particle_index Index of the particle for which to calculate the energy
 * \param [in] box_length The length of each side of the box
 * \param [in] particle_coords Coordinates for all particles of the system
 * \param [in] cutoff The simulation cutoff. Beyond this distance, interactions are not calculated.
 * \param [in] displacement Displacement to apply to the particles
 *
 * \return The pairwise interaction energy of \p particle_index with all other particles
 *         in the system
 */
double get_particle_energy(size_t particle_index,
                           double box_length,
                           Eigen::MatrixXd particle_coords,
                           double cutoff,
                           Eigen::Vector3d displacement);


/*! \brief Calculate the total potential energy of the system using a pairwise
 *         Lennard Jones potential
 *
 * \param [in] box_length The length of each side of the box
 * \param [in] particle_coords Coordinates of all of the particles
 * \param [in] cutoff The simulation cutoff. Beyond this distance, interactions are not calculated.
 *
 * \return The total pairwise potential energy of the system
 */
double calculate_total_pair_energy(double box_length, Eigen::MatrixXd particle_coords, double cutoff);
