#include <cmath>
#include <fstream>
#include <string>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <Eigen/Dense>

// We will be using these a lot so for convenience
// make it so we don't have to put Eigen:: in front of them.
using Eigen::MatrixXd;
using Eigen::Vector3d;


/*! \brief Round all elements of a vector to integers
 *
 * The return type still contains doubles, however
 * it should contain integer values.
 */
Vector3d round_vector(Vector3d v)
{
    return v.array().round();
}


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
double minimum_image_distance(Vector3d r_i, Vector3d r_j, double box_length)
{
    Vector3d rij = r_i - r_j;

    // In python, this is: self.box_length * np.round(rij / self.box_length)
    Vector3d shift = rij / box_length;
    shift = round_vector(shift);
    shift = shift * box_length;

    rij = rij - shift;
    return rij.dot(rij);
}


/*! \brief Calculate the Lennard Jones pairwise potential between two particles based on
 *         a separation distance.
 *
 * \param [in] rij2 The squared distance between two particles.
 *
 * \return The Lennard Jones interaction energy for two particles.
 */
double lennard_jones_potential(double rij2)
{
    double sig_by_r6 = pow(1.0 / rij2, 3);
    double sig_by_r12 = sig_by_r6 * sig_by_r6;
    return 4.0 * (sig_by_r12 - sig_by_r6);
}


/*! \brief Calculate the interaction energy of a particle with its environment
 *         (all other particles in the system)
 *
 * \param [in] particle_index Index of the particle for which to calculate the energy
 * \param [in] box_length The length of each side of the box
 * \param [in] particle_coords Coordinates for all particles of the system
 * \param [in] cutoff The simulation cutoff. Beyond this distance, interactions are not calculated.
 *
 * \return The pairwise interaction energy of \p particle_index with all other particles
 *         in the system
 */
double get_particle_energy(size_t particle_index,
                           double box_length,
                           MatrixXd particle_coords,
                           double cutoff)
{
    const double cutoff2 = cutoff * cutoff;
    const size_t nparticles = particle_coords.rows();

    double e_total = 0.0;
    Vector3d i_position = particle_coords.row(particle_index);

    for(size_t j_particle = 0; j_particle < nparticles; ++j_particle)
    {
        if(particle_index == j_particle)
            continue;

        Vector3d j_position = particle_coords.row(j_particle);
        double rij2 = minimum_image_distance(i_position, j_position, box_length);

        if(rij2 < cutoff2)
        {
            double e_pair = lennard_jones_potential(rij2);
            e_total += e_pair;
        }
    }

    return e_total;
}


/*! \brief Calculate the total potential energy of the system using a pairwise
 *         Lennard Jones potential
 *
 * \param [in] box_length The length of each side of the box
 * \param [in] particle_coords Coordinates of all of the particles
 * \param [in] cutoff The simulation cutoff. Beyond this distance, interactions are not calculated.
 *
 * \return The total pairwise potential energy of the system
 */
double calculate_total_pair_energy(double box_length, MatrixXd particle_coords, double cutoff)
{
    const double cutoff2 = cutoff * cutoff;
    const size_t nparticles = particle_coords.rows();

    double e_total = 0.0;

    for(size_t i_particle = 0; i_particle < nparticles; ++i_particle)
    {
        Vector3d i_position = particle_coords.row(i_particle);
        for(size_t j_particle = 0; j_particle < i_particle; ++j_particle)
        {
            Vector3d j_position = particle_coords.row(j_particle);
            double rij2 = minimum_image_distance(i_position, j_position, box_length);

            if(rij2 < cutoff2)
            {
                double e_pair = lennard_jones_potential(rij2);
                e_total += e_pair;
            }
        }
    }

    return e_total;
}


/*! \brief Read a file containing data about a simulation
 *
 * \param [in] file_path Path to the file to read
 *
 * \return A vector representing the lengths of the sides of the box and a matrix
 *         containing the coordinates of all the particles.
 */
std::pair<Vector3d, MatrixXd> read_file(const std::string & file_path)
{
    std::ifstream data_file(file_path);

    if(! data_file.is_open() )
    {
        throw std::runtime_error("Error opening input file. Are you sure it exists?");
    }

    size_t nparticles;
    Vector3d box_length;

    data_file >> nparticles >> box_length(0) >> box_length(1) >> box_length(2);
    MatrixXd coords(nparticles, 3);

    size_t particle_number; // Dummy variable for reading in the number
                            // It is ignored

    for(size_t i = 0; i < nparticles; ++i)
        data_file >> particle_number >> coords(i, 0) >> coords(i, 1) >> coords(i, 2);

    return std::make_pair(box_length, coords);
}


int main(void)
{
    // Some constants
    const double simulation_cutoff = 3.0;

    auto file_data = read_file("../nist_sample_config1.txt");
    Vector3d box_length = file_data.first;
    MatrixXd coords = file_data.second;

    std::cout << "Read in " << coords.rows() << " particle coordinates\n";
    std::cout << "Box length: " << box_length[0] << " , " << box_length[1] << " , " << box_length[2] << "\n";

    // NOTE: We only use the first entry for box_length (assuming cubic box)
    double box_side = box_length[0];
    double particle_energy = get_particle_energy(0, box_side, coords, simulation_cutoff);
    double total_pair_energy = calculate_total_pair_energy(box_side, coords, simulation_cutoff);
    std::cout << "Particle energy: " << particle_energy << "\n";
    std::cout << "Total pair energy: " << total_pair_energy << "\n";

    return 0;
}
