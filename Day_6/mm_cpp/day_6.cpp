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


Vector3d round_vector(Vector3d v)
{
    return v.array().round();
}


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


double lennard_jones_potential(double rij2)
{
    double sig_by_r6 = pow(1 / rij2, 3);
    double sig_by_r12 = sig_by_r6 * sig_by_r6;
    return 4.0 * (sig_by_r12 - sig_by_r6);
}


double get_particle_energy(size_t particle_index,
                           double box_length,
                           MatrixXd particle_coords,
                           double cutoff,
                           Vector3d displacement)
{
    const double cutoff2 = cutoff * cutoff;
    const size_t nparticles = particle_coords.rows();

    double e_total = 0.0;
    Vector3d i_position = particle_coords.row(particle_index);
    i_position += displacement;

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
