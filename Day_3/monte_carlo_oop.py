import numpy as np


# Molecule
class Molecule:

    def __init__(self, name=None, charge=None, symbols=None, coordinates=None):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.coordinates = coordinates

        self.energy_array = None
        self.total_pair_energy = None  # TODO


    def __str__(self):
        return 'name: ' + str(self.name) + '\nsymbols:' + str(self.symbols)


class Box:

    def __init__(self, box_length, molecules=None):
        self.box_length = box_length
        self.particles = molecules
        self.cutoff = box_length / 2.0
        self.cutoff2 = np.power(self.cutoff, 2)


    @staticmethod
    def lennard_jones_potential(rij2):
        """ Lennard Jones potential implementation"""
        sig_by_r6 = np.power(1 / rij2, 3)
        sig_by_r12 = np.power(sig_by_r6, 2)
        return 4.0 * (sig_by_r12 - sig_by_r6)


    def minimum_image_distance(self, r_i, r_j):
        """Minimum image distance implementation"""
        rij = r_i - r_j
        rij = rij - self.box_length * np.round(rij / self.box_length)
        rij2 = np.dot(rij, rij)
        return rij2

    def calc_total_potential_energy(self):
        """Computation of the total system energy"""
        e_total = 0.0

        for i_particle in range(len(self.particles)):
            for j_particle in range(i_particle):

                r_i = self.particles[i_particle].coordinates
                r_j = self.particles[j_particle].coordinates
                rij2 = self.minimum_image_distance(r_i, r_j)
                if rij2 < self.cutoff2:
                    e_pair = self.lennard_jones_potential(rij2)
                    e_total += e_pair

        self.total_pair_energy = e_total

    def calc_tail_correction(self):
        """Computation of the energy tail correction"""

        volume = np.power(self.box_length, 3)
        sig_by_cutoff3 = np.power(1.0 / self.cutoff, 3)
        sig_by_cutoff9 = np.power(sig_by_cutoff3, 3)
        e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3
        e_correction *= 8.0 / 9.0 * np.pi * len(self.particles) / volume * len(self.particles)

        self.tail_correction = e_correction


    def get_molecule_energy(self, i_particle):

        e_total = 0.0
        i_position = self.particles[i_particle].coordinates

        for j_particle in range(len(self.particles)):
            if i_particle != j_particle:

                j_position = self.particles[j_particle].coordinates
                rij2 = self.minimum_image_distance(i_position, j_position)

                if rij2 < self.cutoff2:
                    e_pair = self.lennard_jones_potential(rij2)
                    e_total += e_pair

        return e_total

    def move_particle(self, max_displacement):

        i_particle = np.random.randint(len(self.particles))
        random_displacement = (2.0 * np.random.rand(3) - 1.0)* max_displacement
        old_position = self.particles[i_particle].coordinates.copy()
        old_energy = self.get_molecule_energy(i_particle)
        self.particles[i_particle].coordinates += random_displacement
        new_energy = self.get_molecule_energy(i_particle)

        return new_energy - old_energy, i_particle, old_position


    def generate_initial_state(self, element, num_particles, method='lattice'):
        """Generate initial state"""
        if method is 'random':
            # Randomly placing particles in a box
            coordinates = (0.5 - np.random.rand(num_particles, 3)) * self.box_length

        elif method is 'lattice':
            # Adding particles in a lattice. There might be many more
            # ways to do this!
            spacing = int(np.cbrt(num_particles) + 1)
            x_vector = np.linspace(0.0, self.box_length, spacing)
            y_vector = np.linspace(0.0, self.box_length, spacing)
            z_vector = np.linspace(0.0, self.box_length, spacing)
            grid  = np.meshgrid(x_vector, y_vector, z_vector)
            stack = np.vstack(grid)
            coordinates = stack.reshape(3, -1).T
            excess = len(coordinates) - num_particles
            coordinates = coordinates[:-excess]
            coordinates *= 0.95

        elif method is 'file':
            # Reading a reference configuration from NIST
            coordinates = np.loadtxt("lj_sample_config_periodic1.txt", skiprows=2, usecols=(1,2,3))

        self.particles = []
        print('coordinates: ', coordinates)
        for coordinate in coordinates:
            mol = Molecule(name=element, coordinates=coordinate)
            print('mol: ', mol)
            self.particles.append(mol)


    def restore_system(self, i_particle, old_position):
        self.particles[i_particle].coordinates = old_position

class OrthorhombicBox(Box):

    def __init__(self, width, *args, **kwargs):
        self.width = width
        super().__init__(*args, **kwargs)

    def move_particle(self, max_displacement):
        """Override parent"""
        pass

class TriclinicBox(Box):

    def __init__(self, width, angle, *args, **kwargs):
        self.width = width
        self.angle = angle
        super().__init__(*args, **kwargs)

    def move_particle(self, max_displacement):
        """Override parent"""
        pass

# --------------------------------------------------------------------------

def accept_or_reject(delta_e, beta):

    if delta_e <= 0.0:
        accept = True
    else:
        random_number = np.random.rand(1)
        p_acc = np.exp(-beta*delta_e)
        if random_number < p_acc:
            accept = True
        else:
            accept = False

    return accept


def adjust_displacement(freq, i_step, n_trials, n_accept, max_displacement):

    if np.mod(i_step + 1, freq) == 0:
        acc_rate = float(n_accept) / float(n_trials)
        if (acc_rate < 0.380):
            max_displacement *= 0.8
        elif (acc_rate > 0.42):
            max_displacement *= 1.2
        n_trials = 0
        n_accept = 0
    return max_displacement


def update_output_file(traj, i_step, freq, particles):

    if np.mod(i_step + 1, freq) == 0:
        traj.write(str(len(particles)) + '\n\n')
        for mol in particles:
            traj.write("%s %10.5f %10.5f %10.5f \n" % (mol.name, mol.coordinates[0],
                                                       mol.coordinates[1], mol.coordinates[2]))


# ---------------
# Parameter setup
# ---------------

reduced_density = 0.9
reduced_temperature = 0.9
num_particles = 100
max_displacement = 0.1
n_steps = 50000
freq = 1000
tune_displacement = True

# -------------------------
# Simulation initialization
# -------------------------

beta = 1 / reduced_temperature
_box_length = np.cbrt(num_particles / reduced_density)

box = Box(_box_length)
box.generate_initial_state(element='C', num_particles=num_particles, method='lattice')
box.calc_total_potential_energy()
box.calc_tail_correction()
box.energy_array = np.zeros(n_steps)


n_trials = 0
n_accept = 0


# --------------------------------
# Metropolis Monte Carlo algorithm
# --------------------------------

traj = open('traj.xyz', 'w')

for i_step in range(n_steps):

    n_trials += 1

    delta_e, i_particle, old_position = box.move_particle(max_displacement)

    accept = accept_or_reject(delta_e, beta)

    if accept:
        box.total_pair_energy += delta_e
        n_accept += 1
    else:
        box.restore_system(i_particle, old_position)

    if tune_displacement:
        max_displacement = adjust_displacement(freq, i_step, n_trials, n_accept, max_displacement)

    box.energy_array[i_step] = (box.total_pair_energy + box.tail_correction) / num_particles

    update_output_file(traj, i_step, freq, box.particles)

    if np.mod(i_step + 1, 1000) == 0:
        print (i_step + 1, box.energy_array[i_step])

traj.close()