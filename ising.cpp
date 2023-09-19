#include <iostream>
#include <random>
#include <time.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include "hdf5.h"

using std::cout;
using std::endl;

typedef std::vector<double> vd;
typedef std::vector<int> vi;
typedef std::vector<vi> vii;
typedef std::vector<vii> viii;

int random_spin(std::default_random_engine &rng) {
    std::uniform_real_distribution<double> dist(0, 1);
    double rand_num = dist(rng);
    // Breaking symmetry on purpose to support convergence
    if (rand_num <= .4) {
        return -1;
    }
    else {
        return 1;
    }
}

int random_site(std::default_random_engine &rng, int size) {
    std::uniform_real_distribution<double> dist(0, 1);
    double rand_num = dist(rng);
    int out = (int)(rand_num*(double)size);
    return out;
}

double random_p(std::default_random_engine &rng) {
    std::uniform_real_distribution<double> dist(0, 1);
    double p = dist(rng);
    return p;
}

vii init_lattice(std::default_random_engine & rng, int size) {
    vii lattice(size, vi(size, 0));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            lattice[i][j] = random_spin(rng);
        }
    }
    return lattice;
}

double energy(vii lattice, double J) {
    int size = lattice.size();
    double E = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            E += -J*lattice[i][j]*(lattice[(i+1+size)%size][j]+lattice[i][(j+1+size)%size]);
        }
    }
    return E;
}

double magnetization(vii lattice) {
    int size = lattice.size();
    int m_int = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            m_int += lattice[i][j];
        }
    }
    double m = (double)m_int;
    return m;
}

double dE(vii lattice, double J, int i, int j) {
    int size = lattice.size();
    int spin = lattice[i][j];
    int n1 = lattice[(i+1+size)%size][j];
    int n2 = lattice[(i-1+size)%(size)][j];
    int n3 = lattice[i][(j+1+size)%size];
    int n4 = lattice[i][(j-1+size)%size];
    double delta = -2*J*spin*(n1+n2+n3+n4);
    return delta;
}

struct MC_out {
    viii animation;
    vd E;
    vd M;
};

MC_out * MC(std::default_random_engine & rng, int size, double J, double T, int steps) {
    MC_out * out = new MC_out;
    // Parameters needed in Markov Chain
    int i, j;
    double delta, p, boltz, E_tmp, M_tmp;
    vd E_ar, M_ar;
    viii animation;
    vii lattice = init_lattice(rng, size);
    E_tmp = energy(lattice, J);
    M_tmp = magnetization(lattice);
    E_ar.push_back(E_tmp);
    M_ar.push_back(M_tmp);
    
    for (int k = 1; k < steps; k++) {
        i = random_site(rng, size);
        j = random_site(rng, size);
        // Flipping spin and checking change in energy
        lattice[i][j] *= -1;
        delta = dE(lattice, J, i, j);
        p = random_p(rng);
        boltz = exp(-delta/T);
        //cout << boltz << endl;
        if (p < boltz) {
            //cout << "Triggered" << endl;
            E_tmp += delta;
            E_ar.push_back(E_tmp);
            M_tmp += 2*lattice[i][j];
            M_ar.push_back(M_tmp);        
            animation.push_back(lattice);
        }
        else {
            lattice[i][j] *= -1;
        }            
    }
    out->animation = animation;
    out->E = E_ar;
    out->M = M_ar;
    return out;
}

int main() {
    // Seed for random number generator
    std::default_random_engine rng(time(NULL));
    int size = 100;
    double J = 1.;
    double T = .1;
    int steps = 100000;

    MC_out * out = MC(rng, size, J, T, steps);
    viii animation = out->animation;
    vd magnetization = out->M;
    vd energy = out->E;
    int n = magnetization.size();

    int nspins = round(pow(size, 2));
    for (int i = 0; i < n; i++) {
        magnetization[i] /= nspins;
        energy[i] /= nspins;
    }

    int n1 = animation.size();
    int n2 = animation[0].size();
    int n3 = animation[0][0].size();
    vi buffer;
    for (int i1 = 0; i1 < n1; i1++) {
        for (int i2 = 0; i2 < n2; i2++) {
            for (int i3 = 0; i3 < n3; i3++) {
                buffer.push_back(animation[i1][i2][i3]);
            }
        }
    }

    // Create .h5 file
    hid_t file_id = H5Fcreate("lattice.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Write lattice data
    hsize_t dims_animation[3] = {n1, n2, n3};
    hid_t dataspace_animation = H5Screate_simple(3, dims_animation, NULL);
    hid_t dataset_animation = H5Dcreate2(file_id, "/animation", H5T_IEEE_F64LE, dataspace_animation, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_animation, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer[0]);    
    H5Dclose(dataset_animation);
    H5Sclose(dataspace_animation);

    // Write magnetization data
    hsize_t dims_magnetization[1] = {n1};
    hid_t dataspace_magnetization = H5Screate_simple(1, dims_magnetization, NULL);
    hid_t dataset_magnetization = H5Dcreate2(file_id, "/magnetization", H5T_IEEE_F64LE, dataspace_magnetization, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_magnetization, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, magnetization.data());
    H5Dclose(dataset_magnetization);
    H5Sclose(dataspace_magnetization); 

    // Write energy data
    hsize_t dims_energy[1] = {n1};
    hid_t dataspace_energy = H5Screate_simple(1, dims_energy, NULL);
    hid_t dataset_energy = H5Dcreate2(file_id, "/energy", H5T_IEEE_F64LE, dataspace_energy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_energy, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, energy.data());
    H5Dclose(dataset_energy);
    H5Sclose(dataspace_energy);

    // Close .h5 file
    H5Fclose(file_id);

    delete out;

    return 0;
}


