#ifndef PARTICLE_H
#define PARTICLE_H
#include <vector>
#include <armadillo>
#include <iostream>
using namespace arma;
class Particle{
    Particle();
    public:
        vector<double> position;
        Particle (vector<double>);

};
#endif // PARTICLE_H
