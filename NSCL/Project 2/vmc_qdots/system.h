#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <armadillo>
#include <iostream>
using namespace arma;
class System{
    System();
    public:
        long number_of_particles;
        int dimensions;
        double w;
        System (long, int, double);

};

#endif // SYSTEM_H
