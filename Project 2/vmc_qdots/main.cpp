#include <iostream>
#include <fstream>
#include <string>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include "system.h"
#include "split.cpp"
using namespace std;

int main(int argc, char *argv[]){
    //argv[0] = infile, argv[1] = outfile
    string infilename = argv[1];
    string outfilename = argv[2];
    ifstream infile(infilename);
    string line;
    getline(infile, line);
    //cout << line << endl;
    //cout << argc << " " << infilename << " " << outfilename << " " << endl;
    vector<string> params = split(line, ',');
    //attempts to print vector
    //copy(x.begin(), x.end(), std::ostream_iterator<char>(std::cout, " "));
    //for (auto const& c : x)
    //    cout << c << ' ';
    long number_of_particles = stoi(params[0]);
    long dimensions = stoi(params[1]);
    double w = stof(params[2]);
    cout << "***System Parameters***" << endl
         << "Number of Particles: " << number_of_particles << endl
         << "Dimensions: " << dimensions << endl
         << "Harmonic Oscillator Frequency w: " << w << endl;
    long mc_cycles = stoi(params[3]);
    double step_length = stof(params[4]);
    double initial_alpha = stof(params[5]);
    double alpha_step = stof(params[6]);
    long alpha_variations = stoi(params[7]);
    double initial_beta = stof(params[8]);
    double beta_step = stof(params[9]);
    long beta_variations = stoi(params[10]);
    bool jastrow_bool;
    if (params[11] == "true"){
        jastrow_bool = true;
    }
    else {
        jastrow_bool = false;
    }
    System system(number_of_particles, dimensions, w);
    infile.close();

    return 0;
}
