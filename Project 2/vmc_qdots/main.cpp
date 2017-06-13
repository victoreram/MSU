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
    vector<string> x = split(line, ',');
    //attempts to print vector
    //copy(x.begin(), x.end(), std::ostream_iterator<char>(std::cout, " "));
    //for (auto const& c : x)
    //    cout << c << ' ';
    long number_of_particles = stoi(x[0]);
    long dimensions = stoi(x[1]);
    double w = stoi(x[2]);
    cout << "***System Parameters***" << endl
         << "Number of Particles: " << number_of_particles << endl
         << "Dimensions: " << dimensions << endl
         << "Harmonic Oscillator Frequency w: " << w << endl;
    System system(number_of_particles, dimensions, w);
    infile.close();

    return 0;
}
