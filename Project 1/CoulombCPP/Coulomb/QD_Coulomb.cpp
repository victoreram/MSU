#include "Coulomb_Functions.hpp"
#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{
  for (unsigned int i = 0; i < 10; i++){
  cout << argv[i] << endl;
  }
  if(argc != 10){ std::cerr << "Wrong Input: should be ./QD_Coulomb  hw  n1  ml1  n2  ml2  n3  ml3  n4  ml4" << std::endl; exit(1); }
  double hw = std::atof(argv[1]);
  int n1 = std::atoi(argv[2]);
  int ml1 = std::atoi(argv[3]);
  int n2 = std::atoi(argv[4]);
  int ml2 = std::atoi(argv[5]);
  int n3 = std::atoi(argv[6]);
  int ml3 = std::atoi(argv[7]);
  int n4 = std::atoi(argv[8]);
  int ml4 = std::atoi(argv[9]);

  double TBME = Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4);
  std::cout << std::setprecision(12);
  std::cout << "< " << n1 << "," << ml1 << " ; " << n2 << "," << ml2 << " || V || " << n3 << "," << ml3 << " ; " << n4 << "," << ml4 << " > = " << TBME << std::endl;
  
  return 0;
}
