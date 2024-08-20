#include <iostream>
#include "ampvv.h"

using namespace std;

int main() {
  // a generic point we used to check the masters
  // double sn = 15.488519990310740605;
  // double un = -8.23989;
  // a generic point we used to check the amplitude
  double sn = 6.863701277323974;
  double un = -1.83989;
  // a difficult point used by LT and EW (point 2):
  //double sn = 5577.776648726018;     /*x2=0.0001793472183635458*/
  //double un =-0.0006513521259708399; /*z2=0.0006513521259708399*/

  double m0m0, m0m1, m1m1, m0m2a, m0m2b, m0m2c, m0m2d;

  ampvv(sn, un, m0m0, m0m1, m1m1, m0m2a, m0m2b, m0m2c, m0m2d);

  cout << m0m0 << "\t" << m0m1 << "\t" << m1m1 << "\t" << m0m2a << "\t" << m0m2b << "\t" << m0m2c << "\t" << m0m2d << endl;

  return 0;
}

