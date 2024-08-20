#include <iostream>
#include "ampvv.h"

using namespace std;

void ZZ_computeAmplitudes(double sn, double un, int order, int debug, double &M0M0, double &M0M1, double &M1M1, double &M0M2a, double &M0M2b, double &M0M2c, double &M0M2d) {
 // double sn = 15.488519990310740605;
 // double un = -8.23989;
//  double sn = 6.863701277323974;
//  double un = -1.83989;i

//  sn=57.19641336981816522753;
//  un=-9.58522923969884403871;

//  double m0m0, m0m1, m1m1, m0m2a, m0m2b, m0m2c, m0m2d;

  ampvv(sn, un, M0M0, M0M1, M1M1, M0M2a, M0M2b, M0M2c, M0M2d);

  //  cout << "libZZ (2014-04-29.1): " << M0M0 << "\t" << M0M1 << "\t" << M1M1 << "\t" << M0M2a << "\t" << M0M2b << "\t" << M0M2c << "\t" << M0M2d << endl;
}

