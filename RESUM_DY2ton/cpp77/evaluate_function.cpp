#include "ampzz/ampvveval.h"
#include <cmath>
#include <list>
#include <fstream>
#include <iomanip>
//#include <iostream>
//#include <ginac/ginac.h>

using namespace std;
//using namespace GiNaC;


void ampvv(double sn, double un,
           double& m0m0fin_N, double& m0m1fin_CF_N, double& m1m1fin_CF2_N,
           double& m0m2fin_CF_N_Nf, double& m0m2fin_CF_N_NfZ2,
           double& m0m2fin_CF2_N, double& m0m2fin_CF)
{

#ifdef DIGITS
  GiNaC::Digits=DIGITS;
#endif

  if (sn < 4)
    throw runtime_error("illegal s value");
  ex b = sqrt(1 - 4/ex(sn));
  numtype x = tonn((1 - b)/(1 + b));
  numtype z = -un;
  if (z < x || z > 1/x)
    throw runtime_error("illegal (s,u) combination");
  numtype y = (1 + x*x - x*z)/x;
  numtype Z = z + 1;
  numtype Y = y + 1;
	
//  cout << " x " << x << " y " << y  << " z " << z << endl;
numtype a[9286];

a[1] = nn(1)+x+x*x+nn(-1)*x*z;
a[2] = nn(1)+x*x+nn(-1)*x*z;
a[3] = nn(1)+nn(-1)*x*z;
a[4] = nn(-1)*x+z;
a[5] = nn(1)+nn(-1)*x*y;
a[6] = nn(-1)*x+y;
a[7] = nn(1)+z;
a[8] = nn(1)+y;
a[9] = nn(1)+nn(-1)*x;
a[10] = nn(1)+x;
a[11] = Y+Z;
a[12] = Y;
a[13] = Z;
a[14] = x;
a[15] = z;
a[16] = y;

#include "./ampzz/amps/aterms.h"
//#include "m11amp.h"

//m0m0fin_N = todouble(a[9279]);
//m0m1fin_CF_N = todouble(a[9280]);
m1m1fin_CF2_N = todouble(a[9285]);
//m0m2fin_CF_N_Nf = todouble(a[9281]);
//m0m2fin_CF2_N = todouble(a[9284]);
//m0m2fin_CF_N_NfZ2 = todouble(a[9283]);
//m0m2fin_CF = todouble(a[9282]);
}
int main() {
	double sn, un;
	ifstream inputFile("/home/vaibhav/work/resummation/RESUM_VP5/inputcpp.input");
	if (!(inputFile >> sn >> un)) {
        cerr << "Error: Unable to read integers from file" << endl;
        inputFile.close();
        return 1; // Exit with error
    }

  double m0m0, m0m1, m1m1, m0m2a, m0m2b, m0m2c, m0m2d;

  ampvv(sn, un, m0m0, m0m1, m1m1, m0m2a, m0m2b, m0m2c, m0m2d);

 // cout <<"m0m0  " <<  m0m0 << "\t" << " m0m1  " << m0m1 << "\t" << "m1m1  " << m1m1 << endl;
  std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << m1m1 << std::endl;

  return 0;
}

