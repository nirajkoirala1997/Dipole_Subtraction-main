#include "ampvveval.h"

// auxiliary functions to split up the code into compiler friendly small pieces
void ampvv01(numtype a[]);
void ampvv02(numtype a[]);
void ampvv03(numtype a[]);
void ampvv04(numtype a[]);
void ampvv05(numtype a[]);
void ampvv06(numtype a[]);
void ampvv07(numtype a[]);
void ampvv08(numtype a[]);
void ampvv09(numtype a[]);
void ampvv10(numtype a[]);
void ampvv11(numtype a[]);
void ampvv12(numtype a[]);
void ampvv13(numtype a[]);
void ampvv14(numtype a[]);
void ampvv15(numtype a[]);
void ampvv16(numtype a[]);
void ampvv17(numtype a[]);
void ampvv18(numtype a[]);
void ampvv19(numtype a[]);
void ampvv20(numtype a[]);
void ampvv21(numtype a[]);
void ampvv22(numtype a[]);
void ampvv23(numtype a[]);
void ampvv24(numtype a[]);
void ampvv25(numtype a[]);
void ampvv26(numtype a[]);
void ampvv27(numtype a[]);
void ampvv28(numtype a[]);
void ampvv29(numtype a[]);
void ampvv30(numtype a[]);
void ampvv31(numtype a[]);
void ampvv32(numtype a[]);
void ampvv33(numtype a[]);
void ampvv34(numtype a[]);
void ampvv35(numtype a[]);
void ampvv36(numtype a[]);
void ampvv37(numtype a[]);
void ampvv38(numtype a[]);
void ampvv39(numtype a[]);
void ampvv40(numtype a[]);
void ampvv41(numtype a[]);
void ampvv42(numtype a[]);
void ampvv43(numtype a[]);
void ampvv44(numtype a[]);
void ampvv45(numtype a[]);
void ampvv46(numtype a[]);
void ampvv47(numtype a[]);
void ampvv48(numtype a[]);
void ampvv49(numtype a[]);
void ampvv50(numtype a[]);
void ampvv51(numtype a[]);
void ampvv52(numtype a[]);
void ampvv53(numtype a[]);
void ampvv54(numtype a[]);
void ampvv55(numtype a[]);
void ampvv56(numtype a[]);
void ampvv57(numtype a[]);
void ampvv58(numtype a[]);
void ampvv59(numtype a[]);
void ampvv60(numtype a[]);
void ampvv61(numtype a[]);
void ampvv62(numtype a[]);
void ampvv63(numtype a[]);
void ampvv64(numtype a[]);
void ampvv65(numtype a[]);
void ampvv66(numtype a[]);
void ampvv67(numtype a[]);
void ampvv68(numtype a[]);
void ampvv69(numtype a[]);
void ampvv70(numtype a[]);
void ampvv71(numtype a[]);
void ampvv72(numtype a[]);
void ampvv73(numtype a[]);
void ampvv74(numtype a[]);
void ampvv75(numtype a[]);
void ampvv76(numtype a[]);
void ampvv77(numtype a[]);
void ampvv78(numtype a[]);
void ampvv79(numtype a[]);
void ampvv80(numtype a[]);
void ampvv81(numtype a[]);
void ampvv82(numtype a[]);
void ampvv83(numtype a[]);
void ampvv84(numtype a[]);
void ampvv85(numtype a[]);
void ampvv86(numtype a[]);
void ampvv87(numtype a[]);
void ampvv88(numtype a[]);
void ampvv89(numtype a[]);
void ampvv90(numtype a[]);
void ampvv91(numtype a[]);
void ampvv92(numtype a[]);
void ampvv93(numtype a[]);

// the main function
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

   

