#-
off statistics,finalstats,allwarnings;
nwrite statistics;

#include ../../main_files/def.h
#include input.h
#include ../../main_files/feyn.h
*#include mandelsterm.h
#include ../../main_files/grfunc.h
*#include ../../main_files/SOn.prc
*#include ../../main_files/SUn.prc
*#include ../../main_files/color.h
*#include ../../main_files/gamma5.h
#include ../../main_files/Camplitude.h
#include ../../main_files/amplitude.h
.sort

l mat  = amp*ampc;
.sort
#call grfunc
.sort
print,mat;
.end


