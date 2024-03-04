cd virtual_standalone
make
time ./runVir
cd ../dipole_standalone/
make
time ./run.x
cd ../PK_standalone/
make
time ./runPK

