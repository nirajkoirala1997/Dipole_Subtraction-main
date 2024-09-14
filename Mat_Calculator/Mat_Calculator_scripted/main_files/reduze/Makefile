# an example for computing reductions with Reduze 2

#REDUZECMD=reduze
# activate this for parallel version
REDUZECMD=mpirun -np 4 /home/pulak/Work/Softwares/Install/reduze/bin/reduze

all:
	${REDUZECMD} jobs_reduce_sectors.yaml

clean:
	rm -rf log/* graphs/* reductions/* sectormappings/* tmp/*
