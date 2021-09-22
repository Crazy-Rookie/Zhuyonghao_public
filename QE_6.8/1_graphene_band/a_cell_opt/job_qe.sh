#! /bin/bash
#PBS -l nodes=node8:ppn=8
#PBS -N zyh_test
#PBS -o stand.log
#PBS -e stand.err
#PBS -q long

cd $PBS_O_WORKDIR

# intel.2017
source /opt/longrun/intel2017/compilers_and_libraries_2017.5.239/linux/bin/compilervars.sh intel64
source /opt/longrun/intel2017/mkl/bin/mklvars.sh intel64
source /opt/longrun/intel2017/impi/2017.4.239/intel64/bin/mpivars.sh

# qe_6.8
export PATH=/data/long08/soft/QE_6.8/bin:$PATH

mpirun -np 8 pw.x < cell_opt.in > out.log
