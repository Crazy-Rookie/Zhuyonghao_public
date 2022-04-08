#!/bin/bash
#SBATCH -p vip_07
#SBATCH -N 1
#SBATCH -n 24

source /public3/soft/other/vasp/cn-module-vasp.5.4.4.sh
module load vasp/intel-17/vasp544-O1-vtst
module unload mpi/intel/17.0.7
module load mpi/intel/18.0.2

# pawpyseed
source activate pawpyseed

VASPCOMMAND='srun vasp_std'

#############################################################
## NOTE: 3_all_elec*0.33=eV; 2_phase_correc*0.33=eV; 1_normal(eV)

## Set Up Parameters
MINB=340     # lowest energy orbital from which electrons can excite
MAXB=410     # highest energy orbital to which electrons can go 
MINT=2000    # starting step
StartTime=2000 # generally =MINIT
MAXT=2002      # ending step
TMSTP=1        # time step
#############################################################

## mkdir cio_npy in ./run_nac/3_all_elec
if [ ! -d ./run_nac/3_all_elec/cio_npy ]; then mkdir ./run_nac/3_all_elec/cio_npy; fi
## mkdir ipr and overlap in ./Ipr_Overlap
if [ ! -d ./Ipr_Overlap/ipr ]; then mkdir ./Ipr_Overlap/ipr; fi
if [ ! -d ./Ipr_Overlap/overlap ]; then mkdir ./Ipr_Overlap/overlap; fi
## mkdir procar_npy in ./procar
if [ ! -d ./procar/procar_npy ]; then mkdir ./procar/procar_npy; fi
#############################################################

printf "\n======= INITIAL PARAMETERS (MINBAND MAXBAND MINTIME MAXTIME TIMESTEP) ========\n" 
printf "           %10d%10d%10d%10d%10d\n\n" $MINB $MAXB $MINT $MAXT $TMSTP

## FIRST RUN
# StartTime = MINT ?
if [ "$StartTime" -eq "$MINT" ];then 
	SUFX=$( printf "%04d" "$MINT" ) &&
	cd ./scf_run/$SUFX &&
	printf "Running VASP at t = $MINT fs and getting initial set of good orbitals\n" &&
	$VASPCOMMAND &&
	
	# workfunction
	vaspkit -task 426 &&
    sleep 2
	mv PLANAR* ../../workfunction/"wf-${SUFX}" &&
	
	# chgcar
	cp ../../chgcar/WF_liner_average.py ./ &&
	python WF_liner_average.py &&
	mv  Liner* ../../chgcar/"chg-${SUFX}" &&
	
	# procar
	cd ../../procar/ &&
	if [ -f ./PROCAR ]; then rm PROCAR; fi
	ln -s ../scf_run/$SUFX/PROCAR &&
	printf "read procar: Getting procar.npy at t = $time fs \n" &&
	if [ -f ./procar.npy ]; then rm procar.npy; fi
	python ReadProcar.py &&
	mv band_energy.npy ./procar_npy/band_energy.$SUFX.npy &&
	mv procar.npy ./procar_npy/procar.$SUFX.npy &&
	rm PROCAR &&
	cd ../scf_run/$SUFX &&
	
	# remove files
	rm -f WF_liner* OUTCAR LOCPOT CHG CHGCAR DOSCAR  &&
	rm -f EIGENVAL IBZKPT PCDAT REPORT XDATCAR &&
	printf "\n=======================================\n" 
	cd ../../
	(( StartTime++ ))
fi
#############################################################

## START CALCULATING THE COUPLiNGS ALONG THE TRAJECTORY
for time in `seq $StartTime $MAXT`
do
# run vasp
	sleep 2
	# define prefix
	SUFX_1=$( printf "%04d" "$time" ) &&
	SUFX_tmp=$((time-1))
	SUFX_0=$( printf "%04d" "$SUFX_tmp" ) &&
	# check and copy WAVECAR
	# ./scf_run/$SUFX_0/WAVECAR
	if [ ! -f ./scf_run/$SUFX_0/WAVECAR ]; then 
		printf "\n--Error! There was no WAVECAR in $SUFX_0! Exiting....--\n"
		exit
	fi
	# ./scf_run/$SUFX_1/WAVECAR
	if [ ! -f ./scf_run/$SUFX_1/WAVECAR ]; then 
		cp ./scf_run/$SUFX_0/WAVECAR ./scf_run/$SUFX_1
	fi
	# run vasp
	cd ./scf_run/$SUFX_1 &&
	printf "Running vasp at t = $time fs \n" &&
	$VASPCOMMAND &&
	sleep 1
	
	# workfunction
	vaspkit -task 426 &&
    sleep 2
	mv PLANAR* ../../workfunction/"wf-${SUFX_1}" &&
	
	# chgcar
	cp ../../chgcar/WF_liner_average.py ./ &&
	python WF_liner_average.py &&
	mv  Liner* ../../chgcar/"chg-${SUFX_1}" &&
	
	# remove files
	rm -f WF_liner* OUTCAR LOCPOT CHG CHGCAR DOSCAR  &&
	rm -f EIGENVAL IBZKPT PCDAT REPORT XDATCAR &&
	printf "\n=======================================\n" 
	cd ../../
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# nac: 1_normal
	cd ./run_nac/1_normal &&
	if [ ! -f WAVECAROLD ]; then ln -s ../../scf_run/$SUFX_0/WAVECAR  WAVECAROLD; fi &&
	if [ ! -f WAVECARNEW ]; then ln -s ../../scf_run/$SUFX_1/WAVECAR  WAVECARNEW; fi &&
	printf "1_normal: Getting NAC at t = $time fs \n" &&
	if [ ! -f ./real$SUFX_1 ]; then ./ovlap_NORM_OS $MINB $MAXB $TMSTP $SUFX_1; fi
	if [ ! -f ./real$SUFX_1 ]; then exit; fi
	if [ -f energy_by_band ]; then  cat energy_by_band >> energy; fi &&
	if [ -f WAVECAROLD ]; then rm WAVECAROLD; fi &&
	if [ -f WAVECARNEW ]; then rm WAVECARNEW; fi &&	
	cd ../../ &&
	printf "\n=======================================\n" 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# nac: 2_phase_correc
	cd ./run_nac/2_phase_correc/ &&
	if [ ! -d ./$SUFX_0 ]; then mkdir ./$SUFX_0; fi
	cd ./$SUFX_0 &&
	if [ ! -f ./WAVECAR ]; then ln -s ../../../scf_run/$SUFX_0/WAVECAR; fi
	cd ../ &&

	if [ ! -d ./$SUFX_1 ]; then mkdir ./$SUFX_1; fi
	cd ./$SUFX_1 &&
	if [ ! -f ./WAVECAR ]; then ln -s ../../../scf_run/$SUFX_1/WAVECAR; fi
	cd ../ &&

	printf "2_phase_correc: Getting NAC at t = $time fs \n" &&
	if [ ! -f ./real$SUFX_1 ]; then python nac.py $MINB $MAXB $time; fi
	if [ ! -f ./real$SUFX_1 ]; then exit; fi
	if [ -d ./$SUFX_0 ]; then rm -r $SUFX_0; fi
	cd ../../ &&
	printf "\n=======================================\n" 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# nac: 3_all_elec
	cd ./run_nac/3_all_elec/ &&
	if [ ! -d ./$SUFX_0 ]; then mkdir ./$SUFX_0; fi
	cd ./$SUFX_0 &&
	if [ ! -f ./WAVECAR ]; then     ln -s ../../../scf_run/$SUFX_0/WAVECAR; fi
	if [ ! -f ./POTCAR ]; then      ln -s ../../../scf_run/$SUFX_0/POTCAR; fi
	if [ ! -f ./vasprun.xml ]; then ln -s ../../../scf_run/$SUFX_0/vasprun.xml; fi
	if [ ! -f ./CONTCAR ]; then     ln -s ../../../scf_run/$SUFX_0/CONTCAR; fi
	cd ../ &&

	if [ ! -d ./$SUFX_1 ]; then mkdir ./$SUFX_1; fi
	cd ./$SUFX_1 &&
	if [ ! -f ./WAVECAR ]; then     ln -s ../../../scf_run/$SUFX_1/WAVECAR; fi
	if [ ! -f ./POTCAR ]; then      ln -s ../../../scf_run/$SUFX_1/POTCAR; fi
	if [ ! -f ./vasprun.xml ]; then ln -s ../../../scf_run/$SUFX_1/vasprun.xml; fi
	if [ ! -f ./CONTCAR ]; then     ln -s ../../../scf_run/$SUFX_1/CONTCAR; fi
	cd ../ &&
	
	printf "3_all_elec: Getting NAC at t = $time fs \n" &&
	sleep 1
	if [ ! -f ./real$SUFX_1 ]; then python paw2.py $MINB $MAXB $time $((MINT++)); fi
	sleep 1
	if [ ! -f ./real$SUFX_1 ]; then exit; fi
	if [ -d ./$SUFX_0 ]; then rm -r $SUFX_0; fi
	cd ../../ &&
	printf "\n=======================================\n" 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# procar: save porcar.npy
	if [ -f ./scf_run/$SUFX_1/PROCAR ]; then
		cd ./procar/ &&
		if [ ! -f PROCAR ]; then ln -s ../scf_run/$SUFX_1/PROCAR; fi
		printf "read procar: Getting procar.npy at t = $time fs \n" &&
		if [ -f ./procar.npy ]; then rm procar.npy; fi
		python ReadProcar.py &&
		if [ ! -f band_energy.npy ]; then exit; fi
		if [ ! -f procar.npy ]; then exit; fi
		mv band_energy.npy ./procar_npy/band_energy.$SUFX_1.npy &&
		mv procar.npy ./procar_npy/procar.$SUFX_1.npy &&
		if [ -f PROCAR ]; then rm PROCAR; fi &&
		cd ../
	fi 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ipr and overlap: ipr_overlap
	if [ -f ./scf_run/$SUFX_1/WAVECAR ]; then
		cd ./Ipr_Overlap/ &&
		if [ ! -f WAVECAR ]; then ln -s ../scf_run/$SUFX_1/WAVECAR; fi
		printf "calc ipr: Getting IPR at t = $time fs \n" &&
		python Ipr_Overlap.py $MINB $MAXB &&
		if [ ! -f ipr.dat ]; then exit; fi     &&
		if [ ! -f overlap.dat ]; then exit; fi &&
		mv ipr.dat ./ipr/ipr.$SUFX_1.dat &&
		mv overlap.dat ./overlap/overlap.$SUFX_1.dat &&
		if [ -f WAVECAR ]; then rm WAVECAR; fi &&
		cd .. &&
		printf "\n=======================================\n" 
	fi
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if [ -f ./scf_run/$SUFX_0/WAVECAR ];     then rm ./scf_run/$SUFX_0/WAVECAR; fi
	if [ -f ./scf_run/$SUFX_0/CONTCAR ];     then rm ./scf_run/$SUFX_0/CONTCAR; fi
	if [ -f ./scf_run/$SUFX_0/PROCAR ];      then rm ./scf_run/$SUFX_0/PROCAR; fi
	if [ -f ./scf_run/$SUFX_0/vasprun.xml ]; then rm ./scf_run/$SUFX_0/vasprun.xml; fi
	if [ -f ./scf_run/$SUFX_0/OSZICAR ];     then rm ./scf_run/$SUFX_0/OSZICAR; fi
	if [ -f ./scf_run/$SUFX_0/KPOINTS ];     then rm ./scf_run/$SUFX_0/KPOINTS; fi
	
done

(( MINT-- ))
