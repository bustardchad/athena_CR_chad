#!/bin/bash
#This file is called submit-script.sh

#SBATCH -J CRTurb
#SBATCH -p development	# default "univ" if not specified
#SBATCH -t 0-00:03:00		# run time in days-hh:mm:ss
#SBATCH -N 1			# number of nodes
#SBATCH -n 1				# TOTAL number of cores
#SBATCH -A TG-PHY210004
#SBATCH -e output/test.err
#SBATCH -o output/test.out
#SBATCH --mail-user=bustardchad@gmail.com
#SBATCH --mail-type=all


#Now list your executable command (or a string of them).
# Example for non-SLURM-compiled code:
#module load fftw-3.3.4  
#module load mpi/intel/mpich-3.1
#module load hdf5-1.8.11 
#module load compile/intel
#mpirun -np 1 /home/bustard/flash/Flash4.2_CR/object_LMC_allTogether_KE/flash4 -par_file wind_flash_refine.p

module load impi/18.0.2 
#module load phdf5/1.10.4
pwd
date


#ibrun /work/05230/bustard/stampede2/athena/athena_master/bin/athena -i athinput.cr_cloud
#ibrun /work/05230/bustard/stampede2/athena/athena_master/bin_diffuseFlux_ionneutral_March28/athena -i params
ibrun /work/05230/bustard/stampede2/athena/athena_master/bin_lcrTest/athena -i params -t 00:58:00
