#!/bin/bash
#SBATCH -J mergeL
#SBATCH -o logs_mergeL.%j.log      # Name of stdout output file
#SBATCH -e logs_mergeL.%j.log      # Name of stdout output file
#SBATCH -p skx-normal
#SBATCH -N 1                     # Total # of nodes 
#SBATCH -n 20                    # Total # of mpi tasks
#SBATCH -t 00:40:00              # Run time (hh:mm:ss)
#SBATCH --mail-user=bustardchad@gmail.com
#SBATCH --mail-type=all         # Send email at begin and end of job
set -e

module purge
module load intel/18.0.2  impi/18.0.2

module load launcher/3.9

# Merges vtk files and plots data

# Plot only every nth, 0 = off
plotn=0

plotdir=analysis/plots/01/

launcher_file="merge.launcher"

# Join vtk executable
join_vtk_exc=/work2/05230/bustard/stampede2/athena/athena-maxbg-fork/vis/vtk/join_vtk.x
# Join particles executable
join_part_exc=/work2/05230/bustard/stampede2/athena/athena-maxbg-fork/vis/particle/join_lis
# Plotting script
plot_script=~/mhd_cloud/analysis/plot2d.py

################################################################################
echo "Creating launcher file for merging of output files for $jobid"
date

# put all id* folders in new out folder
mkdir -p out
mv id* out

mkdir -p $plotdir
mkdir -p logs
mkdir -p merged
rm -rf $launcher_file
touch $launcher_file

nproc=`ls -d out/id*|wc -l`

nfiles=`ls out/id0/*.vtk|wc -l`
nfiles_merged=`ls merged/*.vtk|wc -l`
if [ "$nfiles" -ne "$nfiles_merged" ]; then
    for i in `seq -f '%04g' 0 $((nfiles-1))`; do
        outfile=merged/cloud.${i}.vtk
	      infile=out/id0/*.${i}.vtk
	      if ! [ -a $infile ]; then
	          echo "$infile does not exist!"
	          continue
	      fi
        if ! [ -a $outfile ]; then
            printf "echo \"Merging $i --> $outfile\"" >> $launcher_file
            # Merge vtk files (start several at the same time)
            printf " && $join_vtk_exc -o $outfile out/id*/*.${i}.vtk &> logs/merge_${i}.log" >> $launcher_file
	          if [ "$plotn" -gt 0 ] && [ "$((10#$i % $plotn))" -eq 0 ]; then # force base 10 
		            printf " && python $plot_script $plotdir $outfile &> logs/plot_${i}.log" >> $launcher_file
	          fi
	          printf " && date && echo \"Done with $i\"\n" >> $launcher_file

	          # Merging particle files, too -- if they exist
	          if ls out/id0/*${i}.part01.lis &> /dev/null; then
		            echo "echo \"Merging particles $i\" && cd out && $join_part_exc -p $nproc -i cloud -f $i:$i:1 -o cloud -d ../merged/ -s part01 &> ../logs/merge_part_${i}.log && cd .. && date && echo "Done with merging particles $i"" >> $launcher_file
	          fi
	      fi
    done
fi

## Merge also other levels if they exist
for ilev in `seq 10`; do
    if ls out/*/lev${ilev}/*vtk &> /dev/null; then
        a=`ls out/*/lev${ilev}/*.vtk|head -n1`
        b=`dirname $a`
        cnfiles=`ls ${b}/*.vtk|wc -l`
        mkdir -p merged/lev${ilev}
        for i in `seq -f '%04g' 0 $((cnfiles-1))`; do
            outfile=merged/lev${ilev}/cloud-lev${ilev}.${i}.vtk
            if ! [ -a $outfile ]; then
                printf "echo \"Merging lev${ilev}, $i --> $outfile\"" >> $launcher_file
                printf " && $join_vtk_exc -o $outfile out/id*/lev${ilev}/*.${i}.vtk &> logs/merge_lev${ilev}_${i}.log" >> $launcher_file
	              printf " && date && echo \"Done with lev${ilev}, $i\"\n" >> $launcher_file
            fi
        done
    fi
done

# Launcher bug workaround
echo "" >> $launcher_file 
  
date
echo "Created job script $launcher_file containing `wc -l $launcher_file` lines."
echo "Starting launcher"


export LAUNCHER_PPN=20
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=$launcher_file

$LAUNCHER_DIR/paramrun

# do the plotting
# change to the plots directory
cd plots

# compute clumps
mkdir -p clumpDir
mpirun -np 4 python /work2/05230/bustard/stampede2/athena/athena_master/analysis_forLeadingArm/compute_clumps_Chad_Tcl1e4.py

mpirun -np 20 python dens_slice.py
mpirun -np 20 python pres_slice.py
mpirun -np 20 python temp_slice.py
mpirun -np 20 python vel_slice.py


