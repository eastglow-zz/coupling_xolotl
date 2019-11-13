#!/bin/sh
#SBATCH --job-name=FG64_3D            #Job name
#SBATCH --nodes=4                   #Number of nodes (servers, 32 proc/node)
#SBATCH --ntasks=128                 #Number of tasks/How many parallel jobs to run
#SBATCH --ntasks-per-node=32         #Tasks per node
##SBATCH --cpus-per-task=1           #Number of CPU per task
#SBATCH --mem-per-cpu=2gb           #Memory (120 gig/node)
#SBATCH --time=96:00:00             #Walltime hh:mm:ss
#SBATCH --output=FG3D-%j.out #Output and error log name
#SBATCH --mail-type=ALL             #When to email user: NONE,BEGIN,END,FAIL,REQUEUE,ALL
#SBATCH --mail-user=donguk.kim@ufl.edu       #Email address to send mail to
#SBATCH --qos=michael.tonks                #Allocation group name, add -b for burst job
##SBATCH --partition=hpg2-compute
##SBATCH --array=1-200%10           #Used to submit multiple jobs with one submit

APPPATH=~/gcc_moose/projects/coupling_xolotl
srun --mpi=pmix_v1 $APPPATH/coupling_xolotl-opt -i ./GPM_grain_tracker.i
