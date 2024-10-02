#!/bin/bash
######################################################################################################################################
####################################### No Changes to this part of the batch #########################################################
######################################################################################################################################

#SBATCH --partition rtx3070			### Specify partition name where to run a job.
#SBATCH --time 0-00:01:00			### Job running time limit. Make sure it is not exceeding the partition time limit! Format: D-H:MM:SS
#SBATCH --output job-%J.out			### Output log for running job - %J is the job number variable

##SBATCH --mail-user=user@post.jce.ac.il	### User's email for sending job status
##SBATCH --mail-type=ALL			### Conditions when to send the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --gres=gpu:1       	        	### number of GPUs, ask for more than 1 only if you can parallelize your code for multi GPU
#SBATCH --cpus-per-task=8    			### number of CPU cores
#######################################################################################################################################


######################################################################################################################################
####################################### Start you code below ####  ###################################################################
######################################################################################################################################

#SBATCH --job-name 'hello'		### Name of the job. replace my_job with your desired job name
###SBATCH --nodes 1
###SBATCH --ntasks-per-node=5
#SBATCH --ntasks=5

module load anaconda  			### load anaconda module (must present when working with conda environments)
source activate myenv			### activating environment, environment must be configured before running the job


mpicxx -fopenmp -c main.c -o main.o
/usr/local/cuda/bin/nvcc -gencode arch=compute_86,code=sm_86 -c cudaFunctions.cu -o cudaFunctions.o
mpicxx -fopenmp -o mpiCudaOpemMP  main.o  cudaFunctions.o  -L/usr/local/cuda/lib -L/usr/local/cuda/lib64 -lcudart	
mpirun -np 5 ./mpiCudaOpemMP