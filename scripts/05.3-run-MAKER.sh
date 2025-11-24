#!/bin/env bash

#SBATCH --job-name=MAKER
#SBATCH --partition=pibu_el8
#SBATCH --time=7-0
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/MAKER_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/MAKER_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation

mkdir -p $WORKDIR
cd $WORKDIR

# Set scratch directory for temporary files
export SCRATCH="/scratch/${USER}/MAKER_${SLURM_JOB_ID}"
mkdir -p $SCRATCH

REPEATMASKER_DIR="/data/courses/assembly-annotation-course/CDS_annotation/softwares/RepeatMasker"
export PATH=$PATH:"/data/courses/assembly-annotation-course/CDS_annotation/softwares/RepeatMasker"

module load OpenMPI/4.1.1-GCC-10.3.0
module load AUGUSTUS/3.4.0-foss-2021a

# Create TMP directory in scratch
mkdir -p $SCRATCH/MAKER_TMP

# Run MAKER gene annotation pipeline with MPI support.
mpiexec --oversubscribe -n 50 \ # allows processes to run on the same core if needed
  apptainer exec \
    --bind /data \
    --bind $SCRATCH:/TMP \
    --bind $COURSEDIR \
    --bind $AUGUSTUS_CONFIG_PATH \
    --bind $REPEATMASKER_DIR \
    ${COURSEDIR}/containers/MAKER_3.01.03.sif \
    maker -mpi \ # tells MAKER to run in parallel mode
          --ignore_nfs_tmp \ # forces MAKER to use the user-defined temporary path
          -TMP /TMP \ # specify temporary directory inside the container
          maker_opts.ctl \
          maker_bopts.ctl \
          maker_evm.ctl \
          maker_exe.ctl

# Cleanup scratch if job completes successfully
if [ $? -eq 0 ]; then
    echo "Job completed successfully, cleaning scratch directory..."
    rm -rf $SCRATCH
else
    echo "Job failed, preserving scratch directory: $SCRATCH"
fi