#!/bin/bash

#SBATCH --job ISCU_M140I
#SBATCH --nodes=1
#SBATCH --partition=dept_gpu
#SBATCH --gres=gpu:1


echo Running on `hostname`
echo workdir $PBS_O_WORKDIR
echo ld_library_path $LD_LIBRARY_PATH

cd $SLURM_SUBMIT_DIR

#scratch drive folder to work in
SCRDIR=/scr/${SLURM_JOB_ID}

module load amber/20

#if the scratch drive doesn't exist (it shouldn't) make it.
if [[ ! -e $SCRDIR ]]; then
        mkdir $SCRDIR
fi

chmod +rX $SCRDIR

echo scratch drive ${SCRDIR}

cp $SLURM_SUBMIT_DIR/*.in ${SCRDIR}
cp $SLURM_SUBMIT_DIR/*.prmtop ${SCRDIR}
cp $SLURM_SUBMIT_DIR/*_md2.rst ${SCRDIR}
cp $SLURM_SUBMIT_DIR/*.inpcrd ${SCRDIR}
cp $SLURM_SUBMIT_DIR/dist.RST ${SCRDIR}

cd ${SCRDIR}


#setup to copy files back to working dir on exit
trap "mv *md3.nc $SLURM_SUBMIT_DIR" EXIT

#run the MD, default to name of directory
prefix=complexwfe_M140I
pmemd.cuda -O -i ${prefix}_md3.in -o $SLURM_SUBMIT_DIR/${prefix}_${SLURM_ARRAY_TASK_ID}_md3.out -p ${prefix}.prmtop -c ${prefix}_md2.rst -r ${prefix}_${SLURM_ARRAY_TASK_ID}_md3.rst -x ${prefix}_${SLURM_ARRAY_TASK_ID}_md3.nc -inf $SLURM_SUBMIT_DIR/mdinfo.${SLURM_ARRAY_TASK_ID}
