#!/bin/sh

# embedded options to qsub - start with #PBS
# walltime: defines maximum lifetime of a job
# nodes/ppn: how many nodes? how many cores?

#PBS -q batch
#PBS -l walltime=700:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb


# -- run in the current working (submission) directory --
cd $PBS_O_WORKDIR

chmod g=wx $PBS_JOBNAME

# FILE TO EXECUTE

sleep "$var"

source activate pp_fooof
python pp_src_hh_fooof.py 1> ~/jobs/$PBS_JOBID.out 2> ~/jobs/$PBS_JOBID.err

#matlab -nodisplay -nodesktop -r "pp_mue_src_fooof; exit"  1> ~/jobs/$PBS_JOBID.out 2> ~/jobs/$PBS_JOBID.err

#java -XshowSettings  1> ~/jobs/$PBS_JOBID.out 2> ~/jobs/$PBS_JOBID.err
