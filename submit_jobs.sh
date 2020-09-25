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


matlab -nodisplay -nodesktop -r "pp_src_pupil_power_correlations; exit"  1> ~/jobs/$PBS_JOBID.out 2> ~/jobs/$PBS_JOBID.err

#export _JAVA_OPTIONS=-Djava.io.tmpdir=/mnt/homes/home024/btalluri/tmp_thms/
#java -XshowSettings  1> ~/jobs/$PBS_JOBID.out 2> ~/jobs/$PBS_JOBID.err
