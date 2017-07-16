#!/usr/bin/python
import math
import sys
import os
import re
import time
#======================================================================
# functions
#======================================================================
def write_submission_script(filename, rscript, id, start_point):
    submit_file = open(filename,'w')
    submit_str = """#!/bin/csh
#$ -q long # Specify queue
#$ -N ar.set   # Specify job name
#$ -t 1 # Specify number of tasks in array

fsync $SGE_STDOUT_PATH &

setenv RPATH /afs/crc.nd.edu/user/q/qtenbosc/Movement/Scripts

cd $RPATH

module load bio/R/3.2.3-gcc

Rscript %s %d %d

""" % (rscript,id,start_point)

    submit_file.write(submit_str)
    submit_file.close()

#======================================================================
# Generate submission files
#======================================================================
for i in range(1,11):
    for j in range(1,6):
        submit_file_name = 'submit_relsd_25_%d_%d.sh' % (i,j)
        write_submission_script(submit_file_name,'./ArraySet_RelSD_25.R',i,j)
        print('submitting jobs i:%d j:%d' % (i,j))
        os.system("qsub %s" % (submit_file_name))
        time.sleep(1)
