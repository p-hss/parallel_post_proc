import fileinput
import sys
import os
import subprocess

#++++++++++++ INPUT +++++++++++++++++
param_file='_param_'
job_file='./jobscripts/tub_phil.job' 

#simulation=['Z01', 'Z02', 'Z03', 'Z04', 'Z05', 'Z06', 'Z07', 'Z08']
simulation=['Z00', 'Z01', 'Z03', 'Z04']
#++++++++++++++++++++++++++++++++++++

def replaceline(file,pattern,replaceit):
    for line in fileinput.input(file, inplace=1):
        if pattern in line:
            line = replaceit 
        sys.stdout.write(line)

subprocess.call(["make", "clean"])
subprocess.call(["make"])
for i in range(len(simulation)):
	
	#change job name
	pattern='#$ -N pproc_'
	replaceit='#$ -N pproc_'+simulation[i]+'\n'
	replaceline(job_file,pattern,replaceit)

	#change job output filename
	pattern='./run < param/'
	replaceit='./run < param/'+simulation[i]+'_param.dat'
	replaceline(job_file,pattern,replaceit)

	subprocess.call(["qsub","jobscripts/tub_phil.job"])

