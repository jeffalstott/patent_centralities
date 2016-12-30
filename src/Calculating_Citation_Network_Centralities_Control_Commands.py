
# coding: utf-8

# In[ ]:

from os import path
import sys
python_location = path.dirname(sys.executable)+'/python'
basic_program = open('Calculating_Citation_Network_Centralities.py', 'r').read()


# In[ ]:

data_directory = '../data/'
# job_type = 'citations'
# job_type = 'cooccurrence'
# entity = 'Firm'
job_label = job_type

try:
    job_label += '_'+entity
except NameError:
    pass


# In[ ]:

n_randomizations = 1000

first_rand_id = None

if first_rand_id is not None:
    runs = range(first_rand_id,first_rand_id+n_randomizations)
else:
    from os import listdir
    dirlist = listdir(data_directory+'centralities/%s/controls/%s/'%(job_type,class_system))
    from pylab import *
    unrun_iterations = ones(n_randomizations)

    for f in dirlist:
        if f.startswith('synthetic_control'):
            n = int(f.split('_')[-1][:-3])
            unrun_iterations[n] = 0

    unrun_iterations = where(unrun_iterations)[0]

    runs = unrun_iterations


# In[ ]:

from os import system

for randomization_id in runs:
    header = """#!{3}
#PBS -l nodes=1:ppn=1
#PBS -l walltime=3:00:00
#PBS -l mem=20000m
#PBS -N rand_{0}_{1}_{2}

randomization_id = {0}
print("Randomization number: %i"%randomization_id)
""".format(randomization_id, class_system, job_label, python_location)

    options = """
class_system = '{0}'
target_years = 'all'
n_years = 'all'

data_directory = '{1}'
randomized_control = True
n_years = None
""".format(class_system, data_directory)

    this_program = header+options+basic_program

    this_job_file = 'jobfiles/randomization_{0}_{1}_{2}.py'.format(randomization_id, class_system, job_label)


    f = open(this_job_file, 'w')
    f.write(this_program)
    f.close()

    system('qsub '+this_job_file)

