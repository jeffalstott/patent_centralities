
# coding: utf-8

# In[1]:

from os import path
import sys
python_location = path.dirname(sys.executable)+'/python'
basic_program = open('Calculating_SPNP.py', 'r').read()


# In[2]:

data_directory = '/home/jeffrey_alstott/technoinnovation/patent_centralities/data/'
class_system = 'USPC'


# In[5]:

n_randomizations = 1000

first_rand_id = None

if first_rand_id is not None:
    runs = range(first_rand_id,first_rand_id+n_randomizations)
else:
    from os import listdir
    dirlist = listdir(data_directory+'centralities/controls/%s/'%(class_system))
    from pylab import *
    unrun_iterations = ones(n_randomizations)

    for f in dirlist:
        if f.startswith('synthetic_control'):
            n = int(f.split('_')[-1][:-3])
            unrun_iterations[n] = 0

    unrun_iterations = where(unrun_iterations)[0]

    runs = unrun_iterations


# In[6]:

from os import system

for randomization_id in runs:
    header = """#!{2}
#PBS -l nodes=1:ppn=1
#PBS -l walltime=3:00:00
#PBS -l mem=30000m
#PBS -N rand_{0}_{1}

randomization_id = {0}
print("Randomization number: %i"%randomization_id)
""".format(randomization_id, class_system, python_location)

    options = """
data_directory = '{1}'
class_system = '{0}'
randomized_control = True
save_file_name = 'synthetic_control_%i'%randomization_id
""".format(class_system, data_directory)

    this_program = header+options+basic_program

    this_job_file = 'jobfiles/randomization_{0}_{1}.py'.format(randomization_id, class_system)


    f = open(this_job_file, 'w')
    f.write(this_program)
    f.close()

    system('qsub '+this_job_file)

