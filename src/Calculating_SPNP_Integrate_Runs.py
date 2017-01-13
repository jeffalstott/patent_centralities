
# coding: utf-8

# In[1]:

data_directory = '../data/centralities/'
class_system = 'USPC'
controls_directory = data_directory+'controls/%s/'%class_system


# In[2]:

n_controls = 230


# In[3]:

import pandas as pd
from pylab import *


# In[4]:

import gc
from time import time


# In[5]:

empirical = pd.read_hdf(data_directory+'empirical.h5', 'df')


# In[ ]:

def running_stats(empirical,
                  df_name,
                  file_name,
                  controls_directory=controls_directory,
                  n_controls=n_controls,
                 ):
    M = None
    all_max = None
    all_min = None
    t = time()
    for randomization_id in range(n_controls):

        if not randomization_id%10:
            print(randomization_id)
            print("%.0f seconds"%(time()-t))
            t = time()
        
        f = '%s_%i.h5'%(file_name, randomization_id)
        try:
            x = pd.read_hdf(controls_directory+f, df_name)
        except:
            print("Data not loading for %s. Continuing."%f)
            continue
            

        if M is None:
            M = x
            S = 0
#             all_max = M
#             all_min = M
            p = 0
            k = 0
            continue
        k += 1.0
        M_previous = M
#         M = M_previous.add( x.subtract(M_previous)/k )
#         S = ( x.subtract(M_previous).multiply( x.subtract(M) ) ).add(S)
        M = M_previous+((x-M_previous)/k)
        S += (x-M_previous)*(x-M)
        p += (empirical>x).astype('int')
#         all_max = maximum(all_max, x)
#         all_min = minimum(all_min, x)
        gc.collect()  
#     standard_deviation = sqrt(S/(k-1))

    return M, S, p, k


# In[ ]:

M, S, p, k = running_stats(empirical, 'df', 'synthetic_control', controls_directory)


# In[ ]:

standard_deviation = sqrt(S/(k-1))
empirical_percentile = p/k


# In[ ]:

for df in [M, standard_deviation, empirical_percentile]:
    df['patent_number'] = empirical['patent_number']
    df['filing_year'] = empirical['filing_year']


# In[ ]:

store = pd.HDFStore(data_directory+'summary_statistics.h5',
                        mode='a', table=True)
store.put('/randomized_mean_%s'%(class_system), M, 'table', append=False)
store.put('/randomized_std_%s'%(class_system), standard_deviation, 'table', append=False)

store.put('/empirical_percentile_%s'%(class_system), empirical_percentile, 'table', append=False)
# store.put('/randomized_min_%s'%(class_system), all_min, 'table', append=False)

z_scores = (empirical-M)/standard_deviation

z_scores.values[where(z_scores==inf)]=nan 
z_scores.values[where(z_scores==-inf)]=nan 
#All the cases where the z-scores are inf is where the 1,000 randomized controls said there should be 0 deviation, BUT
#the empirical case was different anyway. In each of these cases, the empirical case was JUST slightly off. Sometimes
#a floating point error, and sometimes off by 1 (the minimal amount for citation counts). We shall treat this as not actually
#deviating, and so it becomes 0/0, which is equal to nan.

store.put('/empirical_z_scores_%s'%(class_system), z_scores, 'table', append=False)

store.put('/empirical', empirical, 'table', append=False)

store.close()


# In[ ]:

cp /home/jeffrey_alstott/technoinnovation/patent_centralities/data/centralities/summary_statistics.h5 /home/jeffrey_alstott/Dropbox/TechnoInnovation/

