
# coding: utf-8

# In[1]:

import pandas as pd
from pylab import *


# In[2]:

data_input = '/home/jeffrey_alstott/Dropbox/Patent_centralities/'
data_output = '/home/jeffrey_alstott/technoinnovation/patent_centralities/data/'


# In[3]:

patents = pd.read_csv(data_input+'PATENT_INFO_1790_2015_no_citcount_no_centrality.csv', 
                 usecols=['patent_number', 'filing_year', 'mainclass_id'])


# In[4]:

def floatify(x):
    from numpy import nan
    try:
        return float(x)
    except ValueError:
        return nan


# In[5]:

# vc = patents['Class'].value_counts()
# from pylab import array
# vc[array(list(map(type, vc.index.values)))==str].sum()


# In[6]:

patents['Class'] = patents['mainclass_id'].apply(floatify)


# In[7]:

patents.dropna(inplace=True)


# In[8]:

del patents['mainclass_id']


# In[9]:

citations = pd.read_csv(data_input+'CITATION_INFO.csv', 
                usecols=['citing_id', 'cited_id', 
                         'filing_year_citing', 'filing_year_cited', 
                        ])


# In[10]:

citations.rename(columns={"citing_id": 'Citing_Patent',
                            "cited_id": 'Cited_Patent',
                            "filing_year_citing": 'Year_Citing_Patent',
                            "filing_year_cited": 'Year_Cited_Patent',
                             }, inplace=True)


# In[11]:

### Drop citations that are erroneous, having a citation that points to a future patent
citations = citations[citations['Year_Citing_Patent']>=citations['Year_Cited_Patent']]


# In[12]:

print(patents.shape)
patents = patents[(patents['patent_number'].isin(citations['Cited_Patent']) | 
                   patents['patent_number'].isin(citations['Citing_Patent'])
                   )]
print(patents.shape)


# In[13]:

patents.to_hdf(data_output+'patents.h5', 'df', complevel=9, complib='blosc')


# In[14]:

patents.set_index('patent_number', inplace=True)


# In[15]:

for patent_type in ['Citing_Patent', 'Cited_Patent']:
    print(patent_type)
    z = empty(citations.shape[0])
    start_ind = 0
    while start_ind<z.shape[0]:
        print(start_ind)
        stop_ind = start_ind+2000000
        z[start_ind:stop_ind] = patents.ix[citations[start_ind:stop_ind][patent_type], 'Class']
        start_ind = stop_ind
    citations['Class_%s'%patent_type] = z


# In[16]:

print(citations.shape)
citations.dropna(inplace=True)
print(citations.shape)


# In[17]:

citations.to_hdf(data_output+'citations.h5', 'df', complevel=9, complib='blosc')

