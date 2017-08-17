
# coding: utf-8

# In[1]:

# data_directory = '../data/'
# randomized_control = False
# randomization_id = 1
class_system = 'USPC'


# # Import Libraries

# In[2]:

import pandas as pd
from pylab import *
import gc
from tqdm import *
import igraph as igraph
import time


# In[3]:

start_time = time.time()


# Define Functions
# ===

# In[4]:

def create_patent_citation_graph(PATENT_INFO, citations):
    ti = time.time()
    print('create igraph object and populate it with patent attributes')
    G = igraph.Graph(directed=True)

    ##### add nodes ####
    print('step 1: add nodes')
    num_nodes = shape(PATENT_INFO)[0]
    G.add_vertices(num_nodes)

    #### add nodes attributes ####
    print('step 2: add node attributes')
    for x in PATENT_INFO.columns:
        G.vs[x] = array(PATENT_INFO[x])
        #if x == 'patent_number': continue
        #else: G.vs[x] = array(PATENT_INFO[x])

    #### add edges ####
    print('step 3: add edges')
    # create a series with 'patent_number' as index and ID from zero to N as values. 
    # We will then use it to translate the edgelist from patent numbers to patent IDs. We will use 'map' for that
    PATENT_INFO['ID'] = range(shape(PATENT_INFO)[0])
    series_patent_ids = PATENT_INFO[['ID','patent_number']].set_index('patent_number').ix[:,0] #s = df.ix[:,0] <-- this create a series out of a dataframe (see: http://stackoverflow.com/questions/15360925/how-to-get-the-first-column-of-a-pandas-dataframe-as-a-series)
    # Now map patent numbers to IDs in our edgelist
    citations['citing_numerical_id'] = citations['Citing_Patent'].map(series_patent_ids) # see: http://stackoverflow.com/questions/25653652/how-to-replace-efficiently-values-on-a-pandas-dataframe
    citations['cited_numerical_id'] = citations['Cited_Patent'].map(series_patent_ids)
    

    # create edgelist with IDs, not patent_numbers (to understand why see:http://igraph.org/python/doc/tutorial/tutorial.html)
    # incidentally, the index of the edgelist dataframe is the edge ID in the igraph object (unless we later remove nodes or edges)
    edgelist = citations[['cited_numerical_id','citing_numerical_id']]

    # now populate our graph G with edges. We will then be able to select subgraphs based on nodes (or edges) attributes
    G.add_edges(edgelist.to_records(index=False))

    del edgelist
    gc.collect()

    tf = time.time()
    time_length = (tf-ti)/60 # unit = minutes
    print("Done! Elapsed time: %f minutes" %time_length) # takes 8.2 minutes to run
    return G


# In[5]:

def topologically_sort_graph(G):
    if not G.is_dag():
        print('Graph is not DAG. If you continue the topological sorting algorithm will be stuck in an endless search. Remove cycles before to continue!')
        raise ValueError
    layers = -1 * ones(G.vcount())
    nodes_in_this_layer = where(array(G.indegree())==0)[0]
    layer = 0

    while nodes_in_this_layer.any():
        layers[nodes_in_this_layer] = layer
        layer += 1
        nodes_in_this_layer = G.neighborhood(vertices=nodes_in_this_layer.tolist(), order=1, mode='OUT')
        nodes_in_this_layer = unique([item for sublist in nodes_in_this_layer for item in sublist[1:]])
    return layers


# In[6]:

def search_path_count_of_graph(G, mode='IN', layer_name='layer'):
    layers = G.vs[layer_name]
    if mode=='IN':
        layer_values = arange(2,max(layers)+1)
    elif mode=='OUT':
        layer_values = arange(max(layers)-2, -1, -1)
    count_paths = array(G.degree(mode=mode))
    for layer in layer_values:
        for n in where(layers==layer)[0]:
            neighbors = G.neighbors(n, mode=mode)
            if neighbors:
            #Each node's count of incoming paths is the sum of its predecessors' count of incoming paths
                count_paths[n] += sum(count_paths[array(neighbors)])
    return count_paths


# In[7]:

### Define functions to generate random controls
def randomize_citations(citations,
                        patent_attributes):
    citations_randomized = citations.copy()
    not_same_year = citations_randomized['Year_Citing_Patent']!=citations_randomized['Year_Cited_Patent']
    ### Take the same-class citations of every class and permute them.
    print("Randomizing Same-Class Citations")
    same_class_ind = citations_randomized['Class_Citing_Patent']==citations_randomized['Class_Cited_Patent']
    cross_class_ind = -same_class_ind 
    same_class_ind = same_class_ind & not_same_year
    grouper = citations_randomized.ix[same_class_ind].groupby(['Year_Citing_Patent', 
                                                               'Year_Cited_Patent', 
                                                               'Class_Citing_Patent', 
                                                              ])[['Citing_Patent', 
                                                                  'Cited_Patent']]
    print("%i groups"%(len(grouper)))
    print("%i groups that can be rewired"%(sum(grouper.size()>1)))
    g = grouper.apply(randomize_citations_helper)
#     g.index = g.index.droplevel(['Year_Citing_Patent','Year_Cited_Patent','Class_Citing_Patent'])

    citations_randomized.ix[same_class_ind, ['Citing_Patent', 
                                             'Cited_Patent']
                            ] = g

    ### Take the cross-class citations and permute them.
    print("Randomizing Cross-Class Citations")        
    cross_class_ind = cross_class_ind & not_same_year
    grouper = citations_randomized.ix[cross_class_ind].groupby(['Year_Citing_Patent', 
                                                               'Year_Cited_Patent', 
                                                              ])[['Citing_Patent', 
                                                                  'Cited_Patent']]
    print("%i groups"%(len(grouper)))
    print("%i groups that can be rewired"%(sum(grouper.size()>1)))
    g = grouper.apply(randomize_citations_helper)
#     g.index = g.index.droplevel(['Year_Citing_Patent','Year_Cited_Patent'])

    citations_randomized.ix[cross_class_ind, ['Citing_Patent', 
                                             'Cited_Patent']
                            ] = g
    
    ### Drop patent attributes (which are now inaccurate for both the citing and cited patent) and bring them in from patent_attributes
    citations_randomized.drop(['Class_Citing_Patent', 'Class_Cited_Patent'], axis=1, inplace=True)
#     citations_randomized = citations_randomized[['Citing_Patent', 'Cited_Patent', 'Same_Class']]

    patent_attributes = patent_attributes[['patent_number', 'Class']].set_index('patent_number')
    citations_randomized = citations_randomized.merge(patent_attributes, 
                    left_on='Citing_Patent', 
                    right_index=True,
                    )

    citations_randomized = citations_randomized.merge(patent_attributes, 
                    left_on='Cited_Patent', 
                    right_index=True,
                    suffixes=('_Citing_Patent','_Cited_Patent'))
    return citations_randomized


def randomize_citations_helper(citing_cited):

#     if all(citing_cited['Year_Citing_Patent']==citing_cited['Year_Cited_Patent']):
#         return citing_cited[['Citing_Patent', 'Cited_Patent']]
    n_Citing = citing_cited.Citing_Patent.nunique()
    n_Cited = citing_cited.Cited_Patent.nunique()
    
    if n_Cited*n_Citing==citing_cited.shape[0]: #The graph is fully connected, and so can't be rewired
        return citing_cited#[['Citing_Patent', 'Cited_Patent']]
    
#     Citing_lookup = pd.Series(index=citing_cited.Citing_Patent.unique(),
#                               data=1+arange(n_Citing))
#     Cited_lookup = pd.Series(index=citing_cited.Cited_Patent.unique(),
#                              data=1+arange(n_Cited))
#     input_to_Birewire = array([Citing_lookup.ix[citing_cited.Citing_Patent].values,
#                                Cited_lookup.ix[citing_cited.Cited_Patent].values + n_Citing]).T
    citing_lookup = citing_cited['Citing_Patent'].astype('category')
    cited_lookup = citing_cited['Cited_Patent'].astype('category')
    input_to_Birewire = array([citing_lookup.cat.codes.values.astype('uint64'),
                               cited_lookup.cat.codes.values.astype('uint64') + n_Citing]).T+1
#     citing_cited.Citing_Patent = Citing_lookup.ix[citing_cited.Citing_Patent].values
#     citing_cited.Cited_Patent = Cited_lookup.ix[citing_cited.Cited_Patent].values
#     citing_cited.Cited_Patent += n_Citing
    import BiRewire
    this_rewiring = BiRewire.Rewiring(data=input_to_Birewire,
                               type_of_array='edgelist_b',
                               type_of_graph='bipartite')
    this_rewiring.rewire(verbose=0)   
    z = this_rewiring.data_rewired-1


#     Citing_lookup = pd.DataFrame(Citing_lookup).reset_index().set_index(0)
#     Cited_lookup = pd.DataFrame(Cited_lookup).reset_index().set_index(0)
#     citing_patents = Citing_lookup.ix[z[:,0]].values.flatten()
#     cited_patents = Cited_lookup.ix[z[:,1]-n_Citing].values.flatten()
    
    citing_patents = citing_lookup.cat.categories.values[z[:,0]]
    cited_patents = cited_lookup.cat.categories.values[z[:,1]-n_Citing]

    rewired_output = pd.DataFrame(index=citing_cited.index,
                                 columns=['Citing_Patent', 'Cited_Patent']
                                  )
    rewired_output['Citing_Patent'] = citing_patents
    rewired_output['Cited_Patent'] = cited_patents
    return rewired_output


# # Load Data 

# In[8]:

patents = pd.read_hdf(data_directory+'patents.h5', 'df')
citations = pd.read_hdf(data_directory+'citations.h5', 'df')


# Randomize Citations
# ===
# If requested

# In[9]:

if randomized_control:
    ti = time.time()
    citations = randomize_citations(citations, patents)
    tf = time.time()
    final_time_length = tf-ti
    print('Done! Randomizing citations took: %f seconds' %final_time_length + '= %f minutes' %(final_time_length/60))


# Compute SPNP
# ===

# In[10]:

G = create_patent_citation_graph(patents, citations)
if not G.is_dag():
    print('Graph is not DAG. If you continue the topological sorting algorithm will be stuck in an endless search. Remove cycles before to continue!')
    raise ValueError
print('Graph is DAG, you can continue!')


# In[11]:

layers = topologically_sort_graph(G)
G.vs['layer'] = layers
print("%i layers"%max(layers))


# In[12]:

count_incoming_paths = array(search_path_count_of_graph(G, mode='IN'))
G.vs['count_incoming_paths'] = count_incoming_paths


# In[13]:

vc = pd.value_counts(G.vs['filing_year'])
n_rows = vc.sort_index().cumsum().ix[1975:].sum()
DF = pd.DataFrame(index=arange(n_rows), 
                  columns=["patent_number",
                           "observation_year",
                          "SPNP_count",
                           "outgoing_path_count_log",
                           "incoming_path_count_log",
                          ]
#                   dtype='uint64'
                 )
size_in_GBs = (prod(DF.shape)*64)*1.25e-10
print("Memory allocated to DF_node_SPNP_over_time: %f GBs" %size_in_GBs)
   
year_list = arange(1975,max(G.vs['filing_year'])+1)
this_year_data_start = 0
for observation_year in tqdm(year_list):
    print(observation_year)
    patents_within_this_year = G.vs.select(filing_year_le=observation_year).indices
    G_subgraph = G.subgraph(patents_within_this_year, 
                            implementation="auto")
    n_row = G_subgraph.vcount()
    
    DF.ix[this_year_data_start:this_year_data_start+n_row-1, 
          'outgoing_path_count_log'] = log(search_path_count_of_graph(G_subgraph, mode='OUT')+1)
    DF.ix[this_year_data_start:this_year_data_start+n_row-1, 
          'incoming_path_count_log'] = log(count_incoming_paths[patents_within_this_year]+1)
#     DF.ix[this_year_data_start:this_year_data_start+n_row-1, 
#           'SPNP_count'] = (log(search_path_count_of_graph(G_subgraph, mode='OUT')+1) +
#                            log(count_incoming_paths[patents_within_this_year]+1)
#                            )
    DF.ix[this_year_data_start:this_year_data_start+n_row-1, 
          'patent_number'] = G_subgraph.vs['patent_number']
    DF.ix[this_year_data_start:this_year_data_start+n_row-1, 
          'observation_year'] = observation_year
    del G_subgraph
    gc.collect()
    this_year_data_start += n_row
DF['observation_year'] = DF['observation_year'].astype('uint16')
DF['patent_number'] = DF['patent_number'].astype('uint32')
DF['outgoing_path_count_log'] = DF['outgoing_path_count_log'].astype('float64')
DF['incoming_path_count_log'] = DF['incoming_path_count_log'].astype('float64')
DF['SPNP_count'] = DF['outgoing_path_count_log']+DF['incoming_path_count_log']


# # Report patents' centrality in t+2, t+3, t+5 and t+8

# In[14]:

DF_patents = DF[DF['observation_year']==2015]
DF_patents.drop(['observation_year'], axis=1, inplace=True)
DF_patents['filing_year'] = G.vs['filing_year']


# In[15]:

ti = time.time()
# for patents filed AFTER 1975 report their centrality 3/5/8 years after filing
# for patents filed BEFORE 1975 report their centrality in 1978/1980/1983. To quickly do this create a "fake filing year"
DF_patents['fake_filing_year'] = DF_patents['filing_year']
DF_patents.ix[DF_patents['filing_year']<1975, 'fake_filing_year']=1975

for horizon in tqdm([2,3,5,8]):
    DF_patents['filing_year+%i'%horizon] = DF_patents['fake_filing_year']+horizon

    # merge to add SPNP after $horizon years from filing
    DF_patents = pd.merge(DF_patents, DF.drop(['incoming_path_count_log','outgoing_path_count_log'], axis=1), 
                          how='left', left_on=['patent_number','filing_year+%i'%horizon], 
                          right_on=['patent_number','observation_year'],
                           suffixes=('', '_filing+%i'%horizon)
                           )
#                            sort=True,, copy=True, indicator=False)
    del DF_patents['filing_year+%i'%horizon]
    del DF_patents['observation_year']
    
    gc.collect()
#     DF_patents.rename(columns={"SPNP_count": 'SPNP_count_t+%i'%horizon}, inplace=True)

del DF_patents['fake_filing_year']
DF_patents.rename(columns={'SPNP_count': 'SPNP_count_2015'}, inplace=True)
DF_patents.rename(columns={'outgoing_path_count_log': 'outgoing_path_count_log_2015'}, inplace=True)
# DF_patents.rename(columns={'incoming_path_count_log': 'incoming_path_count_log_2015'}, inplace=True)

tf = time.time()
final_time_length = tf-ti
print('Done! Reporting patent centrality after x years took: %f seconds' %final_time_length + '= %f minutes' %(final_time_length/60))


# # Compute centrality of cited patents in t-1

# In[16]:

citations = pd.merge(citations, DF,
                      how='left', on=None, left_on=['Cited_Patent','Year_Citing_Patent'], 
                      right_on=['patent_number','observation_year'],
                      )
del citations['patent_number']
del citations['observation_year']
gc.collect()
citations.rename(columns={"SPNP_count": 'SPNP_count_cited_year_of_citation'}, inplace=True)

citations_grouped_by_citing = citations[['SPNP_count_cited_year_of_citation',
                                          'Citing_Patent']].groupby(['Citing_Patent']).agg(['mean',
                                                                                            'std'])

DF_patents = pd.merge(DF_patents, citations_grouped_by_citing['SPNP_count_cited_year_of_citation'], 
                      how='left', left_on='patent_number', 
                      right_index=True)

DF_patents.rename(columns={"mean": 'meanSPNPcited',
                          "std": 'stdSPNPcited'}, inplace=True)


# Some patents are cited the same year they were filed. For these patents we have no information on their SPNP the year before they were cited. We need to compute their SPNP the moment before they were cited. This is done by multiplying the number of incoming paths of the cited patent by the number of citations received in the year they were filed. This is only an approximation of the number of their outgoing paths, but it is not a bad one because they are unlikely to be cited by other patents granted in the same year that have received citations themselves.

# In[17]:

citations['Year_Citing_Patent-1'] = citations['Year_Citing_Patent'] - 1

citations = pd.merge(citations, DF.drop(['incoming_path_count_log', 'outgoing_path_count_log'], axis=1),
                      how='left', on=None, left_on=['Cited_Patent','Year_Citing_Patent-1'], 
                      right_on=['patent_number','observation_year'],
                      )
del citations['patent_number']
del citations['observation_year']
del citations['Year_Citing_Patent-1']

citations.rename(columns={"SPNP_count": 'SPNP_count_cited_1year_before_citation'}, inplace=True)

#### some patents are cited the same year they were filed. For these patents we have no information on their SPNP the year before they were cited. We need to compute their SPNP the moment before they were cited. This is done by multiplying the number of incoming paths of the cited patent by the number of citations received in the year they were filed. This is only an approximation of the number of their outgoing paths, but it is not a bad one because they are unlikely to be cited by other patents granted in the same year that have received citations themselves.

citations_same_year_ind = citations['Year_Citing_Patent']==citations['Year_Cited_Patent']
citations_same_year_count = citations[citations_same_year_ind].groupby(['Cited_Patent']).size()
citations_same_year_count.name = 'citations_at_zero'
citations_same_year = pd.DataFrame(citations_same_year_count)

citations_same_year = pd.merge(citations_same_year, 
                      pd.DataFrame({'count_incoming_paths':G.vs['count_incoming_paths']},
                                  index=DF_patents['patent_number']), 
                      how='left', left_index=True, 
                      right_index=True
                     )

citations_same_year['SPNP_at_Year_Cited_Patent'] = (log(citations_same_year['count_incoming_paths']+1) +
                                                    log(citations_same_year['citations_at_zero']-1+1)
                                                    )

citations.ix[citations_same_year_ind, 
             'SPNP_count_cited_1year_before_citation'] = citations_same_year.ix[citations.ix[citations_same_year_ind,
                                                                                             'Cited_Patent'],
                                                                                'SPNP_at_Year_Cited_Patent'
                                                                               ].values
del citations_same_year_ind, citations_same_year, citations_same_year_count


# In[18]:

citations_grouped_by_citing = citations[['SPNP_count_cited_1year_before_citation',
                                          'Citing_Patent']].groupby(['Citing_Patent']).agg(['mean',
                                                                                            'std'])

DF_patents = pd.merge(DF_patents, citations_grouped_by_citing['SPNP_count_cited_1year_before_citation'], 
                      how='left', left_on='patent_number', 
                      right_index=True)

DF_patents.rename(columns={"mean": 'meanSPNPcited_1year_before',
                          "std": 'stdSPNPcited_1year_before'}, inplace=True)


# In[19]:

citations.rename(columns={"incoming_path_count_log": 'incoming_path_count_log_cited'}, inplace=True)

citations_grouped_by_citing = citations[['incoming_path_count_log_cited',
                                          'Citing_Patent']].groupby(['Citing_Patent']).agg(['mean',
                                                                                            'std'])

DF_patents = pd.merge(DF_patents, citations_grouped_by_citing['incoming_path_count_log_cited'], 
                      how='left', left_on='patent_number', 
                      right_index=True)

DF_patents.rename(columns={"mean": 'mean_incoming_path_count_log_cited',
                          "std": 'std_incoming_path_count_log_cited'}, inplace=True)


# Store data
# ===

# In[20]:

DF_patents['filing_year'] = DF_patents['filing_year'].astype('uint16')
DF_patents['patent_number'] = DF_patents['patent_number'].astype('uint32')

for c in DF_patents.columns:
#     if c.startswith('SPNP_count'):
#         DF_patents[c] = DF_patents[c].fillna(0).astype('uint64')
    if DF_patents[c].dtype =='float':
        DF_patents[c] = DF_patents[c].astype('float32')


# In[21]:

if not randomized_control:
    DF_patents.to_hdf(data_directory+'centralities/empirical.h5', 'df', complevel=9, complib='blosc')
else:
    DF_patents.to_hdf(data_directory+'centralities/controls/%s/synthetic_control_%i.h5'%(class_system,
                                                                                         randomization_id), 'df', complevel=9, complib='blosc')


# In[22]:

final_time_length = time.time()-start_time

print('Done! Total job took: %f seconds' %final_time_length + '= %f minutes' %(final_time_length/60))

