
# coding: utf-8

# In[ ]:

# library of smer functions

# In[ ]:

import numpy as np
import pandas as pd
import h5py
from itertools import product

# data processing

# In[ ]:

def load_mat_file(smer_dir):
    '''load .mat v7.3'''
    #https://stackoverflow.com/questions/19310808/how-to-read-a-v7-3-mat-file-via-h5py
    cdist_list = ['Cdist3.mat','Cdist5.mat','Cdist7.mat','Cdist9.mat']
    remove_last_ind = [0,1,1,0]
    fnames = [smer_dir+cdist for cdist in cdist_list]
    mat = []
    k=0
    for fname in fnames:
        f = h5py.File(fname,'r')
        data = f[list(f.keys())[1]]   
        tmp = []
        i=-1
        for i in range(len(data)-remove_last_ind[k]):
            tmp.append(f[data[i,0]].value)
        mat.append(tmp)
        k+=1
    return mat

# In[ ]:

def load_smer_data(direc = "Data/"):
    '''load smer data'''
    fn = []
    idx_list = ['k3_idx.csv','k5_idx.csv','k7_idx.csv','k9_idx.csv']
    fn = [direc+idxl for idxl in idx_list]
    
    #make list of dictionaries of lists
    di = list()
    for i in range(0,4):
        data = pd.read_csv(fn[i],sep=',',header=None)
        data.set_index([0], drop=True, inplace=True)
        ditmp=[]
        for ii in range(0,data.shape[1]):
            ditmp.append(data.loc[:,ii+1].to_dict())
        di.append(ditmp)
    
    return di, load_mat_file(direc)


# distance functions

# In[ ]:

def seq2int(seq):
    '''sequence to integer'''
    alphabet = ['A','C','G','T','N']
    return np.array([np.where(np.isin(alphabet,x))[0][0] for x in seq])

# In[ ]:

def hamming(in1,in2):
    '''p-distance = hamming/length'''
    return np.count_nonzero(seq2int(in1)-seq2int(in2))

# In[ ]:

def ham(li):
    '''hamming for sequence list'''
    n = len(li)
    mat = []
    for i in range(n):
        v =[]
        for j in range(n):
            d = hamming(li[i],li[j])
            v.append(d)
        mat.append(v)
    return mat

# In[ ]:

def sdist(c1,c2,dm):
    '''get s-distance from precomputed distance matrix'''
    d = 0
    for i in range(len(c1)):
        for j in range(len(c2)):
            d += dm[int(c1[i])-1][int(c2[j])-1]
    return d

# In[ ]:

def sdist_list(li,dm):
    '''s-distance for s-code list'''
    n = len(li)
    mat = []
    for i in range(n):
        v =[]
        for j in range(n):
            d = sdist(li[i],li[j],dm)
            v.append(d) 
        mat.append(v)
    return mat


# basic functions

# In[ ]:

def get_smer_count(kmer_count,dic,k=5,c=64):
    '''kmer counts to smer counts'''
    
    # prerequisites
    k_conversion = np.array([3,5,7,9])
    clust_conversion = np.array([4,8,16,32,64,128,256])
    kmer = np.where(k_conversion == k)[0][0]
    clust = np.where(clust_conversion == c)[0][0]
    
    # smer indices
    idx = list(dic[kmer][clust].values())
    perms = np.unique(idx)
    smer_count = np.zeros((len(kmer_count),len(perms)))
    
    # can probably be vectorized
    for i in range(len(kmer_count)):
        for ii in range(len(kmer_count[0])):
            smer_count[i][idx[ii]-1] += kmer_count[i][ii]
    
    return smer_count    

# In[ ]:

def get_k_from_i(i):
    return 2*(i+1)+1

# In[ ]:

def get_clust_from_i(i):
    return 2**(i+2)

# In[ ]:

def get_i_from_k(k):
    return int((k-1)/2-1)

# In[ ]:

def get_i_from_clust(clust):
    return int(np.log2(clust)-2)

# In[ ]:

def colors(seq,dic,k):
    '''sequence to s-code'''
    cz = []
    for i in range(len(seq)-k+1):
        x = seq[i:i+k]
        color = dic[x]
        cz.append(color)
    return cz

# In[ ]:

def find_initial_points(dic,colo): 
    '''find initial points'''
    intu=[]
    for x in dic:
        if dic[x] == colo:
            intu.append(x)
    return intu

# In[ ]:

def reconstructions(seq,dic,k):
    '''reconstruct sequence variants from s-code'''
    chars = "ACGT"
    intu = find_initial_points(dic,seq[0])
    nexttu = []
    for i in range(1,len(seq)):
        newc = seq[i]
        for tup in intu:
            x = tup[-k:]
            for c in chars:
                y = x[1:]+c
                w = dic[y]
                if w == newc:
                    newtup = tup+c
                    nexttu.append(newtup)
        intu = nexttu[:]
        nexttu = []
    return intu

# In[ ]:

def variants(seq,di,k):
    '''variants for instance of sequence and s-code parameters'''
    laz = colors(seq,di,k)
    re = reconstructions(laz,di,k)
    return re


# other functions

# In[ ]:

def read_fasta_seqs(fname):
    '''read fasta to set'''
    out = []
    with open(fname, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                out.append(str(line).strip('\n'))
    return out

# In[ ]:

def sdist_seq(s1,s2,dm,dic,k):
    '''get s-distance from precomputed distance matrix'''
    c1 = colors(s1,dic,k)
    c2 = colors(s2,dic,k)
    d = 0
    for i in range(len(c1)):
        d += dm[int(c1[i])-1][int(c2[i])-1]
    return d

# In[ ]:

def sdist_kmers(s1,s2,dm,dic):
    '''s-distance for kmers'''
    return dm[int(dic[s1])-1][int(dic[s2])-1]


# kmer functions


# In[ ]:

# http://claresloggett.github.io/python_workshops/improved_kmers.html
def count_kmers(read, k):
    """Count kmer occurrences in a given read.

    Parameters
    ----------
    read : string
        A single DNA sequence.
    k : int
        The value of k for which to count kmers.

    Returns
    -------
    counts : dictionary, {'string': int}
        A dictionary of counts keyed by their individual kmers (strings
        of length k).

    Examples
    --------
    >>> count_kmers("GATGAT", 3)
    {'ATG': 1, 'GAT': 2, 'TGA': 1}
    """
    # Start with an empty dictionary
    counts = {}
    # Calculate how many kmers of length k there are
    num_kmers = len(read) - k + 1
    # Loop over the kmer start positions
    for i in range(num_kmers):
        # Slice the string to get the kmer
        kmer = read[i:i+k]
        # Add the kmer to the dictionary if it's not there
        if kmer not in counts:
            counts[kmer] = 0
        # Increment the count for this kmer
        counts[kmer] += 1
    # Return the final counts
    return counts


# In[ ]:

def count_kmersV2(read, k):
    """
    Update function to fill precomputed dictionary of all possible kmers
    Returns alphabetically sorted list of probabilites
    """
    # list
    perms = [''.join(p) for p in product('ACGT',repeat=k)]
    
    # Calculate how many kmers of length k there are
    num_kmers = len(read) - k + 1
    
    # make empty dictionary
    counts = dict(zip(perms,np.zeros(len(perms))))
    
    # Loop over the kmer start positions
    for i in range(num_kmers):
        
        # Slice the string to get the kmer
        kmer = read[i:i+k]
        counts[kmer] += 1
        
    count_vals = pd.DataFrame.from_dict(counts,orient='index').sort_index().values

    return count_vals #/num_kmers # for simplicity return counts


# In[ ]:

def pandas_expand_kmer_counts(df):
    '''expand kmer counts to include all kmer types'''
    # initialize empty count dictionary
    perms = [''.join(p) for p in product('ACGT',repeat=k)]
    counts = dict(zip(perms,np.zeros(len(perms))))
    for ind,row in df.iterrows():
        counts[ind] = int(row.values[0])
    count_vals = pd.DataFrame.from_dict(counts,orient='index').sort_index()#.values
    return count_vals


# In[ ]:

def get_pandas_smer_count(kmer_count,dic,k=5,c=64):
    '''kmer counts dframe to smer counts'''
    
    # prerequisites dic
    k_conversion = np.array([3,5,7,9])
    clust_conversion = np.array([4,8,16,32,64,128,256])
    kmer = np.where(k_conversion == k)[0][0]
    clust = np.where(clust_conversion == c)[0][0]
    dicc = dic[kmer][clust]
    
    # initialize data container
    perms = np.unique(list(dicc.values()))
    counts = dict(zip(perms,np.zeros(len(perms))))
    #smer_count = np.zeros((len(kmer_count),len(perms)))
    
    # can probably be vectorized
    for ind,row in kmer_count.iterrows():
        counts[dicc[ind]] = int(row.values[0])
    count_vals = pd.DataFrame.from_dict(counts,orient='index').sort_index()#.values
    return count_vals 


# In[ ]:

def dict_keys_from_value(dic,search_val):
    '''for s-cluster print kmers'''
    keys = []
    for key, val in dic.items():    # for name, age in dictionary.iteritems():  (for Python 2.x)
        if val == search_val:
            keys.append(key)
            #print(key)
    return keys

