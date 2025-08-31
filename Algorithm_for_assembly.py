#!/usr/bin/env python
# coding: utf-8

# In[1]:


def overlap(a, b, min_length=3):
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


# In[2]:


import itertools

def scs(ss):
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest


# In[3]:


def scs(ss, min_length=1):
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]
        for i in range(len(ssperm) - 1):
            olen = overlap(ssperm[i], ssperm[i+1], min_length)
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup
    return shortest_sup


# In[4]:


#What is the length of the shortest common superstring of the following strings?
reads = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
result = scs(reads, min_length=1)
print(result)
print(len(result))


# In[35]:


#How many different shortest common superstrings are there for the input strings given in the previous question?
def all_shortest_scs(ss, min_length=1):
    shortest_sup = None
    shortest_set = set()
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]
        for i in range(len(ssperm) - 1):
            olen = overlap(ssperm[i], ssperm[i+1], min_length)
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup
            shortest_set = {sup}
        elif len(sup) == len(shortest_sup):
            shortest_set.add(sup)
    return shortest_set


# In[36]:


strings = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
scs_set = all_shortest_scs(strings)
print(f"Number of distinct shortest common superstrings: {len(scs_set)}")


# In[37]:


get_ipython().system('wget https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ads1_week4_reads.fq')


# In[38]:


from collections import defaultdict

def read_fastq(filename):
    sequences = []
    with open(filename) as f:
        while True:
            f.readline()  # skip name line
            seq = f.readline().strip()
            f.readline()  # skip plus line
            f.readline()  # skip quality line
            if not seq:
                break
            sequences.append(seq)
    return sequences


# In[39]:


reads = read_fastq('ads1_week4_reads.fq.1')


# In[40]:


# De Bruijn Algorithm
def de_bruijn_ize(st, k):
    edges = []
    nodes = set()
    for i in range(len(st) - k + 1):
        edges.append((st[i:i+k-1], st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes, edges


# In[48]:


# greedy
def overlap(a, b, k):
    start = 0
    while True:
        start = a.find(b[:k], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def greedy_scs(reads, k):

    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)


# In[49]:


import itertools

def pick_maximal_overlap(reads, k):
    reada, readb = None, None
    best_olen = 0
    for a, b in itertools.permutations(reads, 2):
        olen = overlap(a, b, k=k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen


# In[50]:


reads = read_fastq('ads1_week4_reads.fq.1')


# In[ ]:


result = greedy_scs(reads.copy(), k=1)
print("Greedy SCS result:", result)


# In[ ]:


print(len(genome))


# In[ ]:


print(f"The total count of adenine {genome.count('A')} and the total count of thymine {genome.count('T')}.")


# In[ ]:




