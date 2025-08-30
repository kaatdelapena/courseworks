#!/usr/bin/env python
# coding: utf-8

# In[8]:


get_ipython().system('wget http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa')


# In[7]:


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
genome = readGenome('lambda_virus.fa')


# In[13]:


len(genome)


# In[14]:


counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
for base in genome:
    counts[base] += 1
print(counts)


# In[15]:


import collections
collections.Counter(genome)


# In[1]:


def naive_match(pattern, text):
    occurrences = []
    for i in range(len(text) - len(pattern) + 1):
        match = True
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences


# In[2]:


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq))


# In[4]:


pattern1 = 'AGGT'
pattern2 = reverse_complement(pattern1)


# In[7]:


count1 = len(naive_match(pattern1, genome))
count2 = len(naive_match(pattern2, genome))


# In[8]:


#How many times does AGGT or its reverse complement ACCT occur in the lambda virus genome?

print(f"Total occurrences of AGGT or ACCT: {count1 + count2}")


# In[9]:


pattern1 = 'TTAA'
pattern2 = reverse_complement(pattern1)
Count1 = len(naive_match(pattern1, genome))
Count2 = len(naive_match(pattern2, genome))


# In[14]:


#How many times does TTAA or its reverse complement occur?

print (f"Total occurences of TTAA and its reverse complement: {Count1}")


# In[15]:


#What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement in the Lambda virus genome?  

def find_leftmost(pattern, genome):
    return genome.find(pattern)
pattern = 'ACTAAGT'
rev_comp = reverse_complement(pattern)

offset_pattern = find_leftmost(pattern, genome)
offset_revcomp = find_leftmost(rev_comp, genome)

# Report the smaller offset (i.e., leftmost)
leftmost_offset = min(offset for offset in [offset_pattern, offset_revcomp] if offset != -1)
print(f"Leftmost offset: {leftmost_offset}")


# In[24]:


#What is the offset of the leftmost occurrence of AGTCGA or its reverse complement in the Lambda virus genome? 
pattern = 'AGTCGA'
rev_comp = reverse_complement(pattern)

offset_pattern = find_leftmost(pattern, genome)
offset_revcomp = find_leftmost(rev_comp, genome)

leftmost_offset = min(offset for offset in [offset_pattern, offset_revcomp] if offset !=-1)
print(f"Leftmost offset: {leftmost_offset}")


# In[10]:


# Implement the naive_2mm function
#How many times does TTCAAGCC occur in the Lambda virus genome when allowing up to 2 mismatches?

def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.strip()
    return genome

def naive_2mm(pattern, text):
    matches = []
    for i in range(len(text) - len(pattern) + 1):
        mismatches = 0
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            matches.append(i)
    return matches

pattern = 'TTCAAGCC'
hits = naive_2mm(pattern, genome)

print(f"Total matches with ≤2 mismatches: {len(hits)}")
print("Offsets of matches:")
print(hits)


# In[11]:


pattern = 'AGGAGGTT'
hits = naive_2mm(pattern, genome)

print(f"Total matches with ≤2 mismatches: {len(hits)}")
print("Offsets of matches:")
print(hits)


# In[14]:


get_ipython().system('wget  https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR037900_1.first1000.fastq')


# In[16]:


#Report which sequencing cycle has the problem -Phred scores

def read_fastq(filename):
    sequences = []
    qualities = []
    with open(filename) as f:
        while True:
            f.readline()
            seq = f.readline().strip()
            f.readline()  
            qual = f.readline().strip()
            if not qual:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


# In[17]:


def phred33_to_q(qual):
    return [ord(c) - 33 for c in qual]


# In[18]:


def find_bad_cycle(qualities):
    num_cycles = len(qualities[0])
    cycle_scores = [0] * num_cycles
    for qual in qualities:
        scores = phred33_to_q(qual)
        for i in range(num_cycles):
            cycle_scores[i] += scores[i]
    avg_scores = [total / len(qualities) for total in cycle_scores]
    return avg_scores.index(min(avg_scores))


# In[21]:


sequences, qualities = read_fastq('ERR037900_1.first1000.fastq.2')
bad_cycle = find_bad_cycle(qualities)
print(f"Problematic sequencing cycle: {bad_cycle}")


# In[ ]:




