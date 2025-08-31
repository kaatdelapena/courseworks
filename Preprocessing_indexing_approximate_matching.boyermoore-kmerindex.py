#!/usr/bin/env python
# coding: utf-8

# In[16]:


get_ipython().system('wget http://d28rh4a8wq0iu5.cloudfront.net/ads1/code/bm_preproc.py')


# In[ ]:





# In[17]:


get_ipython().system('wget http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/chr1.GRCh38.excerpt.fasta')


# In[ ]:





# In[18]:


def naive(p, t):
    occurrences = []
    comparisons = 0
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)): 
            comparisons += 1      #counts comparison
            if t[i+j] != p[j]: 
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences, comparisons


# In[20]:


def boyer_moore_with_counts(p, p_bm, t):
    i = 0
    occurrences = []
    comparisons = 0
    alignments = 0

    while i <= len(t) - len(p):
        alignments += 1
        shift = 1
        mismatched = False
        for j in range(len(p) - 1, -1, -1):
            comparisons += 1
            if p[j] != t[i + j]:
                skip_bc = p_bm.bad_character_rule(j, t[i + j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift

    return occurrences, comparisons, alignments


# In[21]:


def read_fasta_sequence(filename):
    with open(filename) as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence.upper()
t = read_fasta_sequence("chr1.GRCh38.excerpt.fasta")


# In[22]:


#How many alignments and comparisons does the naive exact matching algorithm try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG to the excerpt of human chromosome 1? 
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
t = read_fasta_sequence("chr1.GRCh38.excerpt.fasta")

matches, comparisons = naive(p, t)

print(f"Total matches found: {len(matches)}")
print(f"Number of alignments attempted: {len(t) - len(p) + 1}")
print(f"Total character comparisons: {comparisons}")


# In[23]:


#How many alignments does the Boyer Moore try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG to the excerpt of human chromosome 1?
from bm_preproc import BoyerMoore

p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
t = read_fasta_sequence("chr1.GRCh38.excerpt.fasta")

p_bm = BoyerMoore(p, alphabet='ACGT')
matches, alignments, comparisons = boyer_moore_with_counts(p, p_bm, t)
                                                    
print(f"Matches found: {len(matches)}")
print(f"Alignments attempted: {alignments}")
print(f"Character comparisons: {comparisons}")


# In[24]:


get_ipython().system('wget https://d28rh4a8wq0iu5.cloudfront.net/ads1/code/kmer_index.py')


# In[25]:


def build_kmer_index(text, k):
    index = {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        if kmer in index:
            index[kmer].append(i)
        else:
            index[kmer] = [i]
    return index


# In[26]:


def approximate_match_pigeonhole(pattern, text, max_mismatches, k):
    segment_length = len(pattern) // (max_mismatches + 1)
    matches = set()
    index = build_kmer_index(text, k)

    for i in range(max_mismatches + 1):
        start = i * segment_length
        end = start + k
        kmer = pattern[start:end]

        if kmer in index:
            for offset in index[kmer]:
                # Align pattern so that k-mer matches at correct position
                match_start = offset - start
                if match_start < 0 or match_start + len(pattern) > len(text):
                    continue

                mismatches = 0
                for j in range(len(pattern)):
                    if pattern[j] != text[match_start + j]:
                        mismatches += 1
                        if mismatches > max_mismatches:
                            break
                if mismatches <= max_mismatches:
                    matches.add(match_start)
    return sorted(matches)


# In[27]:


p = 'GGCGCGGTGGCTCACGCCTGTAAT'
t= read_fasta_sequence('chr1.GRCh38.excerpt.fasta')
hits = approximate_match_pigeonhole(p, t, max_mismatches=2, k=8)

print(f"Total approximate matches: {len(hits)}")
print("Offsets:", hits[:10])


# In[28]:


from kmer_index import Index


# In[29]:


index = Index(t, k=8)


# In[33]:


#How many times does the string GGCGCGGTGGCTCACGCCTGTAAT occur and how many index hits with up to 2 substitutions in the excerpt of human chromosome 1?

p = 'GGCGCGGTGGCTCACGCCTGTAAT'
segment_length = len(p) // 3
total_hits = 0

for i in range(3):
    segment = p[i * segment_length: (i+1)*segment_length]
    hits = index.query(segment)
    total_hits += len(hits)

print(f"Segment {i}: {segment}, Hits: {len(hits)}")
print(f"Total index hits across all segments: {total_hits}")


# In[37]:


#Hint 1: Multiple index hits might direct you to the same match multiple times, but be careful not to count a match more than once.

segment_length = 8
candidates = set()
segments = [0, 8, 16] 

for i, start in enumerate(segments):
    kmer = p[start:start + segment_length]
    hits = index.query(kmer)
    for hit in hits:
        candidate_offset = hit - start
        if 0 <= candidate_offset <= len(t) - len(p):
            candidates.add(candidate_offset)
            
matches = []
for offset in candidates:
    mismatches = sum(1 for j in range(len(p)) if p[j] != t[offset + j])
    if mismatches <= 2:
        matches.append(offset)

print(f"Final matches with ≤2 mismatches: {len(matches)}")

