#!/usr/bin/env python
# coding: utf-8

# 
# Introduction
# 
# Scientific Question:
# Cancer is a very prominent disease among mammals; however, since its more relevant in humans in research does this directly 
# translate to humans having more severe/ expressed levels of cancer or do other animals equivalently or express greater 
# levels of expression within cancer cells when sequencing their P53 genes ?
# 
# 
# If P53 is a tumor suppresor gene that is highly expressed in mammals then I would like to sequence both human and whale P53 
# genes to their respective reference genome to see which animal has the cancer gene more highly expressed since they are 
# usually only studied in the context of humans .
# 
# Hypothesis:
# My hypothesis is that through sequencing both humans and whale P53 genes to their reference genome we will get similar results
# on the expression of the cancer levels if not whales having even higher expressions because they have no means of combatting
# against the disease. 
# 
# The data I am using is fasta files from both humans and beluga whales sequenced P53 gene that are from the ncbi website,
# where I will use pairwise sequencing to sequence both animals to the reference gene and compare their differentially expressed
# levels of P53 gene 
# 
# 
# 

# In[1]:


'''
Packages 

pandas : high level fast and flexible data structure package that can be used to work with high dimensional data sets such 
as analyzing fastqc files

numpy: numpy is used for a multitude of mathmatical expressions for working with matrices and arrays which i used for when 
creating the graphs 

matplotlib: python package for plotting library, i used it to visualize my alignment scores in the form of a bar graph

Biopython: is open software for computational analysis of biological data sets from DNA and RNA, which i imported for pairwise
sequencing

Pairwise2: is used for multiple sequence alignment and to generate a cummulative score for the matching sequences between the 
two alignments in which I used for comparing generating the human and whale scores to compare 



'''

import os
import pandas as pd
import numpy as np
import scipy as scip
import matplotlib.pyplot as plt
import Bio
from Bio.Seq import Seq
from Bio import pairwise2
from Bio import SeqIO


# In[2]:


#Loading in Data

#reading in the human dataset
human_file = open('homo_sapien1', 'r')
human_data = ''
#iterating through every line in the fasta file
for line in human_file:
    #joining every line to be read out 
    human_data += ''.join(line)
    
#print(''.join(human_data))

#reading in the human dataset
human_file2 = open('homo_sapien2', 'r')
human_data2 = ''
#iterating through every line in the fasta file
for line in human_file2:
    #joining every line to be read out 
    human_data2 += ''.join(line)
    
#print(''.join(human_data2))

#reading in the animal dataset
animal_file = open('whale_fasta', 'r')
animal_data = ''
#iterating through every line in the fasta file
for line in animal_file:
    #joining every line to be read out     
    animal_data = ''.join(line)

#print(''.join(animal_data))


    


# In[ ]:





# In[3]:


#read in the cancer sample file
seq1 = SeqIO.read("homo_sapien1", "fasta")
#read in the reference data file
seq2 = SeqIO.read("homo_sapien_refseq", "fasta")
#perform alignment with the sample and reference
alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)
#store just the score of the alignment with the alignment_score 
alignments_score = pairwise2.align.globalxx(seq1.seq, seq2.seq, score_only=True)

#print out results
print(alignments_score)
print(pairwise2.format_alignment(*alignments[0])) 


# In[ ]:





# In[4]:


#read in the cancer sample file
seq1 = SeqIO.read("homo_sapien2", "fasta")
#read in the reference data file
seq2 = SeqIO.read("homo_sapien_refseq", "fasta")
#perform alignment with the sample and reference
alignments2 = pairwise2.align.globalxx(seq1.seq, seq2.seq)
#store just the score of the alignment with the alignment_score 
alignments2_score = pairwise2.align.globalxx(seq1.seq, seq2.seq, score_only=True)

#print out results
print(alignments2_score)
print(pairwise2.format_alignment(*alignments2[0])) 


# In[ ]:





# In[5]:


#read in the whale cancer sample file
seq1 = SeqIO.read("whale_fasta", "fasta")
#read in the whale reference data file
seq2 = SeqIO.read("whale_refseq", "fasta")
#perform alignment with the sample and reference
alignments3 = pairwise2.align.globalxx(seq1.seq, seq2.seq)
#store just the score of the whale alignment with the alignment_score 
alignments3_score = pairwise2.align.globalxx(seq1.seq, seq2.seq, score_only=True)

#print out results
print(alignments3_score)
print(pairwise2.format_alignment(*alignments3[0])) 


# In[ ]:





# In[6]:


#store all the alignment scores in a list
alignment_scores = [alignments_score, alignments2_score, alignments3_score]
#name for each score
species = ['homo_sapien1', 'homo_sapien2', 'whale']
#plot the bar with the scores and species
plt.bar(species, alignment_scores)
plt.xticks(species)
#add extra values for reference and easier to visualize
alignment_scores += [0, 1000, 2000]
#labeling the x and y axis
plt.yticks(alignment_scores) 
plt.xlabel('Species')
plt.ylabel('Alignment Scores')


# Results:
# 
# Through comparing the P53 samples of human and whale genes, through the bioinformatics method of pairwise sequencing, where I compared the cancer sample of each species with a refence genome, the data provide adequate findings that display how humans on average display higher levels of the P53 gene and in turn have higher cancer gene expression levels as compared to whales. 
# 
# However the data is not very evident since the results for the levels of gene expression are relatively close, but still display higher expression within the human samples, also a draw back was that it was very difficult to find cancer samples of whales so the scores aren't as normalized as desired. 
# 
# References:
# https://towardsdatascience.com/building-bar-charts-using-matplotlib-c7cf6db3e728
# 
# https://biopython.org/docs/1.75/api/Bio.Align.Applications.html
# 
