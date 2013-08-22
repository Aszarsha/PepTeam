#!/usr/bin/python

#import rpy2,
import os
import csv
import sys
import decimal 
import numpy
from operator import itemgetter
import scipy
import getopt

from scipy.stats import norm
try:
	opts, args = getopt.getopt(sys.argv[1:], "hi:v", ["help", "input="])
except getopt.GetoptError as err:
	# print help information and exit:
	print str(err) # will print something like "option -a not recognized"
	#usage()
	print "usage program --input filename"
	sys.exit(2)
input_file = None
verbose = False
for i, a in opts:
	if i == "-v":
		verbose = True
	elif i in ("-h", "--help"):
		#usage()
		print "usage program --input filename"
		sys.exit()
	elif i in ("-i", "--input"):
		input_file = a
	else:
		assert False, "unhandled option"


def robust_int(s):
	try:
		return int(s)
	except ValueError:
		return s
#input_file = "~/Users/benjamindartigues/Projet_Pepteam/ps05.mapping.prot2pep.t40.scoring.mappings"
#input_file = "./fragments_100.fasta.fastIdx.pepTree.7.mapping.1_0.profiles"

#input_file="~/Downloads/1373548035-KVS-EAE/KVS_EAE_9Core_AR.mapping.prot2pep.scoring.mappings"


PVALTHR=0.05

data= [map(robust_int,x.strip().split()) for x in open(os.path.expanduser(input_file),"r").readlines()]


#mart_annot=[x.strip().split("\t") for x in open("mart_annot.tsv")][1:]
#mart_annot=dict([(x[-1],x) for x in mart_annot])

def display_coverage(keys):
	if type(keys)==str:
		keys=[keys]
	for k in keys:
		print k,"\t"," ".join(map(str,data[k]))


data=dict([(x[0],x[1:]) for x in data])


# Avg prot length is 451 aa 
lengths=dict([(k,len(v)) for k,v in data.items()])


# Decide if we use a running sum of non zero elements or individual stacks 
# Compute the overall average 

all_aucs=[]
for row in data.values() :
	all_aucs.extend(row[1:])

avg,std = scipy.average(all_aucs),scipy.std(all_aucs)


# For each protein, compute the sum and the corresponding Z-scores 
z_scores = {}
all_sums={}
p_values={}

for k,v in data.items() : 
	n_aa = len(v)  # we remove 1 as the prot ID is in the first position 
	this_sum = sum(v)
	this_z = (this_sum - n_aa*avg) / (scipy.sqrt(n_aa)*std)
	z_scores[k]=this_z
	all_sums[k]=this_sum
	p_values[k]=1-norm.cdf(this_z)


# How many significant at PVALTHR ?
len([x for x in p_values.values() if x<=PVALTHR])

# 
z_scores_sorted =sorted(z_scores.items(),key=itemgetter(1)) 

# Minimal p-value ? 
p_values_sorted=sorted(p_values.items(),key=itemgetter(1))

# Perform multiple testing correction: 7127 significant out of 50k 
len([x for x in p_values.values() if x<=PVALTHR/len(p_values)])


# Generate output table 
f=open("significance.tsv","w")
thr=PVALTHR/len(data)
for prot,pval in [(k,v) for k,v in p_values_sorted if v<=thr]:
	#f.write("%s\t%e\t%.4f\t%s\n"%(prot,pval,z_scores[prot],"\t".join(mart_annot.get(prot,[]))))
	f.write("%s\t%e\t%.4f\n"%(prot,pval,z_scores[prot]))
f.close()


# In CSV format 
with open('all_significant_prots.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
	thr=PVALTHR/len(data)
	for prot,pval in [(k,v) for k,v in p_values_sorted if v<=thr]:
		#writer.writerow([prot,pval,z_scores[prot]]+mart_annot.get(prot,[]))
		writer.writerow([prot,pval,z_scores[prot]])
