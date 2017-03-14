##FUNCTION: plot whole genome data for NIPT project
##INPUT: reads files from specified results directory
##OUTPUT: a tabbed pdf file, one sample per tab, with a whole genome view of a specified quantity
import pickle
import matplotlib
matplotlib.use('Agg')
import argparse
import os
import matplotlib.pyplot as plt
import matplotlib.legend as legend
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

parser = argparse.ArgumentParser(description='Plotting Arguments',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('type', default='gcc', type=str,
                   help='suffix of file type to plot')

#specify input directory
parser.add_argument('ref_dir', type=str,
                   help='reference table used for within sample comparison')

args = parser.parse_args()

#get list of files in directory
lst=os.listdir(args.ref_dir)
#initialise genome array
chromsToProcess=range(1,23)

#convert file type to
#name: for pdf file naming
#label: to label axis
#What I'm plotting:
#Z-score per bin, template length (column 9 of bam file), GC corrected coverage, raw coverage (total, not mean), AS score(mean)
files=['.zscore','_2_stats.pickle','.gcc','_ruh.pickle','_qual_stats.pickle']
name=['z-score','length',,'gcc','coverage','quality']
label=['Z Score','Frag Length','Lowess Weighting','GC corrected','Coverage','AS Score']

for item in files:
        if item in args.type:
                type=item

name1=name[files.index(type)]
pos1=name[files.index(type)]

#I've specifically not used a dictionary here as the matching is not exact
dirs=['referr_77','refs_ERR_77','refs_SRR','refs_ERR','whole_batch','HG00','referr','refsrr','inhouse','Coverage','53','batch','haus']
ref_name=['iPSc_77','iPSc_77','cf','iPSc','Sample 4A1','GIAB','iPSc','cf','HD_bwa','HD_isc','Run 53','All of 4A1','HD_bwa_no_14']

#get test set from directory name
for item1 in dirs:
        if item1 in args.ref_dir.split('/')[-2]:
                print(args.ref_dir.split('/')[-2],'set')
                set=ref_name[dirs.index(item1)]
                break

#get ref set from file suffix
set_ref=''
for item1 in dirs:
        if item1 in args.type:
                set_ref=ref_name[dirs.index(item1)]
                print(set_ref,item1,'set_ref')
                break
ref='test'
if not set_ref:
        set_ref=set
        ref='REF'


samples=[]
samples_ref=dict()
print(name1,set,set_ref,ref)
print(args.ref_dir + str(name1)  + '_' + str(set) + '_' +  str(set_ref)   + '_genome_marked.pdf')
with PdfPages(args.ref_dir + str(name1)  + '_' + str(set)  +  '_' + str(set_ref)   +  '_genome_marked.pdf') as pdf:
	#get raw data files with required suffix
	for item in lst:
		if (args.type in item) & ('pdf' not in item) :
			excluded=item.split(".")[0]+'.plot'
			#format to get sample name
			A=item.replace(args.type,'')
			samples.append(A)
			#get name of file with reference bin info
			excluded=item.split(".")[0]+'.plot'
			#add entry for reference bin info to dict
			samples_ref[A]=excluded
			print(item)
#       plot title page with legend.
	ax = plt.axes()
	plt.plot(np.zeros( (len(samples)+1 ,100) ))
	plt.title(set)
	sam=[]


	for item in samples:
		#initialise new page of pdf
		plt.clf()
		plt.figure(figsize=(20,3))
		B=[]
		#initialise array to store no of bins per chr
		chromEnd=[]
		#load up raw data
		lookUpTable2 = pickle.load(open(args.ref_dir + item + args.type,'rb'))
		#initialise dicts to store data on markedBins(excluded) and blinds (no good reference bins found)
		lookUpTable0=dict()
		lookUpTable1=dict()
		print(chromsToProcess)
		length_total=0
		#initialise arrays for each chromosome - do chromosomes separately to make sure we get some no of data points
		#for each chr
		for chrom in chromsToProcess:
			lookUpTable0[str(chrom)]=[np.nan] * len(lookUpTable2[str(chrom)])
			lookUpTable1[str(chrom)]=[np.nan] * len(lookUpTable2[str(chrom)])
			length_total=length_total+len(lookUpTable1[str(chrom)])
			chromEnd.append(length_total)
		#get raw info on marked and blinds
		ref_info=samples_ref[item]
		plotbin = pickle.load(open(args.ref_dir + ref_info,'rb' ))
		A = plotbin['markedBins'];
		#read info into initialised arrays, double check they are both the sam length
		for item2 in A:
			while len(lookUpTable0[str(item2[0])]) <= item2[1]:
				lookUpTable0[str(item2[0])].append(np.nan)
			while len(lookUpTable1[str(item2[0])]) <= item2[1]:
				lookUpTable1[str(item2[0])].append(np.nan)
			while len(lookUpTable2[str(item2[0])]) <= item2[1]:
				lookUpTable2[str(item2[0])].append(np.nan)
			lookUpTable0[str(item2[0])][item2[1]] = lookUpTable2[str(item2[0])][item2[1]]
		C = plotbin['blindsDict']
		for chrom in chromsToProcess:
			for item2 in sorted(C[str(chrom)]):
				while len(lookUpTable1[str(chrom)]) <= item2:
					lookUpTable1[str(chrom)].append(np.nan)
				while len(lookUpTable0[str(chrom)]) <= item2:
					lookUpTable0[str(chrom)].append(np.nan)
				while len(lookUpTable2[str(chrom)]) <= item2:
					lookUpTable2[str(chrom)].append(np.nan)
				lookUpTable1[str(chrom)][item2] = lookUpTable2[str(chrom)][item2]
		A=[]
		B=[]
		C=[]
		#now put all info in dicts to a continuous array-account for different types of data
		for chrom in chromsToProcess:
			for bin in lookUpTable2[str(chrom)]:
				if (type(bin) is int) or (type(bin) is float) or (type(bin) is np.float64):
					A.append(bin)
				else:
					A.append(np.float64(bin[7]))
			for bin in lookUpTable0[str(chrom)]:
				if (type(bin) is int) or (type(bin) is float) or (type(bin) is np.float64):
					B.append(bin)
				else:
					B.append(np.float64(bin[7]))
			for bin in lookUpTable1[str(chrom)]:
				if (type(bin) is int) or (type(bin) is float) or (type(bin) is np.float64):
					C.append(bin)
				else:
					C.append(np.float64(bin[7]))

		print(len(A),len(chromEnd))
		#plot info with labelled horizontal lines to show chr boundaries
		plt.plot(A,marker='+',markersize=1,color='r')
		plt.plot(B,marker='s',color='b',markersize=2)
		plt.plot(C,'o',markersize=1,color='k')
		plt.title('genome_' + str(item))
		plt.annotate(str(1),xy=(1,3),fontsize=8,horizontalalignment='left')
		i=2;
		for vitm in chromEnd[:-1]:
			plt.axvline(x=vitm,color='b',hold='True')
			plt.annotate(str(i),xy=(vitm+1,3),fontsize=8,horizontalalignment='left')
			i=i+1;
#		plt.xlim(0.7,1.2)
#		plt.ylim(-15,15)
		plt.xlabel('Bins')
		#save fig for this sample
		plt.ylabel(label)
		pdf.savefig()

