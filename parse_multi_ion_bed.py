##FUNCTION: parse tmap output to get basechange as well as alt freq before and after filtering 
##INPUT: stdin from bedtools output
##OUTPUT: stdout to required format to match with expected variants (from either dPCR or COSMIC)
import re
import os
import sys
import pysam

#open genome fasta
file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')

for line in sys.stdin:
	#get raw info
	m=line.split('\t')
	alt=m[5].split(',')
	AF=re.search('AF=([\d\.,]+)',m[6]).group(1).split(',')
	SAF=re.search('SAF=([,\d]+)',m[6]).group(1).split(',')
	SAR=re.search('SAR=([,\d]+)',m[6]).group(1).split(',')
	SRF=re.search('SRF=([,\d]+)',m[6]).group(1)
	SRR=re.search('SRR=([,\d]+)',m[6]).group(1)
	tags=m[7].split(':')
	#check for lines with no info
	if len(m[8].split(':'))>1:
		ind=tags.index('AO')
		AO=m[8].split(':')[ind].split(',')
		ind1=tags.index('RO')
		RO=m[8].split(':')[ind1]
		freq0=[]
		total_reads0=float(SRR)+float(SRF)
		#get coverage and frequency of unfiltered reads for each alt (SAF,SAR,SRF,SRR)
		for i in range(len(alt)):
			total_reads0=float(SAR[i])+float(SAF[i])+total_reads0
		for i in range(len(alt)):
			if total_reads0 > 0:
				freq0.append(str( ( float(SAF[i])+float(SAR[i]) ) / total_reads0  ))
			else:
				freq0.append('0')
		#get coverage and frequency of filtered reads for each alt (AO and RO)
		freq=[]
		total_reads=float(RO)
		for i in range(len(alt)):
			total_reads=total_reads+float(AO[i])
		for i in range(len(alt)):
			if total_reads > 0:
				freq.append( float(AO[i])/ total_reads  )
			else:
				freq.append(0.0)
		#get alt with highest freq of filtered reads
		j=freq.index(max(freq))
		#put together info colum: base change, freq(unfiltered),total(unfiltered),freq(filtered),freq(unfiltered)
		info = m2[j] + ',' + str(freq[j]) + ',' + str(total_reads) + ',' + str(freq0[j]) + ',' + str(total_reads0)
		#change . to O for easier regex
		alt1=alt[j].replace('.','O')
		#cut leading or trailing bases that match reference genome (helpful particularly for freebayes based variant callers)
		mref=len(m[4])
		malt=len(alt1)
		trim=0
		#get max no of iterations
		if malt < mref:
			mx=range(malt);
		elif mref < malt:
			mx=range(mref);
			for i in mx:
				if m[4][mref-i] == alt1[malt-i]:
					loc=m[0]+':'+ str(int(m[2])+mref-1)+'-'+str(int(m[2])+mref-1)
					base=file.fetch(region=loc)
					if base.strip() in m[4][mref-i-1]:
						trim=trim+1
					else:
						break
				else:
					break
		#get final reference and alt sequence
		ref=m[4][0:mref-trim]
		alt=alt1[0:malt-trim]
		#make sure ref and alt are the same length for exact checking
		while len(alt) < len(ref):
			alt=alt+'O'
		while len(ref)<len(alt):
			ref=ref+'O'
		for i in range(len(ref)):
			#don't print any anchor bases
			if alt[i] != ref[i]:
				print(('\t').join([m[0],str(int(m[1])+i),str(int(m[2])+i),m[3],ref[i],alt[i],info.strip()]))
