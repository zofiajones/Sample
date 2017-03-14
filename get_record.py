##FUNCTION: Take gene list and query NCBI to get bed file
##INPUT: stdin from get_gene.sh
##OUTPUT: either written to bed file or added to exceptions list

import re
import sys
import urllib.request
import xml.etree.ElementTree as ET

gene=sys.argv[1]
m=[]
for item in sys.argv[2:]:
	m.append(item)

found=0
file1=open('gene.bed','a')
file2=open('not_found.txt','a')

for id1 in m:
	#get NCBI entry
	request=urllib.request.Request('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=' + id1 + '&retmode=xml')
	response=urllib.request.urlopen(request)
	tree = ET.parse(response)
	root = tree.getroot()
	found=0
	gene_ref=''
	#check if provided name matches reference name
	if gene.upper() in root[0][3][0][0].text.upper() or root[0][3][0][0].text.upper() in gene.upper():
		gene_ref=gene;
	else:
	#or else check if provided name matches alternate name
		for item in root[0][3][0].iter(tag='Gene-ref_syn'):
			for alias in item.iter(tag='Gene-ref_syn_E'):
				if gene.upper() in alias.text.upper() or alias.text.upper() in gene.upper():
					gene_ref=gene;
	#find first ie latest record for genomic location for GRCh37
	if gene_ref:
		for item in root[0].iter(tag='Entrezgene_comments'):
			for item0 in item.iter(tag='Gene-commentary'):
				if item0[0].text == '254' and item0[1].text == 'Gene Location History':
					for item2 in item0[2].iter(tag='Gene-commentary_comment'):
						if 'GRCh37' in item2[0][1].text  and item2[0][2].text == 'Reference':
							for item3 in item2[0].iter(tag='Gene-commentary_comment'):
								if item3[0][0].attrib['value'] =='genomic':
									m0=re.search('[\d+XYM]',item3[0][1].text)
									chr=m0.group(0)
									start=str(int(item3[0][4][0][0][0][0].text)-1)
									end=item3[0][4][0][0][0][1].text
									file1.write(chr.strip() + "\t" + start.strip() + "\t" + end.strip() + "\t" + 'GENE=' + root[0][3][0][0].text.strip() + ';NCBI_id=' + id1 + '\n')
									found=1
							break
	else:
		file1.write(gene + '\n')

file1.close()
file2.close()
