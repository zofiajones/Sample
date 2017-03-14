##FUNCTION: get available mutations and create vcf file
##INPUT: list of gene names
##OUTPUT:
## if wildtype:
##	if genomic_location:
##		write to vcf file
##	if not genomic_location:
##		write gene and aa change to send to COSMIC script
##if not wildtype:
##	write cell line description and engineered mutation

import subprocess
import os
import re
import MySQLdb

#Initialise mysql
db=MySQLdb.connect(db="mysql",read_default_file="~/config.cnf")
c=db.cursor()

#open files
file=open('/home/zjones/Panel_beds/genes_1.txt','r')
lines=file.readlines()
file1=open('output_gene.txt','w')
file2=open('output_gene.vcf','w')
file3=open('aa_COSMIC.txt','w')
file4=open('exceptions.txt','w')

for line in lines:
	print('line',line)
	protein=''
	genomic=''
	# get mutation ids
	cmd0='select id,genomic_change,protein_change from mutations where gene like \'' +  line.strip()  +  '\';'
	c.execute(cmd0)
	results=c.fetchall()
	if len(results)>0:
		for entry0 in results:
			genomic=entry0[1]
			protein=entry0[2]
			id0=entry0[0]
			des_wt=''
			des_eng=''
			par=''
			wt=0
			# get cell line id and description
			cmd0='select description , id  from cell_lines where id in (select cell_line_id from mutation_cell_line where mutation_id=\'' + str(id0)   + '\');'
			c.execute(cmd0)
			results=c.fetchall()
			cell_id=[]
			description=[]
			for entry1 in results:
				if 'wildtype' in entry1[0]:
					wt=1;
					des_wt=entry1[0] + ';' + des_wt
					par= entry1[0].split()[0] + ';' + par
				else:
					des_eng=entry1[0] + ';' + des_eng
					cell_id.append(entry1[1])
					description.append(entry1[0])
			if not genomic:
				genomic='NA'
			if not protein:
				protein='NA'
			#get mutation loc
			if 'NA' not in genomic:
				m=genomic.replace('g.','').split(':')
				chr=m[0]
				m0=re.search('(\d+)',m[1])
				loc=m0.group(1)
				if '>' in m[1]:
					n=m[1].split('>')
					m0=re.search('([A-Z]+)',n[0])
					ref=m0.group(1)
					alt=n[1]
				elif 'del' in m[1]:
					n=m[1].replace('\d+','').split('del')
					ref=n[1]
					alt='.'
				elif 'ins' in m[1]:
					n0=m[1].split('_')
					loc=n0[0]
					n=n0[1].replace('\d+','').split('ins')
					ref='.'
					alt=n[1]
				else:
					file4.write(m[1])
				#write to vcf
				file2.write(('\t').join([chr,loc,'.',ref,alt,'.','PASS',protein + ';' + line.strip() + ';' + par ]) + '\n')
			if wt == 0:
				cat_no=''
				#if not wildtype, write to list of engineered cell lines
				for id in cell_id:
					cmd0='select catalog_number from cell_line_catalog_numbers where cell_line_id=\'' + str(id) + '\';'
					c.execute(cmd0)
					results=c.fetchall()
					for item in results:
						cat_no=str(item[0]) + ';' + cat_no
				file1.write(line.strip() + '\t' + genomic + '\t' + protein + '\t' + cat_no   +  '\t' + des_eng + '\n')
				file3.write(line.strip() + '\t' + protein + '\n')




file.close()
file1.close()
file2.close()
file3.close()
