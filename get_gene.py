##FUNCTION: to search NCBI for a list of potential gene IDs
##INPUT: gene name (command line input)
##OUTPUT: stdout of search result

import re
import sys
import urllib.request

gene=str(sys.argv[1])
print(gene)
request=urllib.request.Request('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=human[organism]+AND+' + gene + '[Gene Name]' + '+AND+"has%20ccds"[Properties]')
response=urllib.request.urlopen(request)
html=response.read()
print(html)
