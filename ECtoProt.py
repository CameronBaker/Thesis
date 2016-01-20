import sys
import urllib2
import subprocess
import commands
import time

def motifSearch(motif):
	url = 'http://www.rcsb.org/pdb/rest/search'

	queryText = """
	<orgPdbCompositeQuery version="1.0">
	<queryRefinement>
	<queryRefinementLevel>0</queryRefinementLevel>
	<orgPdbQuery>
	<version>head</version>
	<queryType>org.pdb.query.simple.EnzymeClassificationQuery</queryType>
	<description>Enzyme Classification Search : EC=%s</description>
	<Enzyme_Classification>%s</Enzyme_Classification>
	</orgPdbQuery>
	</queryRefinement>
	</orgPdbCompositeQuery>
	""" % (motif,motif)

	req = urllib2.Request(url, data=queryText)

	f = urllib2.urlopen(req)
	for i in f:
		i = i.strip()
		#p = subprocess.Popen("wget pfam.xfam.org/structure/%s"%i,shell=True)
		time.sleep(0.2)
		print i
	#p.terminate()
	#p.wait()
	#subprocess.Popen(['reset']).wait()

with open(sys.argv[1]) as f:
	lines = (line.rstrip('\n') for line in f)
	for l in lines:
		motifSearch(l)
	exit()
