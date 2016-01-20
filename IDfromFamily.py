import sys
import urllib2
  
def pfamSearch(pfam):
	url = 'http://www.rcsb.org/pdb/rest/search'

	queryText = """
	<orgPdbCompositeQuery version="1.0">
	<queryRefinement>
	<queryRefinementLevel>0</queryRefinementLevel>
	<orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.PfamIdQuery</queryType>
    <description>Pfam Accession Number:  %s</description>
    <pfamID>%s</pfamID>
  	</orgPdbQuery>
	</queryRefinement>
	</orgPdbCompositeQuery>
	""" % (pfam,pfam)

	req = urllib2.Request(url, data=queryText)

	f = urllib2.urlopen(req)
	for i in f:
		i = i.strip()
		print i

pfamSearch(sys.argv[1])