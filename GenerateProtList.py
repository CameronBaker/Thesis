'''
GenerateProtList.py

This file contains the methods specific for generating a list of proteins derived from the 
overlap between Enzyme classification class and pfam family 
It is the final iteration to past supplimentary scripts ECtoProt.py and IDfromFamily.py

To be used across the remainder of the project
Internet connection is required for this script to work

Cameron Baker
'''

import urllib2

def GenerateList(EC,pfam):
	'''
	Given an EC and pfam ID, query the PDB and return the overlap between both datasets
	for a list of proteins to be used later in the pipeline
	'''
	
	url = 'http://www.rcsb.org/pdb/rest/search'
	ECList = []

    #Query the PDB to obtain a list of proteins for a given EC
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
	""" % (EC,EC)

	req = urllib2.Request(url, data=queryText)
	queryResult = urllib2.urlopen(req)
	
	for line in queryResult:
		protein = line.strip()
		ECList.append(protein)

    #Query the PDB to obtain a list of proteins for a given pfam ID
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
	PfamList = []
	
	queryResult = urllib2.urlopen(req)
	
	for line in queryResult:
		protein = line.strip()
		PfamList.append(protein)

	#Obtain the overlap between the proteins for a given EC and pfam ID
	protList = set(PfamList).intersection(ECList)
	protList = list(protList)
	
	return protList
	