from pymol import cmd
import urllib2
import re

import sys
import urllib2
import subprocess
import commands
import time

def FindAndAlign(prot1, prot2):
	'''
	This function locates the active residues for the given proteins as described
	in the catalytic site atlas and aligns them.
	
	Loading of the 2 proteins is not necessary before the running of the script
	Connection to the internet is required to download information from the CSA
	
	ARGS: the 2 proteins in question. order does not matter
	OUTPUT: Only the active sites for the given protein will be selected. Will call
			the built in align function between them.
	RETURNS: -1 for missing protein or the alignment score
	'''
	
	#For each protein
	missing = 0
	for prot in [prot1,prot2]:
	
		#Retrieve the protein and it's CSA file
		cmd.fetch(prot)
		ActiveRes = []
		response = urllib2.urlopen('https://www.ebi.ac.uk/thornton-srv/databases/CSA/SearchResults.php?PDBID=%s'%prot)
		
		#Iterates through the file to parse active site locations
		for line in response:
			#Regex matches active sites for either 1, 2, or 3 digits and appends it to the active site list
			for i in range(1,4):
				exp = "<td>"+'.'*i+"</td>"
				CurRes = re.findall(exp,line)
				for res in CurRes:
					res = res.split('>')[1].split('<')[0]
					if res.isdigit():
						ActiveRes.append(res)
		
		#Only keeps unique residues
		ActiveRes = list(set(ActiveRes))
		
		#Align against the full protein if there are no active residues
		#As long as one protein has an active site it should align well		
		if len(ActiveRes) == 0:
			missing = missing+1
			cmd.set_name(prot,'%s_active'%(prot))
			if missing == 2:
				return -1
		else:
			SelectStatement = "+".join(ActiveRes)
			cmd.create('%s_active'%(prot),'resi %s in %s'%(SelectStatement,prot))
			cmd.delete(prot)
		#Concatenates the active residues into the selection query
		
		#Creates a new object within PyMOL derived from the active site residues
		
		#Deletes the original protein and closes the html connection
		response.close()
	#Aligns the active sites as created by the previous loop and returns the RMSD
	aln = cmd.super('%s_active'%(prot1),'%s_active'%(prot2))
	return aln[0]

cmd.extend("fa",FindAndAlign)

def ComputeOverlap(EC,pfam):
	'''
	This function finds the list of proteins matching a given motif and generates a score matrix by
	matching each of those proteins using the FindAndAlign function, as shown above
	
	The function queries the pdb with the given EC to retrieve the desired proteins
	'''
	
	url = 'http://www.rcsb.org/pdb/rest/search'
	protList = []

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

	f = urllib2.urlopen(req)
	for i in f:
		i = i.strip()
		protList.append(i)

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
	queryList = []
	
	f = urllib2.urlopen(req)
	for i in f:
		i = i.strip()
		queryList.append(i)

	#Obtain the overlap between the proteins for a given EC and pfam ID
	protList = set(protList).intersection(queryList)
	protList = list(protList)
	ScoreMatrix = [[-2 for x in range(len(protList))] for x in range(len(protList))]
	numComplete = 0

        #For each protein, find the distance based overlap between that protein and all others
	for p1 in range(0,len(protList)):
		for p2 in range(p1,len(protList)):
			if protList[p1] == protList[p2]:
				score = 0
			else:
				score = FindAndAlign(protList[p1],protList[p2])
				#cmd.delete('all')
			ScoreMatrix[p1][p2] = score
			ScoreMatrix[p2][p1] = score
			#print protList[p1] + " " + protList[p2] + " "+ str(ScoreMatrix[p1][p2])
		numComplete = numComplete + 1
		print "%d/%d proteins completed"%(numComplete,len(protList))
	
	#Prints the results to a csv file for further use
		
	f = open('%s.(%s).csv'%(pfam,EC),'w')
	header = 'Protein,'
	header = header + ",".join(protList)
	header = header[:-1] + "\n"
	f.write(header)
	
	for i in range(0,len(ScoreMatrix)):
		line = protList[i]
		for num in ScoreMatrix[i]:
			line = line+","+str(num)
		line = line+"\n"
		f.write(line)
	f.close()
	print "Computation finished"
 
cmd.extend("ComputeOverlap",ComputeOverlap)




	
