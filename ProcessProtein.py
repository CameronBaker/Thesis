'''
ProcessProtein.py

This file contains the methods responsible for obtaining
the information required from a given protein.
It includes methods to discern the functional amino acids
from the protein as found in the CSA as well as details specific
to those amino acids such as what residue it is and it's relation 
to others within the functional group

Cameron Baker
'''

from pymol import cmd
from pymol import stored
import urllib2
import re

def ActiveSites(prot):
	'''
	This function implements a command into pymol that allows the user to 
	isolate active sites for a given protein from the pdb file as per locations
	denoted within the catalytic site atlas. It involves the parsing of the 
	locations from the entry for the protein within the atlas.
	
	If active sites are found, a new object is built within pymol containing
	only those residues in their proper orientation
	
	No action if no active sites are found
	
	Internet connection needed
	'''
	cmd.fetch(prot)
	ActiveRes = []
	response = urllib2.urlopen('https://www.ebi.ac.uk/thornton-srv/databases/CSA/SearchResults.php?PDBID=%s'%prot)
	
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
	print "Protein recieved"
	if len(ActiveRes) != 0:
		SelectStatement = "+".join(ActiveRes)
		cmd.create('%s_active'%(prot),'resi %s in %s'%(SelectStatement,prot))
		cmd.delete(prot)
		print "Active sites loaded"
	else:
		print "no active sites found"
		
	response.close()

cmd.extend("ActiveSites", ActiveSites)
print "ActiveSites command added"

def SiteDistances(prot):
	'''
	This function implements a command within pymol that will find the 
	centermost residue within a protein. It works by reading out all
	of the distances between a given residue's alpha carbon and each other 
	residue's alpha carbons. The residue with the smallest average distance 
	across all residues will be the centermost.
	
	Once the centermost residue is determined, the angle between that 
	residue's alpha carbon and each other residue's alpha and beta
	carbons will be determined and returned to the user.
	'''
	stored.list = []
	cmd.iterate("(name ca)","stored.list.append((resi,resn))")
	print stored.list
	
cmd.extend("SiteDistances",SiteDistances)
print "SiteDistances command added"
