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
import urllib2
import re

def ActiveSites(prot):

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
	if len(ActiveRes) != 0:
		SelectStatement = "+".join(ActiveRes)
		cmd.create('%s_active'%(prot),'resi %s in %s'%(SelectStatement,prot))
		cmd.delete(prot)
	else:
		print "no active sites found"
		
	reponse.close()

cmd.extend("ActiveSites", ActiveSites)
print "ActiveSites command added"

