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
	
	#In chain is used to work around the respresentation of residues in CSA
	#The CSA hosts both the active site within the protein and the same active
	#sites in homologous proteins. If these did not match, it will incorrectly select
	#the sites for both the actual protein and the homologous as determined by CSA
	inchain = True
	
	for line in response:
			#Regex matches active sites for either 1, 2, or 3 digits and appends it to the active site list
		for i in range(1,4):
			exp = "<td>"+'.'*i+"</td>"
			CurRes = re.findall(exp,line)
			for res in CurRes:
				res = res.split('>')[1].split('<')[0]
				if res.isdigit() and inchain:
					ActiveRes.append(res)
					inchain = False
				else:
					inchain = True
				
		
		#Only keeps unique residues
		ActiveRes = list(set(ActiveRes))
		
		#Align against the full protein if there are no active residues
		#As long as one protein has an active site it should align well	
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
	
	#Retrieve the list of residues
	stored.list = []
	cmd.iterate("(%s and name ca)"%(prot),"stored.list.append((resi,resn))")
	#print stored.list
	
	#Initialization of variables for finding the centermost residue
	res_num = len(stored.list)
	closest_res_res = 0
	smallest_dist = 100
	for res in stored.list:
		temp_dist = 0
		for res2 in stored.list:
			if res[0] != res2[0]:
				temp_dist = temp_dist + \
							cmd.get_distance("resi %s in %s and name ca"%(res[0],prot),
								         "resi %s in %s and name ca"%(res2[0],prot))
		temp_dist = temp_dist / (res_num - 1)
		if temp_dist < smallest_dist:
			smallest_dist = temp_dist
			closest_res = res
	#print closest_res
	#print smallest_dist
	
	#like protein + profile, get it?
	PROfile = {}
	
	#The following calculates the distance between the anchor residue and others
	#It then sorts the ordering and then calculates angles
	print "anchor_residue,anchor_num,target_residue,target_num,distance,angle"
	for res in stored.list:
		if res != closest_res:
			dist = cmd.get_distance("resi %s in %s and name ca"%(closest_res[0],prot),
							         "resi %s in %s and name ca"%(res[0],prot))
			PROfile[res[0]] = dist
	PROfile = sorted(PROfile.items(), key=lambda x: x[1])
	for res in stored.list:
		for site in range(0,len(PROfile)):
			angle = cmd.get_angle("resi %s in %s and name ca"%(closest_res[0],prot),
							         "resi %s in %s and name ca"%(res[0],prot),
									 "resi %s in %s and name cb"%(res[0],prot))
			if res[0] == PROfile[site][0]:
				PROfile[site] = res,PROfile[site][1],angle
	for res in PROfile:
		print "%s,%s,%s,%s,%s,%s"%(closest_res[1],closest_res[0],res[0][1],res[0][0],res[1],res[2])
	
	
	#Uncomment the following 4 lines for visual representation of distances between centermost residue
	#for res in stored.list:
	#	if res != closest_res:
	#		cmd.distance("resi %s in %s and name ca"%(closest_res[0],prot),
	#							         "resi %s in %s and name ca"%(res[0],prot))
	
cmd.extend("SiteDistances",SiteDistances)
print "SiteDistances command added"
