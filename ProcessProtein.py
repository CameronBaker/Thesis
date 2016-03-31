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

'''
Interesting notes
Duplicate features in 1g66, 1h1w
'''

from pymol import cmd
from pymol import stored
import urllib2
import re
import os

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
		cmd.remove('%s and not chain a'%(prot))
		SelectStatement = "+".join(ActiveRes)
		cmd.create('%s_active'%(prot),'resi %s in %s'%(SelectStatement,prot))
		cmd.delete(prot)
		#print "Active sites loaded"
		return 1
	else:
		cmd.delete(prot)
		print "no active sites found for "+prot
		return 0
	response.close()

cmd.extend("ActiveSites", ActiveSites)
print "ActiveSites pdbid"

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
			try:
				if res[0] != res2[0]:
					temp_dist = temp_dist + \
								cmd.get_distance("resi %s in %s and name ca"%(res[0],prot),
											"resi %s in %s and name ca"%(res2[0],prot))
			except:
				stored_tmp = stored.list
				if res[0] != res2[0]:
					stored.list = []
					cmd.iterate("(resi %s in %s and name ca)"%(res[0],prot),"stored.list.append((ID))")
					if len(stored.list) > 1:
						cmd.remove("id %s"%(stored.list[-1]))
					stored.list = []
					cmd.iterate("(resi %s in %s and name ca)"%(res2[0],prot),"stored.list.append((ID))")
					if len(stored.list) > 1:
						cmd.remove("id %s"%(stored.list[-1]))
					temp_dist = temp_dist + \
								cmd.get_distance("resi %s in %s and name ca"%(res[0],prot),
											"resi %s in %s and name ca"%(res2[0],prot))
				stored.list = stored_tmp
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
	for res in stored.list:
		#if res != closest_res:
		dist = cmd.get_distance("resi %s in %s and name ca"%(closest_res[0],prot),
							         "resi %s in %s and name ca"%(res[0],prot))
		PROfile[res[0]] = dist
	PROfile = sorted(PROfile.items(), key=lambda x: x[1])
	for res in stored.list:
		for site in range(0,len(PROfile)):
			try:
				if res[1] == "GLY":
					angle = cmd.get_angle("resi %s in %s and name ca"%(closest_res[0],prot),
						         "resi %s in %s and name ca"%(res[0],prot),
								 "resi %s in %s and name n"%(res[0],prot))
				else:
					angle = cmd.get_angle("resi %s in %s and name ca"%(closest_res[0],prot),
						         "resi %s in %s and name ca"%(res[0],prot),
								 "resi %s in %s and name cb"%(res[0],prot))
			except:
				print "Duplicate atoms found for a given carbon beta. Deleting duplicate"
				stored.list = []
				cmd.iterate("(resi %s in %s and name ca)"%(res[0],prot),"stored.list.append((ID))")
				if len(stored.list) > 1:
					cmd.remove("id %s"%(stored.list[-1]))
				cmd.iterate("(resi %s in %s and name cb)"%(res[0],prot),"stored.list.append((ID))")
				if len(stored.list) > 1:
					cmd.remove("id %s"%(stored.list[-1]))
				if res[1] == "GLY":
					angle = cmd.get_angle("resi %s in %s and name ca"%(closest_res[0],prot),
						         "resi %s in %s and name ca"%(res[0],prot),
								 "resi %s in %s and name n"%(res[0],prot))
				else:
					angle = cmd.get_angle("resi %s in %s and name ca"%(closest_res[0],prot),
						         "resi %s in %s and name ca"%(res[0],prot),
								 "resi %s in %s and name cb"%(res[0],prot))
								 
			if res[0] == PROfile[site][0]:
				PROfile[site] = res,PROfile[site][1],angle
	
	#Uncomment the following for numbers related to the residues in the structure
	#print "anchor_residue,anchor_num,target_residue,target_num,distance,angle"
	#for res in PROfile:
	#	print "%s,%s,%s,%s,%s,%s"%(closest_res[1],closest_res[0],res[0][1],res[0][0],res[1],res[2])
	
	#Uncomment the following 4 lines for visual representation of distances between centermost residue
	#for res in stored.list:
	#	if res != closest_res:
	#		cmd.distance("resi %s in %s and name ca"%(closest_res[0],prot),
	#							         "resi %s in %s and name ca"%(res[0],prot))
	
	#Comment out the delete command to keep residues after quantification. Good for visualization
	cmd.delete("all")
	return PROfile
	
cmd.extend("SiteDistances",SiteDistances)

def UPfromProt(pdb):
	'''
	This function queries the uniprot database for a specific protein 
	and returns the uniprot families that it belongs to
	'''
	
	try:
		response = urllib2.urlopen('https://www.ebi.ac.uk/thornton-srv/databases/CSA/SearchResults.php?PDBID=%s'%(pdb),timeout=20)
	except:
		response = urllib2.urlopen('https://www.ebi.ac.uk/thornton-srv/databases/CSA/SearchResults.php?PDBID=%s'%(pdb),timeout=20)
	UPList = []
	
	#Parse the Uniprot family out of the webpage
	for line in response:
		exp = "UNIID=[A-Z][0-9]+"
		curLine = re.findall(exp,line)
		if len(curLine) > 0:
			UPList.append(curLine[0].split("=")[1])
	
	if len(UPList) is 0:
		return 0
	return UPList

cmd.extend("Prot2UP",UPfromProt)

def PFfromProt(pdb):
	'''
	The search for pfam ID's requires querying uniprot with a uniprot ID and taking the associated pfam ID's
	for the given uniprot
	'''
	
	UPList = UPfromProt(pdb)
	
	if UPList == 0:
		return 0
	
	pfamID = []
	for id in UPList:
		try:
			response = urllib2.urlopen('http://www.uniprot.org/uniprot/%s'%(id),timeout=20)
		except:
			response = urllib2.urlopen('http://www.uniprot.org/uniprot/%s'%(id),timeout=20)
		
		for line in response:
			exp = "PF\\d+"
			curLine = re.findall(exp,line)
			if len(curLine) > 0:
				for pfam in curLine:
					if pfam not in pfamID:
						pfamID.append(pfam)
	
	return pfamID
	
cmd.extend("Prot2pfam",PFfromProt)

def ECtoProt(EC):
	'''
	This function retrieves the list of proteins for a given EC from the protein data bank
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

	#Parse the query into a usable list
	for line in queryResult:
		protein = line.strip()
		ECList.append(protein.split(":")[0])
	
	return ECList

def adjustRes(anchor, res, distance, angle, protname):
	
	stored.list = []
	
	if res[0] == anchor[0][0]:
		return
		
	cmd.iterate("(resi %s in %s and name ca)"%(res[0],protname),"stored.list.append((ID))")
	if len(stored.list) > 1:
		cmd.remove("id %s"%(stored.list[0]))
	stored.list = []
	cmd.iterate("(resi %s in %s and name ca)"%(anchor[0][0],protname),"stored.list.append((ID))")

	if len(stored.list) > 1:
		cmd.remove("id %s"%(stored.list[0]))
	dist = cmd.get_distance("resi %s in %s and name ca"%(anchor[0][0],protname),
							    "resi %s in %s and name ca"%(res[0],protname))
								
	stored.list = []
	cmd.iterate("(resi %s in %s and name cb)"%(res[0],protname),"stored.list.append((ID))")
	if len(stored.list) > 1:
		cmd.remove("id %s"%(stored.list[0]))
	
	#print angle

	diff = distance - dist
	
	anchor_ca_coord = cmd.get_coords("resi %s in %s and name ca"%(anchor[0][0],protname))
	other_ca_coord = cmd.get_coords("resi %s in %s and name ca"%(res[0],protname))
	
	vx = anchor_ca_coord[0][0] - other_ca_coord[0][0]
	vy = anchor_ca_coord[0][1] - other_ca_coord[0][1]
	vz = anchor_ca_coord[0][2] - other_ca_coord[0][2]
	res_vector = [vx,vy,vz]
	translation_vector = [res_vector[0] * diff,res_vector[1] * diff,res_vector[2] * diff]
	
	cmd.translate(translation_vector,"resi %s in %s"%(res[0],protname))

	if res[1] == "GLY":
		cmd.set_dihedral("resi %s in %s and name ca"%(anchor[0][0],protname),
				 "resi %s in %s and name ca"%(res[0],protname),
				 "resi %s in %s and name ca"%(res[0],protname),
				 "resi %s in %s and name n"%(res[0],protname),angle)
	else:
		cmd.set_dihedral("resi %s in %s and name ca"%(anchor[0][0],protname),
				 "resi %s in %s and name ca"%(res[0],protname),
				 "resi %s in %s and name ca"%(res[0],protname),
				 "resi %s in %s and name cb"%(res[0],protname),angle)

def adjustWrapper(prot):

	ActiveSites(prot)
	profile = SiteDistances("%s_active"%(prot))
	for line in profile:
		print line
	ActiveSites(prot)
	anchor = []
	for res in profile:
		if res[1] == 0:
			anchor = res
	for res in profile:
		if res is not anchor:
			adjustRes(anchor,res[0],res[1],res[2])

cmd.extend("translate",adjustWrapper)
	
def generateMotif(EC,entry_minimum,protein_cap,db,outdir):
	
	cmd.cd(outdir)
	proteinList = ECtoProt(EC)
	UPDictionary = {}
	using_windows = False
	if os.name == "nt":
		using_windows = True
	#cmd.feedback('disable','all','actions')
	#cmd.feedback("disable","all","results")
	
	#Information for the user based on the initial retrieval of proteins
	print "%s proteins found for %s"%(len(proteinList),EC)
	print "Beginning protein grouping"
	count = 0
	
	#Used to contain information on incomplete records (no link between catalytic site atlas and UniProt)
	NullCount = 0
	NullList = []
	
	#Subset the protein list for now
	if protein_cap != "none":
		if len(proteinList) > protein_cap:
			proteinList = proteinList[0:protein_cap]
	
	
	for protein in proteinList:
	
		if db == "pfam":
			UPID = PFfromProt(protein)
		else:
			UPID = UPfromProt(protein)
			
		if UPID is 0:
			count = count + 1
			if count % 5 == 0:
				print "%s/%s"%(count,len(proteinList))
			NullCount = NullCount + 1
			NullList.append(protein)
			continue
		if len(UPID) > 1:
			for ID in UPID:
				if ID not in UPDictionary.keys():
					UPDictionary[ID] = protein
				else:
					UPDictionary[ID] = UPDictionary[ID]+","+protein
		else:
			if UPID[0] not in UPDictionary.keys():
				UPDictionary[UPID[0]] = protein
			else:
				UPDictionary[UPID[0]] = UPDictionary[UPID[0]]+","+protein
		count = count + 1
		if count % 5 == 0:
			print "%s/%s"%(count,len(proteinList))
	
	print "%s proteins with incomplete records"%(NullCount)
	
	for entry in UPDictionary.keys():
		print "%s: %s"%(entry,UPDictionary[entry])
		ProteinSet = UPDictionary[entry].split(",")
		ProteinSet = list(set(ProteinSet))
		ProfileSet = []
		DistDict = {}
		AngleDict = {}
		count = 0
		
		if not int(entry_minimum) < len(ProteinSet):
			print "Not enough proteins for entry %s"%entry
			continue
		
		for protein in ProteinSet:
			#This generates the active sites and quanitifies the distances and angles
			#Only keeps the active site if it successfully generated
			#Successful generation is due to the active sites prescense within the CSA
			count = count + 1
			if ActiveSites(protein) == 1:
				protein = protein + "_active"
				PROfile = SiteDistances(protein)
				ProfileSet.append(PROfile)
			if count % 5 == 0:
				print "%s/%s"%(count,len(ProteinSet))
				
				#Uncomment the following to print out the distances and angles for a given protein
				#for res in PROfile:
				#	print res
		print "%s/%s"%(len(ProteinSet),len(ProteinSet))
		
		if len(ProfileSet) > 0:
		
			#Below is finding the average distance and angle between each residue
			#within the active site
			consensusDic = []
			anchor = []
			ResCount = {}
			for PROfile in ProfileSet:
				for res in PROfile:
					key = res[0]
					if res[1] == 0:
							anchor = res
							
					if key not in DistDict.keys():
						DistDict[key] = res[1]
					else:
						DistDict[key] = (DistDict[key] + res[1]) / 2
						
					if key not in AngleDict.keys():
						AngleDict[key] = res[2]
					else:
						AngleDict[key] = (AngleDict[key] + res[2]) / 2
					
					if key not in ResCount.keys():
						ResCount[key] = 1
					else:
						ResCount[key] = ResCount[key] + 1
			print "------------------------------"
			print "MOTIF STATS for %s_%s:"%(EC,entry)
			print "Different stats: %s across %s proteins"%(str(len(DistDict.keys())),len(ProteinSet))
			print "Anchor residue will have a Distance of 0 and an angle of 90"
			print "Residue Position, Residue, Distance from anchor, Angle, Count"
			
			for key in DistDict:
				consensusDic.append([[key[0],key[1]],DistDict[key],AngleDict[key],ResCount[key]])		
			
			for line in consensusDic:
				print "%s,%s,%s,%s,%s"%(line[0][0],line[0][1],line[1],line[2],line[3])
			
			model = UPDictionary[entry].split(",")[0]
			count = 0
			while True:
				if ActiveSites(model) == 1:
					break
				count = count + 1
				if count >= len(UPDictionary[entry].split(",")):
					print "No active sites found"
					break
				model = UPDictionary[entry].split(",")[count]
			
			if count == len(UPDictionary[entry].split(",")):
				cmd.remove("all")
				continue
			
			modelInfo = SiteDistances("%s_active"%model)
			modelSites = []
			ActiveSites(model)
			
			for site in modelInfo:
				tmp = [site[0][0],site[0][1]]
				modelSites.append(tmp)
			
			ActiveSite = []
			SubstituteSet = []
			ResMax = {}
			ResMap = {}
			for key in consensusDic:
				if key[0][0] not in ResMax:
					ResMax[key[0][0]] = key[3]
					ResMap[key[0][0]] = key[0][1]
				else:
					if ResMax[key[0][0]] < key[3]:
						ResMax[key[0][0]] = key[3]
						ResMap[key[0][0]] = key[0][1]
			
			newCon = []
			for key in consensusDic:
				if key[0][1] == ResMap[key[0][0]]:
					newCon.append(key)
				else:
					if key[0][0] not in SubstituteSet:
						SubstituteSet.append(key[0][0])
						
			consensusDic = newCon
			
			for site in SubstituteSet:
				cmd.remove("resi %s and not name ca+cb"%(site))
			
			#Uncomment this save and save in line 552 for pre and post translational structures.
			#comment out 553
			#cmd.save("%s_active_pre.pdb"%model,"%s_active"%model)
			
			for key in consensusDic:
				if key[0] in modelSites:
					adjustRes(anchor, key[0], key[1], key[2], model)
				
			#cmd.save("%s_active_post.pdb"%model,"%s_active"%model)
			cmd.save("%s_%s.pdb"%(EC,entry),"%s_active"%model)
			cmd.remove("all")
			
			for prot in UPDictionary[entry].split(","):
				if using_windows:
					os.system("del %s.pdb"%prot)
				else:
					os.system("rm %s.pdb"%prot)
			
			print ""
			print "%s_%s written"%(EC,entry)
			print ""
		else:
			print "skipping %s. No usable proteins found"%(entry)
			print ""
			UPDictionary.pop(entry)
	
	
cmd.extend("GenerateMotif",generateMotif)
print "GenerateMotif EC, MinimumProteinsForGroup (number), MaxmimumProteinsTotal (none or a number), pfam/uniprot, OutputDirectory"
