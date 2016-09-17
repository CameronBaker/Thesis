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
Duplicate features in 1g66
'''

from pymol import cmd
from pymol import stored
import urllib2
import re
import os
from requests import exceptions

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
	Sites = {}
	response = urllib2.urlopen('https://www.ebi.ac.uk/thornton-srv/databases/CSA/SearchResults.php?PDBID=%s'%prot)
	
	#In chain is used to work around the respresentation of residues in CSA
	#The CSA hosts both the active site within the protein and the same active
	#sites in homologous proteins. If these did not match, it will incorrectly select
	#the sites for both the actual protein and the homologous as determined by CSA
	inchain = True
	
	ActiveSite = []
	i = 1
	for line in response:
			#Regex matches active sites for either 1, 2, or 3 digits and appends it to the active site list
		exp = "<td>[0-9]+</td>"
		CurRes = re.findall(exp,line)
		if CurRes > 0:
			for s in line.split("<th>UniProtKB Number</th>"):
				cursite = "SITE "+str(i)
				#print s
				
				ActiveSite = []
				for res in CurRes:
					res = res.split('>')[1].split('<')[0]
					if res.isdigit() and inchain:
						ActiveSite.append(res)
						inchain = False
					else:
						inchain = True
				ActiveSite = list(set(ActiveSite))
				if ActiveSite not in Sites.values():
					Sites[cursite] = ActiveSite
					i = i + 1
				
		
		#Only keeps unique residues
		
		
		#Align against the full protein if there are no active residues
		#As long as one protein has an active site it should align well
		
	statements = []
	response.close()
	for value in Sites.values():
		if len(value) != 0:
			SelectStatement = "+".join(value)
			statements.append(SelectStatement)
	if len(statements) > 0:
		return statements
	else:
		print "no active sites found for "+prot
		return 0

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
	print "here: "+prot
	#Retrieve the list of residues
	stored.list = []
	cmd.iterate("(%s and name ca)"%(prot),"stored.list.append((resi,resn))")
	#print stored.list
	print prot
	#Initialization of variables for finding the centermost residue
	res_num = len(stored.list)
	closest_res = 0
	smallest_dist = 100
	if res_num == 1:
		closest_res = stored.list[0]
	else:
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
		try:
			response = urllib2.urlopen('https://www.ebi.ac.uk/thornton-srv/databases/CSA/SearchResults.php?PDBID=%s'%(pdb),timeout=20)
		except:
			print "Timeout"
			return 0
	UPList = []
	
	#Parse the Uniprot family out of the webpage
	try:
		for line in response:
			exp = "UNIID=[A-Z][0-9]+"
			curLine = re.findall(exp,line)
			if len(curLine) > 0:
				UPList.append(curLine[0].split("=")[1])
	
		if len(UPList) is 0:
			return 0
	except:
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
			try: 
				response = urllib2.urlopen('http://www.uniprot.org/uniprot/%s'%(id),timeout=20)
			except:
				return 0
		
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
	'''
	This method is responsible for translating the distance and angle
	between the anchor and another amino acid into the correct distance
	'''
	
	stored.list = []
	if res[0] == anchor[0][0]:
		return
	
	#A good chunk of the following is responsible for filtering out duplicate atoms from pdb structures
	cmd.iterate("(resi %s in %s and name ca)"%(res[0],protname),"stored.list.append((ID))")
	if len(stored.list) > 1:
		cmd.remove("id %s"%(stored.list[0]))
	stored.list = []
	cmd.iterate("(resi %s in %s and name ca)"%(anchor[0][0],protname),"stored.list.append((ID))")

	if len(stored.list) > 1:
		cmd.remove("id %s"%(stored.list[0]))
		
	#Retrieves the distance between the anchor and other residues
	dist = cmd.get_distance("resi %s in %s and name ca"%(anchor[0][0],protname),
							    "resi %s in %s and name ca"%(res[0],protname))
	#More residue cleaning
	stored.list = []
	cmd.iterate("(resi %s in %s and name cb)"%(res[0],protname),"stored.list.append((ID))")
	if len(stored.list) > 1:
		cmd.remove("id %s"%(stored.list[0]))

	#finds the difference between the average distance and the distance found between the residues
	diff = distance - dist
	
	#Retrieves the 3D coordiantes of the alpha carbons in both residues
	anchor_ca_coord = cmd.get_coords("resi %s in %s and name ca"%(anchor[0][0],protname))
	other_ca_coord = cmd.get_coords("resi %s in %s and name ca"%(res[0],protname))
	
	#Finds the translation vector
	vx = anchor_ca_coord[0][0] - other_ca_coord[0][0]
	vy = anchor_ca_coord[0][1] - other_ca_coord[0][1]
	vz = anchor_ca_coord[0][2] - other_ca_coord[0][2]
	res_vector = [vx,vy,vz]
	translation_vector = [res_vector[0] * diff,res_vector[1] * diff,res_vector[2] * diff]
	
	#Apply that vector to each atom in the other residue
	cmd.translate(translation_vector,"resi %s in %s"%(res[0],protname))

	#Sets the dihedral angle to be the consensus
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
	
def generateMotif(EC,entry_minimum,protein_cap,db,dump,outdir):
	
	cmd.cd(outdir)
	proteinList = ECtoProt(EC)
	UPDictionary = {}
	using_windows = False
	if os.name == "nt":
		using_windows = True
		
	if dump == "true":
		outfile = open('output.txt','w')
	
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
		if len(proteinList) > int(protein_cap):
			proteinList = proteinList[0:int(protein_cap)]
	
	#The following loop populates the proteins into their respective groups
	#Groups are determined by pfam or uniprot ID, as shown on a couple lines down
	for protein in proteinList:
	
		if db == "pfam":
			UPID = PFfromProt(protein)
		else:
			UPID = UPfromProt(protein)
			
		print UPID
		if UPID is 0:
			count = count + 1
			if count % 5 == 0:
				print "%s/%s"%(count,len(proteinList))
			NullCount = NullCount + 1
			NullList.append(protein)
			continue
		if len(UPID) > 0:
			for ID in UPID:
				if ID not in UPDictionary.keys():
					UPDictionary[ID] = protein
				else:
					UPDictionary[ID] = UPDictionary[ID]+","+protein
		else:
			NullCount = NullCount + 1
			NullList.append(protein)
		count = count + 1
		if count % 5 == 0:
			print "%s/%s"%(count,len(proteinList))
	
	print "%s proteins with incomplete records"%(NullCount)
	
	
	for entry in UPDictionary.keys():
		print "%s: %s"%(entry,UPDictionary[entry])
		if dump == "true":
			outfile.write("%s: %s"%(entry,UPDictionary[entry]))
			outfile.write('\n')
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
			selects = ActiveSites(protein)
			if selects != 0:
				for select in selects:
					cmd.fetch(protein)
					cmd.remove('%s and not chain a'%(protein))
					cmd.create('%s_active'%(protein),'resi %s in %s'%(select,protein))
					if dump == "true":
						cmd.save("%s.pdb"%protein,"%s"%protein)
					cmd.delete(protein)
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
			anchorset = []
			ResCount = {}
			for PROfile in ProfileSet:
				for res in PROfile:
					key = res[0]
					
					if res[1] == 0:
						anchorset.append(res)
							
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
						
			print ""	
				
			print "------------------------------"
			print "MOTIF STATS for %s_%s:"%(EC,entry)
			print "Different stats: %s across %s proteins"%(str(len(DistDict.keys())),len(ProteinSet))
			print "Anchor residue will have a Distance of 0 and an angle of 90"
			print "Residue Position, Residue, Distance from anchor, Angle, Count"
			
			for key in DistDict:
				if len(anchorset) > 0:
					anchormax = 0
					for res in anchorset:
						if ResCount[key] > anchormax:
							anchormax = ResCount[key]
							anchor = res
				consensusDic.append([[key[0],key[1]],DistDict[key],AngleDict[key],ResCount[key]])		
			
			for line in consensusDic:
				print "%s,%s,%s,%s,%s"%(line[0][0],line[0][1],line[1],line[2],line[3])
			
			model = UPDictionary[entry].split(",")[0]
			count = 0
			while True:
				selects = ActiveSites(model)
				if len(selects) != 0:
					
					stored.list = []
					cmd.fetch(model)
					cmd.remove('%s and not chain a'%(model))
					cmd.create('%s_active'%(model),'resi %s in %s'%(select[0],model))
					cmd.delete(model)
					cmd.iterate("(%s and name ca)"%(model),"stored.list.append((resi,resn))")
					if anchor[0] not in stored.list:
						model = UPDictionary[entry].split(",")[count]
						count = count + 1
						cmd.remove('all')
						continue
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
			
			if dump == "true":
				cmd.save("%s_%s_pre.pdb"%(EC,entry),"%s_active"%model)
			
			for key in consensusDic:
				if key[0] in modelSites:
					adjustRes(anchor, key[0], key[1], key[2], model)
			
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
	
	if dump == "true":
		nullstr = ""
		for prot in NullList:
			cmd.fetch(prot)
			cmd.remove(prot)
			nullstr = prot+","+nullstr
		nullstr = nullstr[:-1]
		outfile.write("None: %s"%nullstr)
		outfile.close()
	
	cmd.remove("all")
	
	print "done!"
	
cmd.extend("GenerateMotif",generateMotif)
print "GenerateMotif EC, MinProt, MaxProt, pfam/uniprot, DataDump, OutputDirectory"

def EvalMotif(motif, db, dir):
	'''
	This method takes data dumped motifs, active sites, and protein structures 
	and measures alignments
	'''

	if str(db) not in "PFUP":
		print "Enter the database argument as UP for UniProt or PF for PFam"
		return
	cmd.cd(dir)
	file = open("output.txt",'r')
	outfile = open("%s_%s_test.txt"%(motif,db),'w')
	MotifDict = {}
	
	
	for line in file:
		print line
		mname = line.split(' ')[0]
		mname = mname[:-1]
		if mname != "None":
			mname = motif + "_" + mname
		protlist = line.split(' ')[1]
		if "\n" in protlist:
			protlist = protlist[:protlist.index("\n")]
		#protlist = protlist[:-2]
		MotifDict[mname] = protlist

	print MotifDict
	for key in MotifDict:
		if key == "None":
			continue
		try:
			cmd.load("%s.pdb"%key)
		except:
			print "No motifs found for %s"%(key)
			continue
		for OKey in MotifDict:
			protList = MotifDict[OKey].split(',')
			for prot in protList:
				if OKey == "None":
					cmd.load("%s.pdb"%prot)
					aln = cmd.align("%s"%key,"%s"%prot)
					print "%s,%s,%s,%f"%(key,OKey,prot,aln[0])
					outfile.write("%s,%s,%s,%f\n"%(key,OKey,prot,aln[0]))
					cmd.remove("%s"%prot)
				else:
					try:
						cmd.load("%s_active.pdb"%prot)
						aln = cmd.align("%s"%key,"%s_active"%prot)
						print "%s,%s,%s,%f"%(key,OKey,prot,aln[0])
						outfile.write("%s,%s,%s,%f\n"%(key,OKey,prot,aln[0]))
						cmd.remove("%s_active"%prot)
					except:
						print "No motifs found for %s"%(key)
						try:
							cmd.remove("%s_active.pdb"%prot)
						except:
							print "It's not your fault"
			
	cmd.remove("all")
	print "done!"
	file.close()
	outfile.close()

cmd.extend("EvalMotif",EvalMotif)
print "EvalMotif EC, db, DumpDirectory"

def run():
	ec_nums = ['2.4.2.4', '1.13.11.2', '3.4.21.100', '2.7.1.11', '2.4.2.14', '3.2.1.31', '4.1.1.41', '2.7.1.15', '3.2.1.35', '1.3.1.34', '4.1.1.48', '4.1.1.49', '2.3.1.15', '5.3.99.3', '1.1.5.2', '2.3.1.16', '4.2.1.94', '2.7.1.90', '1.1.1.158', '2.1.1.72', '3.1.11.2', '1.14.99.3', '1.3.99.1', '4.2.1.40', '4.2.2.1', '1.14.99.1', '3.4.21.62', '3.1.2.22', '2.3.1.74', '3.4.11.5', '3.2.1.133', '3.2.1.96', '1.6.99.1', '1.20.4.1', '2.4.1.38', '3.4.21.68', '2.5.1.7', '3.8.1.6', '3.8.1.5', '5.4.2.1', '2.5.1.3', '3.1.3.1', '3.1.3.2', '1.17.4.1', '2.1.3.3', '1.17.4.2', '2.6.1.62', '2.5.1.9', '1.1.99.28', '4.1.1.50', '3.6.1.7', '3.4.22.2', '4.2.1.20', '1.16.3.1', '2.7.1.69', '1.2.1.2', '5.1.3.20', '1.10.3.2', '3.1.1.29', '5.4.99.5', '4.1.99.4', '4.1.99.3', '5.4.99.1', '3.5.1.5', '3.1.1.6', '3.1.1.7', '3.1.1.4', '1.7.1.1', '4.2.1.24', '3.1.1.1', '4.4.1.1', '5.1.1.3', '3.6.1.29', '5.1.1.7', '2.4.2.1', '3.4.23.4', '1.1.3.38', '2.3.2.10', '2.6.1.57', '3.4.23.1', '3.2.1.10', '3.4.19.12', '3.2.1.14', '3.2.1.17', '3.4.21.4', '3.1.4.11', '3.4.21.7', '4.1.2.18', '5.3.4.1', '1.8.1.9', '2.7.7.12', '3.8.1.2', '2.4.1.1', '3.2.1.91', '2.1.4.2', '4.1.1.23', '3.1.25.1', '4.2.1.11', '2.3.1.118', '2.3.1.28', '4.6.1.13', '2.1.1.13', '3.1.3.57', '3.1.30.2', '1.1.99.18', '4.2.1.52', '5.1.3.14', '3.4.22.27', '3.1.3.33', '6.3.1.1', '5.3.3.10', '5.1.3.13', '5.3.3.8', '3.4.22.39', '3.4.21.87', '2.1.2.8', '5.3.3.1', '5.4.99.11', '6.3.3.3', '2.1.2.2', '5.4.3.8', '3.1.27.5', '1.3.1.20', '3.2.1.45', '3.5.1.1', '1.14.13.39', '6.3.5.4', '3.1.27.3', '1.3.1.26', '4.1.3.7', '3.4.23.24', '4.1.3.3', '3.4.23.20', '4.1.3.1', '3.4.23.22', '3.1.1.47', '2.8.2.4', '5.5.1.6', '3.5.2.6', '3.1.1.43', '3.1.3.48', '3.1.1.41', '1.1.1.252', '1.5.1.34', '3.2.1.68', '4.1.1.39', '3.2.2.22', '2.4.1.25', '1.18.1.2', '6.3.2.3', '4.2.1.47', '2.5.1.19', '2.1.1.45', '2.1.4.1', '2.3.1.47', '2.1.1.48', '5.3.1.1', '2.6.1.1', '2.5.1.15', '5.3.1.5', '1.6.1.2', '3.5.1.77', '1.1.1.28', '3.4.21.92', '4.4.1.14', '4.1.2.13', '2.7.7.3', '3.2.1.129', '3.2.1.73', '2.7.1.40', '3.1.1.3', '5.3.99.2', '3.4.16.5', '3.4.22.40', '3.1.3.46', '2.3.1.41', '4.2.1.92', '2.2.1.2', '1.5.1.28', '1.11.1.6', '1.11.1.5', '3.2.1.113', '1.11.1.1', '1.2.7.1', '3.4.22.53', '3.1.4.4', '3.2.1.8', '3.1.27.4', '3.1.3.16', '2.4.1.207', '3.2.1.1', '1.1.1.35', '3.4.17.11', '1.1.1.37', '3.2.1.78', '1.2.1.12', '1.2.1.11', '6.3.2.1', '4.1.2.25', '4.2.1.3', '6.3.2.4', '4.2.1.1', '1.1.3.15', '3.1.27.1', '6.3.2.9', '1.11.1.10', '1.11.1.11', '3.2.3.1', '1.3.1.24', '2.3.1.9', '6.3.4.5', '2.6.1.21', '5.2.1.8', '4.2.2.5', '3.4.13.21', '4.99.1.1', '5.4.2.9', '3.1.1.61', '2.4.2.22', '2.3.1.54', '3.4.23.25', '2.3.3.9', '1.1.1.49', '3.7.1.8', '3.4.22.37', '2.7.7.48', '2.3.3.1', '2.4.1.44', '3.2.2.20', '2.7.4.6', '1.5.8.2', '1.1.2.3', '1.1.1.29', '4.1.2.39', '2.1.1.63', '1.14.17.3', '1.1.1.22', '1.1.1.21', '4.1.1.20', '1.1.1.27', '2.3.1.87', '2.7.4.3', '1.4.1.1', '2.1.1.113', '3.7.1.9', '3.4.23.23', '3.5.1.45', '3.2.1.52', '1.14.15.1', '3.5.3.15', '3.4.22.16', '2.7.6.3', '2.5.1.2', '3.5.4.4', '2.6.1.16', '3.4.23.15', '6.1.1.11', '2.7.1.38', '3.5.99.6', '1.2.1.59', '5.5.1.5', '4.1.1.64', '2.4.2.36', '6.1.1.19', '6.1.1.18', '2.7.1.37', '3.1.21.1', '2.4.2.31', '2.4.1.19', '3.4.21.26', '2.3.1.39', '1.2.1.8', '3.1.1.11', '3.1.1.13', '3.4.23.21', '3.5.1.11', '5.5.1.1', '2.4.1.10', '6.3.4.4', '3.4.21.43', '3.1.4.8', '4.2.1.10', '3.4.11.1', '2.3.1.97', '3.2.1.89', '4.1.2.20', '4.2.1.17', '1.14.13.7', '3.4.24.69', '2.3.2.13', '1.1.1.1', '1.1.1.3', '3.1.3.11', '1.1.1.8', '2.5.1.61', '3.2.1.21', '2.7.2.1', '4.2.3.3', '1.5.1.3', '4.2.3.9', '1.5.1.6', '4.2.1.60', '5.99.1.2', '5.99.1.3', '1.15.1.1', '4.2.99.18', '5.1.3.3', '5.1.3.1', '3.4.11.19', '3.1.27.10', '3.3.2.8', '1.14.19.2']
	for ec in ec_nums:
		generateMotif(ec,1,1000,'uniprot','false','C:/Users/Cameron/Documents/Pymol')
		generateMotif(ec,1,1000,'pfam','false','C:/Users/Cameron/Documents/Pymol')
cmd.extend("run",run)