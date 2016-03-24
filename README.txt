ProcessProtein.py contains all of the scripts related to distance based motif characterization

To run, place ProcessProtein.py within you PyMOL extensions folder at the following (or wherever you have it installed)

C:\Program Files\PyMOL\PyMOL\modules\pmg_tk\startup

Assuming no error, a message will pop up in PyMOL indicating the arguments for GenerateMotif. They are as follows

EC
	the ec number to base the motifs

MinProt
	takes a number
	the minimum number of proteins in a cluster to base a motif off of
	example: 1.2.3.4 with uniprot ID UP5678 may only have 4 proteins associated with it

MaxProtTotal
	takes none or a number
	the maximum proteins allowed in the analysis
	example: 2.7.11.1 has over 1600 proteins associated with it, maybe set the cap lower to make run time reasonable

db
	pfam or uniprot
	defaults to uniprot unless pfam is typed

OutDir
	The output directory for the motifs. Windows does not like writing to C:\Program Files so I recommend
	making a folder elsewhere or running PyMOL with administrator privileges 


Other noteworthy commands for PyMOL include 

ActiveSites protein
	Retrieves the active sites for the given protein and saves that as a PyMOL object
