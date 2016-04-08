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

DataDump
	true or false
	defaults to false if true is not typed
	Writes a lot more information, namely keeping pdb files of all active sites and proteins without active
	sites. Also keeps the full proteins of the proteins without database designation
	A motif must be generated using this function to be used in evaluation

OutDir
	The output directory for the motifs. Windows does not like writing to C:\Program Files so I recommend
	making a folder elsewhere or running PyMOL with administrator privileges 


Other noteworthy commands for PyMOL include 

ActiveSites protein
	Retrieves the active sites for the given protein and saves that as a PyMOL object

EvalMotif EC,  DumpDirectory
	EC - EC number
	DumpDirectory - the directory where all of the information for a given motif was dumped to
			see DataDump argument in GenerateMotif
	
	The output of EvalMotif is a CSV denoting alignments between all motifs and all proteins within the 
	EC class. If there were no motifs found or not enough proteins to generate a motif for a given class,
	that motif will be skipped. 

	The CSV format is: 
	motif, motif designation of the target, the pdbid of the target, the RMSD of the alignment

	Ex.
	1.1.1.1_P11766,1.1.1.1_P11766,2FZW,0.281117
	1.1.1.1_P11766,1.1.1.1_P39462,1JVB,1.651022

stats.R
	This script contains a function that will flag motifs whose alignments
	are not significantly different then others as a way to check for false positives

	It also prints boxplots relating to each motif depending on alignment scores