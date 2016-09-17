The current goal for updated GenerateMotifs is to allow multiple sites to come from a single protein 
and to resolve issues caused by certain ECs

In the current version, all amino acids found within active sites as per the CSA are exported as a single
active site. Multiple active sites represented as a single site is not representative of the nature of an
enzyme.

To combat this, I have made changes within ActiveSites and the 443 - 462 region of GenerateMotif to export
each found active site as a potential contributor to the motif. This part is functional as far as I know.

The next step is to find the most common residues within the list of amino acids belonging to an active sites.
By filtering out outliers in terms of amino acid counts, we can filter out the access active sites found by
ActiveSites. This is the current place to work from, located in lines 506 - 545. 

From there, the algorithm looks for a protein with only the amino acids located within the consensus site
and adjusts those amino acids from there. Small adjustments may need to be made there (only in places
with ActiveSites but still)

The second major problem that needs to be resolved pops up when running motifs for 1.3.11.2. PyMOL throws 
a critical error and stops running. I think it is something related to atom selection but idk. If this 
problem is resolved everything should run smoothly

I also made a run function for PyMOL that will auto generate the motifs

Changes aren't too hard. It will just take print statements and time