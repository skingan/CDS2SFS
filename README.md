# CDS2SFS
Calculates site frequency spectra for synonymous and nonsynonymous sites in coding sequence alignment

# Author:
Sarah B. Kingan
University of Rochester
17 March 2014
last update 26 March 2015
	
# Input: 
one or more files containing multiple alignments of CDS in fasta format.

# Formatting requirements for alignment:
1. sequence must be coding (length is multiple of 3, no premature stop codons)
2. sequences may include terminal stop codon or not, but must be uniform across sequences
3. all sequences must be the same length
4. single outgroup sequence identified with 'outroup' in header line
5. files must have ".fa" extension


#Output:
			
1. STOUT
tab delimited text with the following fields (where N is the number of ingroup samples):

		filename_Syn S1 S2 ... SN
		
		filename_NonSyn S1 S2 ... SN
		
2. STDERR
progress reports, alignment information, details about excluded sites

# Example Calls		

EXAMPLE CALL ON SINGLE FILE: CDS2SFS.pl mydata.fa > mydata_sfs.out 2> mydata_sfs.err

EXAMPLE CALL ON ALL FILES IN DIRECTORY: CDSS2SFS.pl > mydata_sfs.out 2> mydata_sfs.err


# Exclusion of sites/handling of missing data
1. codons with any character other than AGCT in any sample are excluded in all samples
2. codons with sites that have more than 2 states in the ingroup sequences are excluded (violation of infinite sites model)
3. codons where ancestral state of the ingroup sequences cannot be inferred are ommited 
4. codons where the ingroup has more than two states are excluded

# Revision History

21 April 2014, modified to print out per gene stats - Sarah B. Kingan

22 January 2015, modified to be generalized to any dataset - Anthony J. Geneva

23 January 2015, modified to take either a single file or all files in directory - Sarah B. Kingan

27 January 2015, modified to accept files with no terminal stops and to include codons where neither ingroup codon matches outgroup codon but ancestral codon can be inferred. - Sarah B. Kingan

26 March 2015, commented and reformated for submission to github - Sarah B. Kingan
