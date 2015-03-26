#!/usr/bin/perl -w
use strict;
use warnings;

use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::AlignIO;
use Cwd;

####################################################################################################
#
#		Sarah B. Kingan
#		University of Rochester
#		17 March 2014
#		last update 26 March 2015
#	
#		Input: one or more files containing multiple alignments of CDS in fasta format.
#
#		Output:
#			STOUT
#				tab delimited text with the following fields (where N is the number of ingroup samples):
#					filename_Syn S1 S2 ... SN
#					filename_NonSyn S1 S2 ... SN
#			STDERR
#				progress reports, alignment information, details about excluded sites
#		
#		EXAMPLE CALL ON SINGLE FILE: CDS2SFS.pl mydata.fa > mydata_sfs.out 2> mydata_sfs.err
#		EXAMPLE CALL ON ALL FILES IN DIRECTORY: CDSS2SFS.pl > mydata_sfs.out 2> mydata_sfs.err
#
#		Formatting requirements for alignment:
#			sequence must be coding (length is multiple of 3, no premature stop codons)
#			sequences may include terminal stop codon or not, but must be uniform across sequences
#			all sequences must be the same length
#			single outgroup sequence identified with 'outroup' in header line
#			files must have ".fa" extension
#
#		Exclusion of sites/handling of missing data
#			codons with any character other than AGCT in any sample are excluded in all samples
#			codons with sites that have more than 2 states in the ingroup sequences are excluded (violation of infinite sites model)
#			codons where ancestral state of the ingroup sequences cannot be inferred are ommited 
#			codons where the ingroup has more than two states are excluded
#
#		21 April 2014, modified to print out per gene stats - Sarah B. Kingan
#		22 January 2015, modified to be generalized to any dataset - Anthony J. Geneva
#		23 January 2015, modified to take either a single file or all files in directory - Sarah B. Kingan
#		27 January 2015, modified to accept files with no terminal stops and to include 
#			codons where neither ingroup codon matches outgroup codon but ancestral codon
#			can be inferred. - Sarah B. Kingan
#		26 March 2015, commented and reformated for submission to github - Sarah B. Kingan
#
####################################################################################################


my $usage = "CDS2SFS.pl <OPT:file.fa>\n";


# ARRAY OF FILES TO ANALYZE
# MAY BE SINGLE FILE IF SPECIFIED AT COMMAND LINE
my @sorted_fasta_files;
my @fasta_files_list;
if ($ARGV[0]) {
	my $infile = $ARGV[0];
	push(@sorted_fasta_files, $infile)
}
else {
	@fasta_files_list = glob("*.fa*");
	@sorted_fasta_files = sort @fasta_files_list;
}



# MAIN LOOP FOR FILES ARRAY
my $file_count = 0;
for (my $f = 0; $f < scalar(@sorted_fasta_files); $f++) {
	$file_count++;
	my %seq_hash;
	my $infile = $sorted_fasta_files[$f];
	my $seqioobj = Bio::SeqIO->new(-file => "<$infile", -format => 'fasta');
	while (my $seqobj = $seqioobj->next_seq) {
		my $seq_name = $seqobj->display_id;
		$seq_hash{$seq_name} = $seqobj->seq;
	}
		
	
# QUALITY CONTROL 
#		checking seqs are same length
#		checking that seq lengths are multiples of 3
#		checking if seqs contain premature stop codons
#		checkin if all or no sequences have terminal stop codons
	print STDERR "Quality control details on $infile:\n";	
	my $QC = seqQC(\%seq_hash);
	if ($QC == '1') {
		print STDERR "\t$infile fails QC!\n";
		next;
	}
	my $terminal_stop = terminalStop(\%seq_hash);
	
	my $total_seq_number = keys(%seq_hash);
	if ($terminal_stop == $total_seq_number) { # all sequences have terminal stop: OK!
	}
	elsif ($terminal_stop == 0) { # no sequences have terminal stop: OK!
	}
	else { # some sequences have terminal stop: NO WAY! You make a mistake!
		print STDERR "\tnot all seqs end with stop codon!\n";
		next;
	}	
	print STDERR "\t$infile passes the QC!\n";
	
	
# ALIGNMENT DETAILS
# 		make ingroup seq hash
# 		identify outgroup seq
	print STDERR "Alignment details:\n";
	my %ingroup_seq_hash;
	my $outgroup_key;
	foreach my $key (sort keys %seq_hash) {
		unless ($key =~ /(\S)*outgroup(\S)*/){
			$ingroup_seq_hash{$key} = $seq_hash{$key};
		}
		if ($key =~ /(\S)*outgroup(\S)*/){
			$outgroup_key = $key;
			print STDERR "\t$outgroup_key is the outroup.\n";
		}
	}	
	my @ingroup_exemplar = keys %ingroup_seq_hash;
	my $new_hash_size = keys %ingroup_seq_hash;
	print STDERR "\talignment includes $new_hash_size ingroup individuals.\n";
	print STDERR "\t$outgroup_key is the outgroup.\n";
	print STDERR "\t$ingroup_exemplar[1] used as ingroup exemplar.\n";

# INTIALIZE SFS	
	my $ingroup_N = keys %ingroup_seq_hash;
	$ingroup_N--;
	my @cumulative_S_sfs = (0)x$ingroup_N;
	my @cumulative_N_sfs = (0)x$ingroup_N;
	
	
# IDENTIFY "BAD SITES" THAT VIOLATE INFINITE SITES MODEL FOR INGROUP OR HAVE MISSING DATA IN FULL DATASET
	my @bad_site_array = bad_sites($ingroup_exemplar[1], %seq_hash);
	
	
# CALCULATE SFS
#		sfs for synonymous (S) and nonsynonymous (N) changes,
#		including fixed class (divergent sites)
	
	my %unique_codon_hash;
	my $codon;
	my $unique_codon_count;
	my @unique_codon_array;
	my @unique_codon_freqs;
	my @aa_array;
	my $variable_position_count;
	my @S_N_counts;
	# ARRAY LENGTH IS EQUAL TO SAMPLE SIZE OF INGROUP SEQUENCES - 1
	my @S_sfs = (0)x$ingroup_N;
	my @N_sfs = (0)x$ingroup_N;
	my $dS = 0;
	my $dN = 0;
	# OUTGROUP
	my $outgroup_codon;
	my $outgroup_aa;
	my $ancestral_codon;
	my $ancestral_aa;
	my $original_outgroup_codon;
	my $original_outgroup_aa;
	my $which_ancestor;
	my $variable_position_count_from_outgroup;
	my $variable_position_count_from_ancestor;
	my $use_ancestor;
	my @dS_dN_counts;
	my $bad_site;
	my $count_of_violations = 0;
	my $codon_count = 0;
	my $file_name = $infile;
		
# LOOP FOR EACH CODON IN ALIGNMENT
# PICK A SINGLE SAMPLE FROM INGROUP SEQUENCES
	for (my $i = 0; $i <length($ingroup_seq_hash{$ingroup_exemplar[1]}); $i+=3) {
		my $codon_start_position = $i+1;
		%unique_codon_hash = ();
		$unique_codon_count = 0;
		$variable_position_count = 0;
		$outgroup_codon = '';
		$bad_site = 0;
		$which_ancestor = ''; # "outgroup" "inferred" or "NA"
	
# identify codons with IS violations or missing data
		foreach my $site (@bad_site_array) {
			my $diff = $site - $i;
			if (($diff < 3) && ($diff >=0)) {
				$bad_site++;
			}
		}
# analyse sites that do not have violations
		if ($bad_site == 0) {	
		
# get filename identifier		
			$file_name =~ s/.fa(\S)*//;
			
# get outgroup codon
			$outgroup_codon = substr($seq_hash{$outgroup_key},$i,3);
			$outgroup_aa = dna2aa($outgroup_codon);
# make hash where keys are unique ingroup codons sequence and values are the frequency (count) for that codon
			foreach my $key (keys %ingroup_seq_hash) {
				$codon = substr($ingroup_seq_hash{$key},$i,3);
				if (exists($unique_codon_hash{$codon})) {
					$unique_codon_hash{$codon}++;
				}
				else {
					$unique_codon_hash{$codon} = 1;
				}
			}
				
# make array of unique ingroup codons
# make array of frequency of those codons
# make array of amino acid translation for each codon
			@unique_codon_array = keys(%unique_codon_hash);
			@unique_codon_freqs = values(%unique_codon_hash);
			@aa_array = ();
			foreach my $element (@unique_codon_array) {
				push(@aa_array, dna2aa($element));
			}

# get count of unique in group codons
# get count of the number of sites that differ between ingroup codons
			$unique_codon_count = scalar(keys(%unique_codon_hash));
			$variable_position_count = count_variable_positions(@unique_codon_array);
	
# INVARIANT INGROUP CODONS: DIVERGENCE ONLY
# count substitutions
			if ($unique_codon_count == 1) {	
				$codon_count++;
# outgroup is different
				if ($outgroup_codon ne $unique_codon_array[0]) {
					$variable_position_count_from_outgroup = count_diffs($outgroup_codon, $unique_codon_array[0]);
# outgroup is different by one position...
					if ($variable_position_count_from_outgroup == 1) {
						if (dna2aa($outgroup_codon) eq dna2aa($unique_codon_array[0])) {
							$dS++;
						}
						else {
							$dN++;
						}
					}
# ...by 2 positions
					if ($variable_position_count_from_outgroup == 2) { 
						@dS_dN_counts = two_changes($outgroup_codon, $unique_codon_array[0]);
						$dS += $dS_dN_counts[0];
						$dN += $dS_dN_counts[1];
					}
# ...by 3 positions
					if ($variable_position_count_from_outgroup == 3) { 
						@dS_dN_counts = three_changes($outgroup_codon, $unique_codon_array[0]);
						$dS += $dS_dN_counts[0];
						$dN += $dS_dN_counts[1];
					}
				}
			}	
				
					
# POLYMORPHIC SITES: TWO INGROUP CODONS
			if ($unique_codon_count == 2) {
			
# determine what the ancestral codon is: equal to outgroup, inferred with parsimony, or NA
				if (($outgroup_codon eq $unique_codon_array[0]) || ($outgroup_codon eq $unique_codon_array[1])) {
					$which_ancestor = "outgroup";
				}
				else {
					$ancestral_codon = inferAncestor($outgroup_codon, $unique_codon_array[0], $unique_codon_array[1]);
					if ($ancestral_codon !~ /error/) {
						$which_ancestor = "inferred";
						$original_outgroup_codon = $outgroup_codon;
						$original_outgroup_aa = $outgroup_aa;
						$outgroup_codon = $ancestral_codon;
						$outgroup_aa = dna2aa($outgroup_codon);
					}
					else {
						$which_ancestor = "NA";
					}
				}
# codons whose ancestral state can be inferred
				if ($which_ancestor ne "NA") {
					$codon_count++;
# SFS FOR INGROUP CODONS
# codons differ by 1 position...
					if ($variable_position_count == 1) {
						for (my $i = 0; $i < scalar@unique_codon_array; $i++) {
							if ($outgroup_codon ne $unique_codon_array[$i]) {
								if ($outgroup_aa eq $aa_array[$i]) { # synonymous change
									$S_sfs[$unique_codon_freqs[$i]-1]++;
								}
								else { 								# non synonymous change
									$N_sfs[$unique_codon_freqs[$i]-1]++;
								}
							}
						}
					}
# ...by 2 positions...
					elsif ($variable_position_count == 2) {
						for (my $i = 0; $i < scalar@unique_codon_array; $i++) {
							if ($outgroup_codon ne $unique_codon_array[$i]) {
								@S_N_counts = two_changes(@unique_codon_array);
								$S_sfs[$unique_codon_freqs[$i]-1] += $S_N_counts[0];
								$N_sfs[$unique_codon_freqs[$i]-1] += $S_N_counts[1];			
							}
						}
					}
# ...by 3 positions...
					elsif ($variable_position_count == 3) {
						for (my $i = 0; $i < scalar@unique_codon_array; $i++) {
							if ($outgroup_codon ne $unique_codon_array[$i]) {
								@S_N_counts = three_changes(@unique_codon_array);
								$S_sfs[$unique_codon_freqs[$i]-1] += $S_N_counts[0];
								$N_sfs[$unique_codon_freqs[$i]-1] += $S_N_counts[1];			
							}
						}
					}
				}
# DIVERGENT SITES ALONG OUTGROUP BRANCH
				if ($which_ancestor eq "inferred") {
					$codon_count++;
					$variable_position_count_from_outgroup = count_diffs($original_outgroup_codon, $ancestral_codon);
# outgroup is different by one position...
					if ($variable_position_count_from_outgroup == 1) {
						if (dna2aa($original_outgroup_codon) eq dna2aa($ancestral_codon)) {
							$dS++;
						}
						else {
							$dN++;
						}
					}
# ...by 2 positions
					if ($variable_position_count_from_outgroup == 2) { 
						@dS_dN_counts = two_changes($original_outgroup_codon, $ancestral_codon);
						$dS += $dS_dN_counts[0];
						$dN += $dS_dN_counts[1];
					}
# ...by 3 positions
					if ($variable_position_count_from_outgroup == 3) { 
						@dS_dN_counts = three_changes($original_outgroup_codon, $ancestral_codon);
						$dS += $dS_dN_counts[0];
						$dN += $dS_dN_counts[1];
					}	
				}
# ANCESTOR CANNOT BE INFERRED
				if ($which_ancestor eq "NA") {
					print STDERR "#####\t$codon_start_position: no ingroup codon matches outgroup, ancestor CANNOT be inferred\n";
					$count_of_violations++;
				}
			}
# MORE THAN TWO INGROUP CODONS, SKIP THIS CODON!
			elsif ($unique_codon_count > 2) {
				print STDERR "#####\t$codon_start_position: more than 2 ingroup codons\n";
				$count_of_violations++;
			}		
		}
# SITES THAT VIOLATE INFINITE SITES IN INGROUP SEQUENCES OR HAVE MISSING DATA IN FULL DATASET
		else {
			print STDERR "#####\t$codon_start_position: violates IS or has missing data\n";
			$count_of_violations++;
		}
	}
	

# PRINT ACCOUNTING OF CODONS IN ALIGNMENT
	print STDERR "Analysis details:\n";
	print STDERR "\tTotal excluded codons: ", $count_of_violations, "\n";
	print STDERR "\tTotal analyzed codons: ", $codon_count, "\n";

# PRINT RESULTS
	print $file_name, "_Syn\t";
	print join("\t", @S_sfs), "\t", $dS, "\t","\n";
	print $file_name, "_NonSyn\t";
	print join("\t", @N_sfs), "\t", $dN,"\n";
	
}	

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#																	<
#>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<SUBROUTINES<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#																	<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##########################
sub dna2aa {
##########################
	my($seq) = @_;
	my $DNAseq_obj = Bio::Seq->new(-seq => $seq);
	my $aa_obj = $DNAseq_obj->translate;
	my $aa_seq = $aa_obj->seq;
	return $aa_seq;
}

# check that seqs of equal length with no internal stop codons and lengths are multiple of 3
##########################
sub seqQC {
##########################
	my $hash_ref = shift;
	my %hash = %$hash_ref;
	my $N_error = 0;
	my $length_error = 0;
	my $stop_error = 0;
	my $triplet_error = 0;
	my $N = scalar(keys %hash);
	my @length_array;
	my $aa;
	my $mod;
# are the seqs the same length?
# are there internal stop codons?
	foreach my $key (keys %hash) {
		push (@length_array, length($hash{$key}));
		$aa = dna2aa($hash{$key});
		if (substr($aa, 0, -1) =~ /\*/) {
			print STDERR "\t", $key, " has a premature stop!\n";
			$stop_error++;
		}
	}
	for (my $i = 0; $i < (scalar@length_array-1); $i++) {
		for (my $j = 1; $j < scalar@length_array; $j++) {
			if ($length_array[$i] ne $length_array[$j]) {
				$length_error++;
			}
		}
	}
	if ($length_error > 0) {
		print STDERR "\tnot all sequences are the same length!\n";
	}
	if ($length_error == 0) {
		$mod = ($length_array[0] % 3);
		if ($mod != 0) {
			$triplet_error++;
			print STDERR "\tsequence length not multiple of 3!\n";
		}
	}
	if (($N_error > 0) || ($length_error > 0) || 
		($stop_error > 0) || ($triplet_error > 0)) {
		return '1';
	}
}

##########################
sub terminalStop {
##########################
	my $hash_ref = shift;
	my %hash = %$hash_ref;
	my $aa;
	my $stop_count = 0;
	foreach my $key (keys %hash) {
		$aa = dna2aa($hash{$key});
		if (substr($aa, -1, 1) eq '*') {
			$stop_count++;
		}
	}
	return $stop_count;
}

# find sites with missing data or that violate infinite sites model
###########################
sub bad_sites {
###########################
	my $in_key = shift;
	my %hash = @_;
	my @site_array = ();
	my @unique_site_array = ();
	my %seen = ();
	my @uniq = ();
# ANCHOR BY SINGLE WITHIN-SPECIES SAMPLE
	for (my $i = 0; $i<length($hash{$in_key}); $i++) { 
		my @all_base_array = ();
		my @ingroup_base_array = ();
		foreach my $seq (keys %hash) {
			push(@all_base_array, substr($hash{$seq}, $i, 1));
			unless ($seq =~ /outgroup/) {
				push(@ingroup_base_array, substr($hash{$seq}, $i, 1));
			}
		}
		%seen = ();
		@uniq = ();
# check for missing data, all seqs
		foreach my $base (@all_base_array) {
			push(@uniq, $base) unless $seen{$base}++;
			if ($base !~ /A|C|G|T|a|c|t|g/) {
				push(@site_array, $i);
			}
		}
		%seen = ();
		@uniq = ();
		foreach my $base (@ingroup_base_array) {
			push(@uniq, $base) unless $seen{$base}++;
		}	
		if (scalar(@uniq) > 2) {
			push(@site_array, $i);
		}
	}
	my %tmp_hash;
	foreach my $site (@site_array) {
		$tmp_hash{$site} = '1';
	}
	@unique_site_array = sort keys %tmp_hash;
	return @unique_site_array;
}

# rounds numbers in array
##########################
sub round_array {
##########################
	my (@array) = @_;
	my $rounded;
	my @rounded_array;
	foreach my $element (@array) {
		$rounded = int($element + 0.5);
		push(@rounded_array, $rounded);
	}
	return @rounded_array;
}

# count number of differences between codon pair
##########################
sub count_diffs {
##########################
	my (@array) = @_; # codon pair
	my $count = 0;
	for(my $i = 0; $i<length$array[0]; $i++) {
		if (substr($array[0], $i, 1) ne substr($array[1], $i, 1)) {
			$count++;
		}
	}
	return $count;
}

# count number of differences between array of codons
##########################
sub count_variable_positions {
##########################
	my (@array) = @_;
	my %hash;
	my $count = 0;
	for (my $i = 0; $i<3; $i++) {
		%hash = ();
		for (my $j = 0; $j<scalar(@array); $j++) {
			 $hash{substr($array[$j],$i,1)} = 1;
		}
		if (scalar(keys(%hash))>1) {
			$count++;
		}
	}
	return $count;
}

# input array of 3 codons, the first of which is for outgroup
# return ancestral codon for ingroup
##########################
sub inferAncestor {
##########################
	my (@array) = @_;
# OUTGROUP
	my $outgroup_codon = shift(@array);
	my @ingroup_codons = @array;
	my $anc_codon = '';
	my @ingroup_base_array;
	my @unique_ingroup_base_array;
	my $base;
	for (my $i = 0; $i<3; $i++) {
		@ingroup_base_array = ();
		foreach my $codon (@ingroup_codons) {
			$base = substr($codon,$i,1);
			push(@ingroup_base_array, $base);
		}
		@unique_ingroup_base_array = uniq(@ingroup_base_array);
		if (scalar@unique_ingroup_base_array == 1) {
			$anc_codon .= $ingroup_base_array[0];
		}
		elsif (($ingroup_base_array[0] eq substr($outgroup_codon, $i, 1)) ||
				($ingroup_base_array[1] eq substr($outgroup_codon, $i, 1))) {
			$anc_codon .= substr($outgroup_codon, $i, 1);
		}
		else {
			$anc_codon = 'error';
			return $anc_codon;
		}
	}
	return $anc_codon;
}

# input array
# return array og unique sorted elements
##########################
sub uniq {
##########################
	my @input = @_;
	my %hash;
	foreach my $a (@input) {
		$hash{$a} = '1';
	}
	my @output = (sort keys %hash);
	return @output;
}



# input array of codons, the first is the ancestor, the second is outgroup
# return array of codons encoded as 0 and 1 where 1 is derived
##########################
sub nt2binary {
##########################
	my (@array) = @_;
	my $anc_codon = shift(@array);
	my @binary_array = ('000');
	foreach my $codon (@array) {
		my $binary_codon = '';
		for (my $i=0; $i<3; $i++) {
			my $anc_base = substr($anc_codon, $i, 1);
			if (substr($codon, $i, 1) eq $anc_base) {
				$binary_codon .= '0';
			}
			else {
				$binary_codon .= '1';
			}
		}
		push (@binary_array, $binary_codon);
	}
	return @binary_array;
}

# count S and N between pair of codons that differ at two positions
##########################
sub two_changes {
##########################
	my (@array) = @_; # codon pair
	my $codon0 = $array[0];
	my $codon1 = $array[1];
	my @positions;					# positions that differ between codons
	my $S = 0;
	my $N = 0;
	my @annotation_array;
	my $stop;
 
# determine positions that differ between codons    
	for (my $i = 0; $i < 3; $i++) {
		if (substr($codon0, $i, 1) ne substr($codon1, $i, 1)) {
			push(@positions, $i);
		}
	}
	
# make intermediate codons and arrays of terminal and intermediate codons
	my $intermedA = $codon0;
	my $intermedB = $codon0;	
	substr($intermedA, $positions[0], 1) = substr($codon1, $positions[0], 1);
	substr($intermedB, $positions[1], 1) = substr($codon1, $positions[1], 1);
	my @codons = ($codon0, $codon1);
	my @intermeds = ($intermedA, $intermedB);

# ignore paths through stops
	$stop = 0;
	foreach my $c (@intermeds) {
		if (dna2aa($c) eq '*') {
			$stop++;
		}
	}
	if ($stop == 0) {
	
# determine types of changes for both 2-step paths
		foreach my $i (@codons) {
			foreach my $j (@intermeds) {
				if (dna2aa($i) eq dna2aa($j)) {
					$S++;
				}
				else {
					$N++;
				}
			}
		}
	}
	
# return array of changes
	if (($S == 0) && ($N == 0)) {
		@annotation_array = (0,0)
	}
	else {
		@annotation_array = (2*$S/($S+$N), 2*$N/($S+$N));
	}
	return @annotation_array;
}
		
# count S and N between pair of codons that differ at three positions
##########################
sub three_changes {
##########################
	my (@array) = @_; # array of two codon seqs that differ at three bases
	my @codon = (0,0,0,0);
	$codon[0] = $array[0];
	$codon[3] = $array[1];
	my $S = 0;
	my $N = 0;
	my $stop;
	my @annotation_array;
# each path of three changes	
	for (my $i = 0; $i < 3; $i++) {
		for (my $j = 0; $j < 3; $j++) {
			if ( $i != $j ) {
				$stop = 0;
# populate codon array with intermediates
				$codon[1] = $codon[0];
				substr($codon[1], $i, 1) = substr($codon[3], $i, 1);
				$codon[2] = $codon[1];
				substr($codon[2], $j, 1) = substr($codon[3], $j, 1);
# ignore paths through stops
				foreach my $c (@codon[1..2]) {
					if (dna2aa($c) eq '*') {
						$stop++;
					}
				}
				if ($stop == 0) {
# count S and N changes
					for (my $k = 0; $k < 3; $k++) {
						if (dna2aa($codon[$k]) eq dna2aa($codon[$k+1])) {
							$S++;
						}
						else {
							$N++;
						}	
					}
				}
			}
		}
	}
# return counts of changes
	if (($S == 0) && ($N == 0)) {
		@annotation_array = (0,0)
	}
	else {
		@annotation_array = (3*$S/($S+$N), 3*$N/($S+$N));
	}
	return @annotation_array;
}				

# calculate TD
###########################
sub sfs2TD {
###########################
	my @sfs = @_;
# SAMPLE SIZE
	my $n = scalar(@sfs) + 1;
	my $S = 0;
	my $a1 = 0;
	my $a2 = 0;
	my $b1;
	my $b2;
	my $e1;
	my $e2;
	my $D;
# total number of segregating sites
	for (my $i = 0; $i < scalar(@sfs); $i++) {
		$S += $sfs[$i];
	}
	if ($S == 0) {
		$D = "NAN";
	}
	else {
# constants a1 and a2
		for (my $i = 1; $i < $n; $i++) {
			$a1 += (1/$i);
			$a2 += (1/($i*$i));
		}
# constants b1, b2, e1, e2
		$b1 = ($n + 1)/(3*($n-1));
		$b2 = (2*($n*$n + $n + 3))/(9*$n*($n-1));
		$e1 = ($b1 - (1/$a1))/$a1;
		$e2 = ($b2 - ($n+2)/($a1*$n) + $a2/($a1*$a1))/($a1*$a1 + $a2);
# Tajima's D
		for (my $i = 1; $i < scalar(@sfs); $i++) {
			$D += $sfs[$i-1] * (((2*$i)*($n-$i))/($n*($n-1)) - 1/$a1)/sqrt($e1*$S + $e2*$S*($S-1));
		}
	}
	return $D;
}

exit;
