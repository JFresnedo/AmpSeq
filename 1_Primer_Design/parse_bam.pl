#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Sam;

my ($bamfile, $genome_file, $marker_list, $maxpvalue)=@ARGV;

# high level API
my $sam = Bio::DB::Sam->new(-bam  =>$bamfile,
                             -fasta=>$genome_file,
                             );
open OUT, ">selected_sites";
print OUT "marker\trefallele\tnum_tags\tpvalue\tphase\treject_code\ta1\tcount\tlseq\tlseq_len\tlseq_pct\trseq\trseq_len\trseq_pct\ta2\tcount\tlseq\tlseq_len\tlseq_pct\trseq\trseq_len\trseq_pct\tall_cigars\tall_records\n";
open ERR, ">error";
$,="\t";
my $min_edge_size = 10;
my $max_match_length = 20;



open IN, "$marker_list";
<IN>;
LOOP1:while (<IN>) 
{
	chomp;
	my @data= split "\t";
	my ($marker, $pvalue, $phase) = @data[0, 2, 7];
	next LOOP1 if ($pvalue>$maxpvalue);
	$marker=~s/\s//g;
	my $chr = 0;
	my $pos = 0;
	if ($marker=~/^S([^_]+)_(\d+)/) 
	{
		$chr = "chr".$1;
		$pos = $2;
	}
	else
	{
		next LOOP1;
	}
	my $segment = $sam->segment($chr,$pos,$pos);
	unless (defined $segment) 
	{
		print "Warning: marker $marker cannot find tags ($chr $pos)! skipped\n.";
		next LOOP1;
	}
	my $iterator = $segment->features(-iterator=>1);
	my $allrecords = "";
	my $allcigars = "";
	my $refallele = "";
	my $num_tags = 0;
	my %allelerecords = ();
	my $reject_code =0;
	#my $badstr = "";
	LOOP2:while (my $align = $iterator->next_seq) 
	{
		$num_tags++;
		my $cigar = $align->cigar_str;
		$allcigars .= $cigar. "; ";
		my $r_start = $align->start;
		my $r_end = $align->end;
		my ($ref,$matches,$query) = $align->padded_alignment;

		if ($cigar=~/^(\d+)S/) 
		{
			my $cliplen = $1;
			$cigar=~s/${cliplen}S//;
			$ref=substr ($ref, $cliplen);
			$matches=substr ($matches, $cliplen);
			$query=substr ($query, $cliplen);
		}

		if ($cigar=~/(\d+)S$/) 
		{
			my $cliplen = $1;
			$cigar=~s/${cliplen}S//;
			$ref=substr ($ref, $cliplen);
			$matches=substr ($matches, $cliplen);
			$query=substr ($query, $cliplen);
		}


		if (($cigar=~/^\d+M$/) && ($pos>=$r_start) && ($pos<=$r_end))
		{
			my $allele = substr ($query, $pos-$r_start, 1);
			next LOOP2 if (! (defined $allele));
			$refallele = substr ($ref, $pos-$r_start, 1);
			my $qleft = substr ($query, 0, $pos-$r_start);
			my $qright = substr ($query, $pos-$r_start+1);

			my $rleft = substr ($ref,  0, $pos-$r_start);
			my $rright = substr ($ref, $pos-$r_start+1);

			my $lpct =0;
			my $rpct =0;
			if (length($qleft)>=$min_edge_size) 
			{
				my $t1 = $qleft;
				my $t2 = $rleft;
				if (length($t1)>$max_match_length) 
				{
					$t1 = substr($t1, length($t1) - $max_match_length);
					$t2 = substr($t2, length($t2) - $max_match_length);
				}
				$lpct = get_pct($t1, $t2);
			}
			if (length($qright)>=$min_edge_size) 
			{
				my $t1 = $qright;
				my $t2 = $rright;
				if (length($t1)>$max_match_length) 
				{
					$t1 = substr($t1, 0, $max_match_length);
					$t2 = substr($t2, 0, $max_match_length);
				}
				$rpct = get_pct($t1, $t2);
			}
			if (exists $allelerecords{$allele})
			{
				$allelerecords{$allele}{"count"}++;
			}
			else
			{
				$allelerecords{$allele}{"count"}= 1;
			}

			if (exists $allelerecords{$allele}{"lpct"}) 
			{
				if ($lpct >=$allelerecords{$allele}{"lpct"}  ) 
				{
					##lseq keep the best left end pct match
					$allelerecords{$allele}{"lpct"}  = $lpct;
					##lseq keep the best left end sequence
					if (($lpct >$allelerecords{$allele}{"lpct"}) || ((length $qleft) > (length $allelerecords{$allele}{"lseq"}) ))
					{
						$allelerecords{$allele}{"lseq"}  = $qleft;
					}	
				}
			}
			else
			{
				$allelerecords{$allele}{"lpct"}  = $lpct;
				$allelerecords{$allele}{"lseq"}  = $qleft;
			}

			if (exists $allelerecords{$allele}{"rpct"}) 
			{
				if ($rpct >=$allelerecords{$allele}{"rpct"}  ) 
				{
					##lseq keep the best left end pct match
					$allelerecords{$allele}{"rpct"}  = $rpct;
					##lseq keep the best left end sequence
					if (($rpct >$allelerecords{$allele}{"rpct"}) || ((length $qright) > (length $allelerecords{$allele}{"rseq"}) ))
					{
						$allelerecords{$allele}{"rseq"}  = $qright;
					}	
				}
			}
			else
			{
				$allelerecords{$allele}{"rpct"}  = $rpct;
				$allelerecords{$allele}{"rseq"}  = $qright;
			}
			$allrecords .=$allele. ":". ($pos-$r_start). ":$lpct-". ($r_end-$pos). ":$rpct ";
		}
		#else
		#{
		#	$badstr .="$cigar:" . ($r_end-$r_start+1). "-".(length $query)."; ";
		#}
		#my $cigar_array = $align->cigar_array;		
		#my $cigar_segments = @{$cigar_array};

		#print $marker, $cigar,  length($ref), $r_end-$r_start+1, "x${ref}x", "y${matches}y", "z${query}z";
		#for (my $i=0; $i<$cigar_segments; $i++) 
		#{
			#print $cigar_array->[$i]->[1], "";
		#}
		#print "\n";
		
		#$cigar_array->[0]->[1],"\n";
		#print $cigar, $r_start, $r_end, $matches, "\n";
	}
	my @alleles = keys %allelerecords;
	if (@alleles>2) 
	{
		$reject_code = 1;
	}
	elsif (@alleles==1)
	{
		$reject_code = 2;
		$alleles[1] = $alleles[0];
	}
	elsif (@alleles==0)
	{
		$reject_code = 6;
	}
	if ($reject_code<=2) 
	{
		if ($alleles[0] eq $refallele) 
		{
			##do nothing
		}
		elsif ($alleles[1] eq $refallele)  
		{
			##switch
			my $t = $alleles[0];
			$alleles[0] = $alleles[1];
			$alleles[1] = $t;	
		}
		else
		{
			$reject_code = 3;
		}
	}
	
	print OUT $marker, $refallele,  $num_tags, $pvalue, $phase, "";
	if (($reject_code==0) || ($reject_code==1) || ($reject_code==2))
	{
		my ($a1, $a2);
		$a1 = $alleles[0];
		if (@alleles>1) 
		{
			$a2 = $alleles[1];
		}
		else
		{
			$a2 = $alleles[0];
		}

		if (($allelerecords{$a1}{"lpct"} <= 90) || ($allelerecords{$a1}{"rpct"} <= 90) )
		{
			$reject_code = 4;
		}
		if (($allelerecords{$a2}{"lpct"} <= 90) || ($allelerecords{$a2}{"rpct"} <= 90) )
		{
			$reject_code = 5;
		}
		print ERR $reject_code, $allelerecords{$a2}{"rpct"}, "\n";
		print OUT $reject_code, $a1, $allelerecords{$a1}{"count"},
		  $allelerecords{$a1}{"lseq"},
		  length $allelerecords{$a1}{"lseq"},
		  $allelerecords{$a1}{"lpct"},
		   $allelerecords{$a1}{"rseq"},
		  length $allelerecords{$a1}{"rseq"},
		  $allelerecords{$a1}{"rpct"}, 
		  $a2, $allelerecords{$a2}{"count"},
		  $allelerecords{$a2}{"lseq"},
		  length $allelerecords{$a2}{"lseq"},
		  $allelerecords{$a2}{"lpct"},
		   $allelerecords{$a2}{"rseq"},
		  length $allelerecords{$a2}{"rseq"},
		  $allelerecords{$a2}{"rpct"}, "$allcigars\t$allrecords\n";
	}
	else
	{
		print OUT $reject_code, "";
		for (my $i=0; $i<16; $i++) 
		{
			print OUT "\t";
		}
		print OUT "$allcigars\t$allrecords\n";
	}

}

#### reject code
#1. more than 2 alleles
#2. less than 2 alleles
#3. no reference allele
#4. ref allele flanking failed
#5. alt allele flanking failed

sub get_pct
{
	my $a= shift;
	my $b= shift;
	my $len1 = length $a;
	my $len2 = length $b;
	unless ($len1 == $len2) 
	{
		return 0;
	}
	my $match =0;
	for (my $i=0; $i<$len1; $i++) 
	{
		if (substr($a, $i, 1) eq substr($b, $i, 1)) 
		{
			$match ++;
		}
	}
	return int($match/$len1*100+0.5);

}
l
