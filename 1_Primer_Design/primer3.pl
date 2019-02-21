#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my ($genome_file,  $markerlist, $primer_len)=@ARGV;


my $primer3_path = "/usr/local/bin/primer3_core";
my $primer3_config = "/usr/local/bin/primer3/src/primer3_config/";
my $resultfile = "primer_results.txt";
my $temp_primer3_output = "primer3_results";

open IN, "$markerlist";
my %markers = ();
my %genome = ();
<IN>;
while (<IN>) 
{
	my @data = split "\t";
	if ($data[4] <= 2) 
	{
		$data[0] =~/^S([^_]+)_(\d+)$/;
		my $chr=$1;
		my $pos = $2;
		push @{$markers{$chr}}, $pos; 

	}
}
close IN;
my $in  = Bio::SeqIO->new(-file => $genome_file ,
                     -format => 'Fasta');
while ( my $seq = $in->next_seq() ) 
{
  my $chr = $seq->display_id;
  $chr=~s/chr//;
  if (exists $markers{$chr}) 
	{
	  $genome{$chr} = $seq->seq();
	}
}
open OUT, "|$primer3_path > $temp_primer3_output";
foreach my $chr(keys %markers)
{
	my @markers = @{$markers{$chr}};
	foreach my $pos (@markers) 
	{
		my $primer1 = substr($genome{$chr}, $pos-$primer_len-1, $primer_len);
		my $primer2 = revdnacomp(substr($genome{$chr}, $pos, $primer_len));
		print OUT "SEQUENCE_ID=S${chr}_${pos}_left\n";
		print OUT "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$primer3_config\n";
		print OUT "PRIMER_TASK=check_primers\n";
		print OUT "PRIMER_EXPLAIN_FLAG=1\n";
		print OUT "SEQUENCE_PRIMER=$primer1\n";
		print OUT "PRIMER_PICK_ANYWAY=1\n";
		print OUT "=\n";
		print OUT "SEQUENCE_ID=S${chr}_${pos}_right\n";
		print OUT "PRIMER_TASK=check_primers\n";
		print OUT "PRIMER_EXPLAIN_FLAG=1\n";
		print OUT "SEQUENCE_PRIMER=$primer2\n";
		print OUT "PRIMER_PICK_ANYWAY=1\n";
		print OUT "=\n";
	}
}
close OUT;
open IN, $temp_primer3_output;
my %primer_info = ();
my $marker_id ="";
my $primer_id="";
while (<IN>) 
{
	chomp;
	if (/^SEQUENCE_ID=(\S+)_(left|right)/) 
	{
		$marker_id= $1;
		$primer_id = $2;
	}
	else
	{
		my ($key, $value)=split "=";
		$primer_info{$marker_id}{$primer_id}{$key} = $value;
 	}
}
close IN;

open IN, "$markerlist";
<IN>;
open OUT, ">$resultfile";
print OUT "marker\trefallele\tnum_tags\tpvalue\tphase\treject_code\ta1\tcount\tlseq\tlseq_len\tlseq_pct\trseq\trseq_len\trseq_pct\ta2\tcount\tlseq\tlseq_len\tlseq_pct\trseq\trseq_len\trseq_pct\tall_cigars\tall_records\tlprimer\tltemp\tlproblem\trprimer\trtemp\trproblem\n";
$,="\t";
LOOP:while (<IN>) 
{
	chomp;
	next LOOP unless (/\w/);
	my @data = split "\t";
	my $marker_id= $data[0];
	print OUT @data;

	my $lprimer = "";
	my $rprimer = "";
	
	my $lproblem = "";
	my $rproblem = "";

	my $ltemp = "";
	my $rtemp = "";

	
	
	if ($data[4] <=2) 
	{
		if (exists $primer_info{$marker_id}{"left"}{"SEQUENCE_PRIMER"}) {
		$lprimer = $primer_info{$marker_id}{"left"}{"SEQUENCE_PRIMER"};
		}
		else
		{
			$lprimer = "";
		}

		if (exists $primer_info{$marker_id}{"right"}{"SEQUENCE_PRIMER"}) {
		$rprimer = $primer_info{$marker_id}{"right"}{"SEQUENCE_PRIMER"};
		}
		else
		{
			$rprimer = "";
		}

		if (exists $primer_info{$marker_id}{"left"}{"PRIMER_LEFT_0_TM"}) {
		$ltemp = $primer_info{$marker_id}{"left"}{"PRIMER_LEFT_0_TM"};
		}
		else
		{
			$ltemp = "";
		}

		if (exists $primer_info{$marker_id}{"right"}{"PRIMER_LEFT_0_TM"}) {
		$rtemp = $primer_info{$marker_id}{"right"}{"PRIMER_LEFT_0_TM"};
		}
		else
		{
			$rtemp = "";
		}

		if (exists $primer_info{$marker_id}{"left"}{"PRIMER_LEFT_0_PROBLEMS"}) {
		$lproblem = $primer_info{$marker_id}{"left"}{"PRIMER_LEFT_0_PROBLEMS"};
		}
		else
		{
			$lproblem = "";
		}

		if (exists $primer_info{$marker_id}{"right"}{"PRIMER_LEFT_0_PROBLEMS"}) {
		$rproblem = $primer_info{$marker_id}{"right"}{"PRIMER_LEFT_0_PROBLEMS"};
		}
		else
		{
			$rproblem = "";
		}

		print OUT "", $lprimer, $ltemp, $lproblem,$rprimer, $rtemp, $rproblem, "\n"
		 
	}
	else
	{
		print OUT "", "", "", "", "", "", "", "\n";
	}
}
close OUT;

sub revdnacomp {
  my $dna = shift;
  my $revcomp = reverse($dna);

  $revcomp =~ tr/ACGTacgt/TGCAtgca/;

  return $revcomp;
}
