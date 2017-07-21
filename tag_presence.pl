#!/usr/bin/perl
use strict;
use warnings;

my $directory = $ARGV[0];
my $taglistfile = $ARGV[1];
my $checkseqlen = 40;
if ($ARGV[2]) 
{
	$checkseqlen = $ARGV[2];
}

$directory=~s/\/$//;
unless (-d $directory) 
{
	print "File directory $directory does not exist!";
	exit;
}

unless (-e $taglistfile) 
{
	print "Tag file $taglistfile does not exist!";
	exit;
}

my $minQ = 10;

my %seq2id = ();
my %matrix = ();
my %sampletotal = ();
my @samples = ();
my @ids = ();

$minQ = $minQ +33;
open OUT1, ">matrix.txt";
open OUT2, ">matrix_nomalized.txt";

open (IN, "$taglistfile") || die "Cannot open file $taglistfile";
<IN>;
LOOP:while (<IN>) 
{
	next LOOP unless (/\w/);
	chomp;
	my ($id, $seq, $leftp) = split "\t";	
	$id=~s/\s+//g;
	$seq=~s/\s+//g;
	$leftp=~s/\s+//g;

	
	if ($seq=~/($leftp.+)/) 
	{
		$seq = $1;
	}
	else
	{
		print "$id Tag_seq does not match left primer! Exit.\n";
		exit;
	}
    next LOOP if ((length $seq)< $checkseqlen);

	$seq=substr($seq, 0, $checkseqlen);
	$seq2id{$seq} = $id;
	push @ids, $id;
}

opendir(my $dh, $directory) || die "can't opendir $directory: $!";
my @files = grep { /\.gz$/ && -f "$directory/$_" } readdir($dh);
closedir $dh;

LOOP2:foreach my $file (@files) 
{
	print "Reading $file ...\n";
	my $sample = $file;
	$sample=~s/\.fastq\.gz//;
	push @samples, $sample;
	my $total=0;
	open ( IN,  "gunzip -c $directory/$file|") or die $!;
	LOOP3:while (<IN>) 
	{
		my $seq =<IN>;
		<IN>; 
		my $qual=<IN>;
		$seq=substr($seq, 0, $checkseqlen);
		$qual=substr($qual, 0, $checkseqlen);
		my @ascii = grep { $_ < $minQ } unpack("C*", $qual);
		if (@ascii>0) 
		{
			next LOOP3;
		}
		$total ++;

		my $id = $seq2id{substr($seq, 0, $checkseqlen)};
		if ($id) 
		{
			$matrix{$id}{$sample} ++;
		}
	}
	$sampletotal{$sample} = $total;
}

print OUT1 join "\t",  ("ID", @samples); print OUT1 "\n";
print OUT2 join "\t",  ("ID", @samples); print OUT2 "\n";
foreach my $id(@ids) 
{
	print OUT1 $id;
	print OUT2 $id;
	foreach  my $sample(@samples) 
	{
		my $c = $matrix{$id}{$sample};
		unless ($c) 
		{
			$c=0;
		}
		print OUT1 "\t", $c;

		my $total = $sampletotal{$sample};
		if ($total==0) 
		{		
			print OUT2 "\t", 0;
		}
		else
		{
			print OUT2 "\t",  sprintf("%.0f", 1000000 * $c/$total) ;
		}
		
	}
	print OUT1 "\n";
	print OUT2 "\n";
}
close OUT1;
close OUT2;
