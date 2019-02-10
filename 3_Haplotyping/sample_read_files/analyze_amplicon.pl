#!/usr/local/bin/perl
use strict;
use warnings;
use Getopt::Std;
use diagnostics;
my %opts;

##v2
unless (getopts("s:k:i:o:l:b:c:r:a:w:v:t:m:d:x:j:h", \%opts))
{
        printhelp();
        print "Error: some options are not set properly!\n";
        exit;
}

if (defined $opts{"h"}) 
{
	 printhelp();
	 exit;
}

########################################
#base commands, default:AmpSeq_Container
########################################
our $trimmomatic_cmd = "java -jar /usr/local/bin/Trimmomatic-0.38/trimmomatic-0.38.jar";
our $trimmomatic_adapter = "/usr/local/bin/Trimmomatic-0.38/adapters/NexteraPE-PE.fa";
our $flash_cmd = "flash";
our $split_cmd = "fastx_barcode_splitter.pl";
our $collapse_cmd = "fastx_collapser -Q33 ";
our $blastdb_cmd = "makeblastdb";
our $blastn_cmd = "blastn";
our $bwa_cmd = "bwa";
our $samtools_cmd = "samtools";
our $bcf_cmd = "bcftools";
our $vcfutils_cmd = "vcfutils.pl"; ### part of the samtools package
our $clustalw_cmd = "clustalw2";
our $clustalo_cmd = "clustalo";
our $vcfconcatcmd = "vcf-concat"; ### part of vcftools
our $picard_path = "/usr/local/bin/picard/build/libs/picard.jar";
our $gatk_path = "/usr/local/bin/gatk/";


########################################
#set up parameters
########################################
my $threads = 8; 
my $jobs = 1;
if ($opts{"t"}) 
{
	$threads = $opts{"t"};
	unless ($threads=~/^\d+$/) 
	{
		printhelp();
		print "Error: -t must be a number!";
	}
}

if ($opts{"j"}) 
{
	$jobs = $opts{"t"};
	unless ($jobs=~/^\d+$/) 
	{
		printhelp();
		print "Error: -j must be a number!";
	}
}


#####trimmomatic parameter
my $quality_trimmer = 20;
my $minimum_length = 100;
if ($opts{"d"}) 
{
	my $t = $opts{"d"};
	if ($t=~/^(\d+):(\d+)$/)
	{
		$quality_trimmer = $1;
		$minimum_length = $2;
	}
	else
	{
		printhelp();
		print "Error: -d must be two numbers separated by ':'. E.g. 20:100. They are: quality score:minimum_read_length!\n";
		exit;
	}
}


###contig modes 1. contig; 2. concatenate; 3. mixed
my $contigmode =1;  ### contig the overlap reads
my $concatenate_read_length;
my $min_overlap_length = 20;
my $max_overlap_length = 250;
my $contiging_mismatch = 0.05;
##minimum and maximum overlap length when contiging
if ($opts{"v"}) 
{
	my $t = $opts{"v"};
	if ($t=~/^(\d+):(\d+):([0-9.]+)$/)
	{
		$min_overlap_length = $1;
		$max_overlap_length = $2;
		$contiging_mismatch = $3;
		$flash_cmd = "$flash_cmd -t $threads --allow-outies -m $min_overlap_length -M $max_overlap_length -x $contiging_mismatch ";
	}
	else
	{
		printhelp();
		print "Error: -v must be three numbers separated by ':'. E.g. 20:250:0.05.  The three numbers are: minimum overlap; maximum overlap; mismatch rate!\n";
		exit;
	}
}


if ($opts{"w"}) 
{
	$contigmode = $opts{"w"};
	unless ($contigmode=~/^[123]$/) 
	{
		printhelp();
		print "Error: the -w (contiging mode) can only be 1,2,3. You set the value as '$contigmode'\n";
		exit;
	}
	### specify the contiging parameter

	###specify the concatenate parameter
	if ($contigmode=~/^[23]$/)  
	{	
		if ($opts{"l"}) 
		{
			$concatenate_read_length = $opts{"l"};
			unless ($concatenate_read_length=~/^\d+$/) 
			{
				printhelp();
				print "Error: -l must be a number!";
				exit;
			}
		}
		else
		{
			printhelp();
			print "Error: You specified mode $contigmode for non-overlapping reads. You must specify '-l' truncated reads length\n"; 
			exit;
		}
	}
}


my $min_copies = 5; ##default minimum copy of reads per haplotype
if ($opts{"c"}) 
{
	$min_copies = $opts{"c"};
}

######genotyping step
my $haplotype_reads_ratio = 0.2;
if ($opts{"a"}) 
{
	$haplotype_reads_ratio = $opts{"a"};
}


########################################
#set up working directories and others
########################################

our %compnt = ("A"=>"T",
                        "T"=>"A",
                        "G"=>"C",
                        "C"=>"G",
                        "N"=>"N",
                        );


my $outdir = "out";
if ($opts{"o"}) 
{
	$outdir = $opts{"o"};
}
unless (-d $outdir) 
{
	mkdir $outdir;
}

my $statdir = "$outdir/stats";
unless (-d $statdir ) 
{
	mkdir $statdir ;
}


my $trimmdir = "$outdir/trimmed";
unless (-d $trimmdir ) 
{
	mkdir $trimmdir ;
}

my $assembledir = "$outdir/assemble";
unless (-d $assembledir ) 
{
	mkdir $assembledir ;
}

my $splitdir = "$outdir/split";
unless (-d $splitdir ) 
{
	mkdir $splitdir ;
}

my $collapsedir = "$outdir/collapse";
unless (-d $collapsedir ) 
{
	mkdir $collapsedir ;
}

my $mergedir = "$outdir/merged";
unless (-d $mergedir ) 
{
	mkdir $mergedir ;
}

my $tbtdir = "$outdir/haplotype2sample_raw";
unless (-d $tbtdir ) 
{
	mkdir $tbtdir ;
}

my $haplotype2fastadir = "$outdir/haplotype2fasta";
unless (-d $haplotype2fastadir ) 
{
	mkdir $haplotype2fastadir ;
}

my $samplegene2fastadir = "$outdir/tmp";
unless (-d $samplegene2fastadir ) 
{
	mkdir $samplegene2fastadir;
}

my $vcfdir = "$outdir/vcf";
unless (-d $vcfdir ) 
{
	mkdir $vcfdir;
}

open LOG, ">$outdir/runlog";

my $samplefile = "";
my $keyfile= "";
if ($opts{"s"}) 
{
	$samplefile = $opts{"s"};
}
else
{
	printhelp();
    print "Error: sample file is not specified!\n";
    exit;
}

if ($opts{"k"}) 
{
	$keyfile = $opts{"k"};
}
else
{
	printhelp();
    print "Error: key file is not specified!\n";
    exit;
}

my $hapfilter_aln_len = 0;
my $hapfilter_aln_len_pct = 0;
my $hapfilter_aln_pct = 0;
my $refseqfile;
my $refseqs;
if ($opts{"r"}) 
{
	$refseqfile = $opts{"r"};
	unless (-e $refseqfile)
	{
		print "Reference sequence (-r) is defined but cannot find the reference sequence $refseqfile\n";
		exit;
	}
	$refseqs = fasta2hash($refseqfile);
}

if ($opts{"x"}) 
{
	unless ($refseqfile)
	{
		print "The -x parameter requires the reference sequence file defined by -r. \n";
		exit;
	}
	if ($opts{"i"} && ($opts{"i"}=~/5/))
	{
		printhelp();
		print "Error: As -x is set. You cannot skip step 5. Remove 5 from -i and try again.\n";
		exit;	
	}
	if ($opts{"x"} =~/^(\d+):(\d+):(\d+)$/)
	{
		$hapfilter_aln_len = $1;
		$hapfilter_aln_len_pct = $2;
		$hapfilter_aln_pct = $3;
	}
	else
	{
		printhelp();
		print "Error: the parameter for -x is not right. It should by number:number like 25:90.\n";
		exit;	
	}
	open HAPCLEANLOG, ">$statdir/filter_hap.log";

	print HAPCLEANLOG "gene\tstatus\thap_id\tsequence(${hapfilter_aln_len}_bp)\tmatchlen\tpct\n";

}


my $runmsa = "";
if (defined $opts{"m"}) 
{
	my @t = split ":", $opts{"m"};
	if ($t[0]=~/clustalw/i) 
	{
		$runmsa = "$clustalw_cmd -INFILE=xxxxx -TYPE=DNA -OUTORDER=INPUT -OUTPUT=$t[1] -QUIET >& /dev/null";
	}
	elsif ($t[0]=~/clustalo/i) 
	{
		my $my_wrap = "";
		if ($t[2])
		{
			$my_wrap = " --wrap=$t[2] ";
		}
		$runmsa = "$clustalo_cmd --threads=$threads -i xxxxx --seqtype=DNA --output-order=input-order --outfmt=$t[1] -o xxxxx.alignment $my_wrap --force > /dev/null 2>&1 &";
	}
	else
	{
		printhelp();
		print "Error: the parameter for -m is not right. It should by clustalw:FORMAT or clustalo:FORMAT:length.\n";
		exit;
	}
}


open FASTQCLEANLOG, ">$statdir/filter_fastq.log";
print FASTQCLEANLOG "gene\tsample\tkept\trejected\n";

open (IN, "$samplefile")  || die "Error: cannot open sample file $samplefile\n";
my @samples;
my %sample2file1;
my %sample2file2;
my @genes;
my %gene2seq;

my $paired = 1;
LOOP1:while (<IN>) 
{
	next LOOP1 unless (/\w/);
	chomp;
	s/\s+$//;
	my ($sample, $file1, $file2) = split "\t";
	unless ($file1=~/\w/) 
	{
		next LOOP1;
	}


	if ((defined $sample) && ($sample=~/\w/) )
	{
		$sample=~s/\s//g;
		if (exists $sample2file1{$sample}) 
		{
			print "Error: Sample name $sample is duplicated in the file $samplefile. Exit now!\n";
			exit;
		}
		push @samples, $sample;
		if (-e $file1) 
		{
			$sample2file1{$sample} = $file1;
		}
		else
		{
			print "Error: $file1 in Sample file $samplefile does not exist! Exit now!\n";
			exit;
		}

		if ((defined $file2)&& ($file2=~/\w/) )
		{
			if (-e $file2) 
			{
				$sample2file2{$sample} = $file2;
			}
			else
			{
				print "$file2 in Sample file $samplefile does not exist! Exit now!\n";
				exit;
			}
		}
		else
		{
			$paired =0;
		}		 
	}
}
close IN;

if ($paired==0) 
{
	print "Error: single end sequencing data is not supported now! Make sure that in your sample file, there are two files for each sample!\n";
	exit;
}
open (IN, "$keyfile")  || die "Error: cannot open key file $keyfile\n";
while (<IN>)
{
	chomp;
	s/\s+$//;
	my ($key, $seq) = split "\t";


	if ((defined $key) && (defined $seq) && ($key=~/\w/) && ($seq=~/\w/)) 
	{
		$key=~s/\s//g;
		$seq=~s/\s//g;
		$seq=uc $seq;
		if ($seq=~/[^ACGT]/) 
		{
			print "In key file $keyfile, the gene $key sequence '$seq' has non 'ACGT' charactors! Exit now!\n";
			exit;
		}

		if ($key=~/(\S+)_([FR])$/) 
		{
			my $gene = $1;
			my $strand = $2;
			if (exists $gene2seq{$gene}{$strand}) 
			{
				print "Error: the key $key in $keyfile is duplicated!\n";
				exit;
			}
			$gene2seq{$gene}{$strand} = $seq;
			if ($strand eq "F") 
			{
				push @genes, $gene;
			}
		}
		else
		{
			print "Error: in keyfile $keyfile, the primer name $key is not correct!\n";
			exit;
		}
	}
}
close IN;

foreach my $gene (@genes) 
{
	unless ((exists $gene2seq{$gene}{"F"}) && (exists $gene2seq{$gene}{"F"})) 
	{
		print "Gene $gene miss either F or R primer\n"; 
		exit;
	}
}

####################################
## trim adapter sequences 
####################################
if (((! $opts{"i"})) || (($opts{"i"}) && (! ($opts{"i"}=~/1/))))
{
	print "Step 1. Trimming adapter and low quality reads.\n";
	print "  Trimming quality: $quality_trimmer; Minimum read length: $minimum_length\n\n";
	
	open STAT, ">$statdir/trim.stat.txt";
	LOOPi1:foreach my $sample (@samples) 
	{
		my $file1 = $sample2file1{$sample};
		my $file2 = $sample2file2{$sample};

		unless (-e $file1) 
		{
			print "Error: file for trimming $file1 does not exist!\n";
			exit;
		}

		unless (-e $file2) 
		{
			print "Error: file for trimming $file2 does not exist!\n";
			exit;
		}

		my $filesize = -s $file1;
		if ($filesize ==0 ) 
		{
			print "Warning: sample $sample file $file1 is empty; skipped!\n";
			##create empty files
			system ("touch $trimmdir/$sample.r1.fastq");
			system ("touch $trimmdir/$sample.r2.fastq");
			next LOOPi1;
		}

		my $trimmomatic_command = "$trimmomatic_cmd PE -phred33 $file1 $file2 $trimmdir/$sample.r1.fastq $trimmdir/$sample.u1.fastq $trimmdir/$sample.r2.fastq $trimmdir/$sample.u2.fastq ".
								"ILLUMINACLIP:${trimmomatic_adapter}:2:30:10 LEADING:$quality_trimmer TRAILING:$quality_trimmer SLIDINGWINDOW:4:15 MINLEN:$minimum_length";
		print "Trimming $sample ... ";
		print LOG $trimmomatic_command, "\n";
		my $trim_out = `$trimmomatic_command  2>&1`;
		 
		if ($trim_out =~/(Input Read Pairs.+)/)
		{
			$1=~/Input Read Pairs:\s+(\d+)\s+Both Surviving:\s+(\d+)\s+\(([0-9.]+)/;
			print  "Total $1; Survived: $2; Percent $3%\n";
			print  STAT "$sample\tTotal $1\tSurvived: $2\tPercent $3%\n";
			print LOG "$trim_out\n";
		}
		else
		{
			print LOG "$trim_out\n";
			close LOG;
			print "Error: Trimmomatic on $file1 $file2 failed!\n";
			print "It is also possible that Trimmomatic output format is changed! Please contat Qi Sun with the log file.\n";
			exit;
		}
		
	}
	print "Read trimming finished.\n\n";
	close STAT;
}
else
{
	print "##Step 1. Read trimming skipped!##\n\n";
}





####################################
## assemble read pair

####################################
if (((! $opts{"i"})) || (($opts{"i"}) && (! ($opts{"i"}=~/2/))))
{
	open STAT, ">$statdir/assemble.stat.txt";
########switch between concatenate or assemble the paired end, if concatenate, truncate read to defined length
	if ($contigmode eq "1") 
	{
		my $min_overlap_length = 20;
		print "Step 2. Contiging (mode 1) the overlapping read pairs. Discard the non-overlapping reads.\n";
		print "  Range of overlap: $min_overlap_length bp - $max_overlap_length bp; Mismatch rate: $contiging_mismatch\n";
	}
	elsif ($contigmode eq "2") 
	{
		print "Step 2. Concatenating (mode 2) the paired reads into one sequence. \n";
		print "  Truncating reads to $concatenate_read_length bp; Reverse complement the 2nd reads before concatenating.\n";	
	}
	else
	{
		print "Step 2. Mixed mode (mode 3) of contiging and concatenating. \n";
		print "  Range of overlap: $min_overlap_length bp - $max_overlap_length bp; Mismatch rate: $contiging_mismatch\n";
		print "  Concatenate reads with no overalop. Truncate reads to $concatenate_read_length bp; Reverse complement the 2nd reads before concatenating\n";	
	}

	LOOPi2:foreach my $sample (@samples) 
	{
		my $file1 = "$trimmdir/$sample.r1.fastq";
		my $file2 = "$trimmdir/$sample.r2.fastq";

		unless (-e $file1) 
		{
			print "Error: file for assemble $file1 does not exist!\n";
			exit;
		}

		unless (-e $file2) 
		{
			print "Error: file for assemble $file2 does not exist!\n";
			exit;
		}

		my $filesize = -s $file1;
		if ($filesize ==0 ) 
		{
			print "Warning: at assembly step: sample $sample file $file1 is empty; skipped!\n";
			##create empty files
			if ($contigmode=~/^[13]$/)
			{
				system ("touch $assembledir/$sample.extendedFrags.fastq");
			}
			if ($contigmode=~/^[23]$/)
			{
				system ("touch $assembledir/$sample.concat.fastq");
			}
			next LOOPi2;
		}
		if ($contigmode=~/^[13]$/) 
		{
			my $flash_command = " $flash_cmd -d $assembledir -o $sample $file1 $file2";
			print "Contiging read pair for $sample ...";
			print LOG $flash_command, "\n";
		
			my $flash_out = `$flash_command  2>&1`;
			print LOG $flash_out, "\n";
			if ($flash_out =~/Percent combined:\s*([0-9.]+)/)
			{
				my $pctassembled = $1;
				$flash_out =~/Total pairs:\s+(\d+)/;
				my $totalpairs = $1; 
				print " Total:  $totalpairs; Form contig: ${pctassembled}%\n";
				print STAT "$sample\tTotal: $totalpairs\tAssembled: ${pctassembled}%\n";
			}
			else
			{
				
				print "Flash on $file1 $file2 failed\n";
				close LOG;
				exit;
			}
		}
		if ($contigmode=~/^[23]$/)
		{
			################################
			#concatenate the reads, no assembly
			print "Concatenating reads in $sample ...\n";
			my $tofile = "$assembledir/$sample.concat.fastq";

			if ($contigmode=~/^[3]$/) 
			{
				$file1= "$assembledir/$sample.notCombined_1.fastq";
				$file2= "$assembledir/$sample.notCombined_2.fastq";
			}

			open (IN1, "$file1")|| die "Cannot read $file1\n";
			open (IN2, "$file2")|| die "Cannot read $file2\n";
			open (OUT_C, ">$tofile") || die "Cannot write to $tofile\n";
			
			my ($tln1, $ln1, $qln1, $ln2, $qln2);
			while ($tln1 = <IN1>) 
			{
				chomp $tln1;
				$ln1 = <IN1>;
				<IN1>;
				$qln1 = <IN1>;
				<IN2>;
				$ln2 = <IN2>;
				<IN2>;
				$qln2 = <IN2>;
				unless ($ln2) 
				{
					print "Error: $file1 and $file2 have unequal number of reads!\n";
					exit;
				}
				$ln1=~s/\s//g;
				$ln2=~s/\s//g;
				$qln1=~s/\s//g;
				$qln2=~s/\s//g;
				$ln1 = substr ($ln1.("N"x $concatenate_read_length), 0, $concatenate_read_length);
				$ln2 = substr ($ln2.("N"x $concatenate_read_length), 0, $concatenate_read_length);
				$qln1 = substr ($qln1.("!"x $concatenate_read_length), 0, $concatenate_read_length);
				$qln2 = substr ($qln2.("!"x $concatenate_read_length), 0, $concatenate_read_length);
				
				($ln2, $qln2) = revcom($ln2, $qln2);
				$ln1 = $ln1.$ln2;
				$qln1 = $qln1.$qln2;
				print OUT_C "$tln1\n$ln1\n+\n$qln1\n"; 
			}
			close OUT_C;
			###############################
		}

	}
	close STAT;
	print "Read assembly/concatenate finished.\n\n";
}
else
{
	print "##Step 2. Read assembly/concatenate skipped!##\n\n";
}


####################################
##  splitting file by genes 
####################################
if (((! $opts{"i"})) || (($opts{"i"}) && (! ($opts{"i"}=~/3/))))
{
	print "Step 3. Splitting sequences files by locus using primer F.\n";
	open STAT, ">$statdir/split.stat.txt";
	## user provided primer length to split gene
	my $barcode_length;
	if ($opts{"b"}) 
	{
		$barcode_length = $opts{"b"};
		unless ($barcode_length=~/^\d+$/) 
		{
			printhelp();
			print "Error: -b must be a number!";
		}
	}

	###make barcode file
	## get mininum length
	my $min_len = 99;
	my $max_len = 0;
	foreach my $gene (@genes) 
	{
		my $barcode = $gene2seq{$gene}{"F"};
		my $blen = length $barcode;
		if ($min_len > $blen) 
		{
			$min_len =$blen;
		}
		if ($max_len < $blen) 
		{
			$max_len =$blen;
		}
	}
	
	print "  Forward primer length:  $min_len to $max_len\n";
	if (defined $barcode_length) 
	{
		if ($barcode_length > $min_len) 
		{
			print STAT "The provided PCR primer barcode length $barcode_length longer than the shortest primer F length $min_len. The size $min_len is used.\n";
			print  "  The provided PCR primer barcode length $barcode_length longer than the shortest primer F length $min_len. The size $min_len is used.\n";
			$barcode_length = $min_len;
		}
		else
		{
			print STAT "Use $barcode_length nt to split genes in  sequence file\n";
			print "  Use $barcode_length bp to split locus in  sequence file\n";
		}
	}
	else
	{
		$barcode_length = $min_len;
		print STAT "Use smallest primer length  $min_len bp to split locus in  sequence file\n";
		print "  Use shortest primer length  $min_len bp to split locus in  sequence file\n";
	}



	open (B, ">$splitdir/barcode.txt") || die "Error: cannot make barcode file. Writing to $splitdir/barcode.txt failed\n";
	my %check_reduncancy = ();
	foreach my $gene (@genes) 
	{
		my $barcode = substr ($gene2seq{$gene}{"F"}, 0, $barcode_length);
		if (exists $check_reduncancy{$barcode}) 
		{
			print "Error: $barcode is duplicated, cannot uniquely identify a gene!\n";
			exit;
		}
		$check_reduncancy{$barcode} = "";
		print B "$gene\t$barcode\n", 
	}
	close B;

	my %genestat = ();
	LOOPi3:foreach my $sample (@samples) 
	{
		my $command_contig = "";
		my $command_concat = "";
		my $split_out = "";


		if ((-e "$assembledir/$sample.extendedFrags.fastq") && ((-s "$assembledir/$sample.extendedFrags.fastq")>0))
		{
			$command_contig = "$split_cmd --bcfile $splitdir/barcode.txt --prefix $splitdir/${sample}_tig_ --bol --mismatches 0 < $assembledir/$sample.extendedFrags.fastq";
			$split_out .= `$command_contig  2>&1`;
		}
		if ((-e "$assembledir/$sample.concat.fastq") && ((-s "$assembledir/$sample.concat.fastq")>0))
		{
			$command_concat = "$split_cmd --bcfile $splitdir/barcode.txt --prefix $splitdir/${sample}_cat_ --bol --mismatches 0 < $assembledir/$sample.concat.fastq";
			$split_out .= `$command_concat  2>&1`;
		}

		if ($split_out eq "") 
		{
			next LOOPi3;
		}


		unless (($command_contig=~/\w/) ||  ($command_concat=~/\w/))
		{
			print "Error: contigged/concatenated reads file from sample '$sample' does not exist.!\n";
			exit;;
		}

		
		if ($split_out =~/Count/)
		{
			my @total_counts = $split_out =~/total\s+(\d+)/g;
			my @unmatched_counts = $split_out =~/unmatched\s+(\d+)/g;	
			my $total_count = join ",", @total_counts;
			my $unmatched_count = join ",", @unmatched_counts;
			print "$sample\tTotal: $total_count\tUnmatched: $unmatched_count\n";
			print STAT "$sample\tTotal: $total_count\tUnmatched: $unmatched_count\n";
			print LOG "$split_out\n";
			
			open SPOUT, '<', \$split_out or die $!;
			<SPOUT>;
			while (<SPOUT>) 
			{
				chomp;
				my ($g, $c, $l) = split "\t";
				$genestat{$sample}{$g} .= $c. ",";
			}
			close SPOUT;

		}
		else
		{
			print LOG "$split_out\n";
			close LOG;
			print "Error: split reads by gene on sample $sample failed!\n";
			print "It is also possible that fastx_barcode_splitter.pl output format is changed! Please contact Qi Sun with the log file.\n";
			exit;
		}
	}
	close STAT;
	##making gstat file
	open GSTAT, ">$statdir/gene_sample_read_count.txt";
	print GSTAT "Gene\t";
	print GSTAT (join "\t", @samples);
	print GSTAT "\n";
	foreach my $gene (@genes) 
	{
		print GSTAT $gene;
		foreach  my $sample (@samples) 
		{
			my $count = $genestat{$sample}{$gene};
			unless ($count) 
			{
				$count =0;
			}
			$count=~s/,$//;
			print GSTAT "\t$count";
		}
		print GSTAT "\n";
	}
	close GSTAT;
	print "File splitting finished!\n\n";
}
else
{
	print "##Step 3. File splitting skipped!##\n\n";
}


####################################
## collapse identical sequences in each sample 
####################################
if (((! $opts{"i"})) || (($opts{"i"}) && (! ($opts{"i"}=~/4/))))
{
	print "Step 4. Collapsing identical reads in each sample\n";
	
	foreach my $gene (@genes) 
	{
		LOOPi4:foreach my $sample (@samples) 
		{
			print "Collapsing reads for sample $sample ...\n";
			system ("touch $collapsedir/${sample}_$gene");
			if ((-e "$splitdir/${sample}_tig_$gene") && ((-s "$splitdir/${sample}_tig_$gene")>0) )
			{
				system ("$collapse_cmd -i $splitdir/${sample}_tig_$gene -o $collapsedir/${sample}_$gene");

			}


			if((-e "$splitdir/${sample}_cat_$gene")&& ((-s "$splitdir/${sample}_cat_$gene")>0))
			{
				system ("$collapse_cmd -i $splitdir/${sample}_cat_$gene -o $collapsedir/tmpfile");
				open IN, "$collapsedir/tmpfile";
				open OUT, ">>$collapsedir/${sample}_$gene";
				while (<IN>)
				{
					if (/^>/) 
					{
						s/>/>C/;
					}
					print OUT $_;
				}
				close IN;
				close OUT;
			}
		}	
	}
	print "Collapsing identical reads in each sample finished.\n\n";
}
else
{
	print "##Step 4. Collapsing identical sequences skipped!##\n\n";
}


####################################
## collapse identical sequences for same locus across all samples
####################################
if (((! $opts{"i"})) || (($opts{"i"}) && (! ($opts{"i"}=~/5/))))
{
	print "Step 5. Collapsing reads across all samples.\n";
	
	LOOPi5:foreach my $gene (@genes) 
	{
		my @merging_cmd_tig = ();
		my @merging_cmd_cat = ();
		print "Collapsing gene $gene ...\n";
		LOOPi5:foreach my $sample (@samples) 
		{
			if (-e "$splitdir/${sample}_tig_$gene") 
			{
				push @merging_cmd_tig, "$splitdir/${sample}_tig_$gene";
			}

			if (-e "$splitdir/${sample}_cat_$gene") 
			{
				push @merging_cmd_cat, "$splitdir/${sample}_cat_$gene";
			}

		}

		if (@merging_cmd_tig>0) 
		{
			concatfile (\@merging_cmd_tig, "$mergedir/$gene.tig.fastq");
		}
		else
		{
			system ("touch $mergedir/$gene.tig.fastq");

		}

		if (@merging_cmd_cat>0) 
		{
			concatfile (\@merging_cmd_cat, "$mergedir/$gene.cat.fastq");
		}
		else
		{
			system ("touch $mergedir/$gene.cat.fastq");

		}


		system ("touch $mergedir/$gene");
		if ((-e "$mergedir/$gene.tig.fastq") && ((-s "$mergedir/$gene.tig.fastq")>0) )
		{
			system ("$collapse_cmd -i $mergedir/$gene.tig.fastq -o $mergedir/$gene");
		}

		
		if ((-e "$mergedir/$gene.cat.fastq") && ((-s "$mergedir/$gene.cat.fastq")>0) )
		{
			system ("$collapse_cmd -i $mergedir/$gene.cat.fastq -o $mergedir/tmpfile");
			open IN, "$mergedir/tmpfile";
			open OUT, ">$mergedir/$gene";
			while (<IN>)
			{
				if (/^>/) 
				{
					s/>/>C/;
				}
				print OUT $_;
			}
			close IN;
			close OUT;
		}

		if ($hapfilter_aln_len>0) 
		{
			if (defined $refseqs->{$gene}) 
			{		
				clean_hap_file($gene, "$mergedir/$gene", $refseqs->{$gene}, $hapfilter_aln_len, $hapfilter_aln_len_pct, $hapfilter_aln_pct);
			}
			else
			{
				print "Warning: $gene does not have a reference sequence. The haplotype for this gene is not cleaned.\n";
			}

		}
	}

	print "Collapsing identical sequences per locus finished.\n\n";
}
else
{
	print "##Step 5. Collapsing identical sequences for same locus across sample skipped!##\n\n";
}

####################################
## get haplotype by sample table for each locus
####################################
if (((! $opts{"i"})) || (($opts{"i"}) && (! ($opts{"i"}=~/6/))))
{
	print "Step 6. Making haplotype-by-sample matrix for each locus.\n";
	
	LOOPi6:foreach my $gene (@genes)
	{
		print "process gene $gene ... \n";
		my $all_haplotype_file = "$mergedir/$gene";
		unless (-e $all_haplotype_file) 
		{
			print "Warning: $gene does not have a file. Looking for $all_haplotype_file\n";
			next LOOPi6;
		}
		my $filesize = -s $all_haplotype_file;
		if ($filesize ==0) 
		{
			print "Warning: $gene file $all_haplotype_file is empty.\n";
			next LOOPi6;
		}
		open(INPUT, $all_haplotype_file) || die "ERROR: can't read input file: $!";
		open (OUT, ">$tbtdir/$gene.haplotypes.txt") || die "cannot write to $tbtdir/$gene.haplotypes.txt";
		my %seq2haplotype = ();
		my %haplotype2copies = ();
		my $haplotypeid = "";
		my $copies =0 ;
		my $line;
		my $seqstr = "";
		while($line=<INPUT>)
		{
			if($line=~/^>/)
			{
				if (($seqstr=~/\w/) && ($copies > $min_copies))
				{
						$seq2haplotype{$seqstr} = $haplotypeid;
						print OUT $haplotypeid, "\t", $copies, "\t", $seqstr, "\n";
						
				}
				$seqstr = "";
				$copies = 0;
				if ($line=~/(C?\d+)-(\d+)/)
				{

					$haplotypeid = $1;
					$copies = $2;
					if ($copies>$min_copies) 
					{
						$haplotype2copies{$haplotypeid} = $copies;
					}				
				}
				else
				{
					$haplotypeid = "unknown";
					$copies = 0;
				}
			}
			elsif ($line=~/\w/) 
			{
				if ($copies>$min_copies)
				{
					$line=~s/\s//g;
					$seqstr .= $line;
				}
				
			}		
        }
		if (($seqstr=~/\w/) && ($copies > $min_copies))
		{
			$seq2haplotype{$seqstr} = $haplotypeid;
			print OUT $haplotypeid, "\t", $copies, "\t", $seqstr, "\n";
						
		}
		close INPUT;
		close OUT;

		my @haplotypes = reverse sort {$haplotype2copies{$a}<=>$haplotype2copies{$b}} keys %haplotype2copies;

		

		###make tbt table
		my %tbt =();
		LOOPi6_2:foreach my $sample (@samples) 
		{
			my $samplegenefile = "$collapsedir/${sample}_$gene";


##############
			unless (-e $samplegenefile) 
			{
				print "Warning: $sample $gene does not have a file. Looking for $samplegenefile\n";
				next LOOPi6_2;
			}
			my $filesize = -s $samplegenefile;
			if ($filesize ==0) 
			{
				print "Warning: $samplegenefile file is empty.\n";
				next LOOPi6_2;
			}
			open(INPUT, $samplegenefile) || die "ERROR: can't read input file: $!";
			my $haplotypeid = "";
			my $copies;
			my $line;
			my $seqstr = "";
			LOOPi6_3:while($line=<INPUT>)
			{
				if($line=~/^>/)
				{
					if ($seqstr=~/\w/) 
					{
						if (exists $seq2haplotype{$seqstr}) 
						{
							$tbt{$seq2haplotype{$seqstr}}{$sample} = $copies;
						}
					}

					
					$seqstr = "";
					if ($line=~/(C?\d+)-(\d+)/)
					{
						$haplotypeid = $1;
						$copies = $2;
					}
					else
					{
						$haplotypeid = "unknown";
						$copies = 0;
					}

				}
				elsif ($line=~/\w/) 
				{
					$line=~s/\s//g;
					 $seqstr.= $line;
				}		
			}
			close INPUT;
			if ($seqstr=~/\w/) 
			{
				if (exists $seq2haplotype{$seqstr}) 
				{
					$tbt{$seq2haplotype{$seqstr}}{$sample} = $copies;
				}
			}
		}
		open (OUT, ">$tbtdir/$gene.haplotype2sample.txt") || die "cannot write to $tbtdir/$gene.haplotype2sample";
		print OUT "Haplotype\t";
		print OUT (join "\t", @samples);
		print OUT "\n";
		LOOPi6_3:foreach  my $haplotype (@haplotypes) 
		{
			next LOOPi6_3 unless (exists $haplotype2copies{$haplotype});
			print OUT $haplotype;
			foreach my $sample (@samples) 
			{
				my $count = 0;
				if (exists $tbt{$haplotype}{$sample}) 
				{
					$count = $tbt{$haplotype}{$sample};
				}
				print OUT "\t", $count;
			}
			print OUT "\n";
		}
		close OUT;


	}
	print "Making haplotype-by-sample table finished.\n\n";
}
else
{
	print "##Step 6. Making haplotype-by-sample matrix skipped!##\n\n";
}

####################################
## get genotypes
####################################
if (((! $opts{"i"})) || (($opts{"i"}) && (! ($opts{"i"}=~/7/))))
{
	my %sample2gt; 
	my %haplotypefrq;
	print "Step 7. Identifying top 2 haplotypes from each sample.\n";
	LOOPi7:foreach my $gene (@genes)
	{
		my $haplotype2samplefile = "$tbtdir/$gene.haplotype2sample.txt";
		my $haplotype2seqfile = "$tbtdir/$gene.haplotypes.txt";
		unless (-e $haplotype2samplefile) 
		{
			next LOOPi7;
		}

		my $filesize = -s $haplotype2samplefile;
		if  ($filesize == 0) 
		{
			next LOOPi7;
		}

		### read in the haplotype2sample file
		my %tbt = ();
		open (IN, "$haplotype2samplefile") || die "Cannot open $haplotype2samplefile \n";
		my $line = <IN>;
		chomp $line;
		my @thessamples = split "\t", $line;
		shift @thessamples;
		while (<IN>) 
		{
			chomp;
			my ($haplotype, @counts) = split "\t";
			for (my $i=0; $i<=$#counts; $i++) 
			{
				$tbt{$thessamples[$i]}{$haplotype} = $counts[$i];
			}
		}
		close IN;

		LOOPi7_1:foreach  my $sample(@samples) 
		{
			if (exists $tbt{$sample}) 
			{
				my %haplotype2counts = %{$tbt{$sample}};
				my @tophaplotypes = reverse sort {$haplotype2counts{$a} <=> $haplotype2counts{$b}} keys %haplotype2counts;
				my ($a1, $a2);
				if (@tophaplotypes>=2) 
				{
					($a1, $a2) = @tophaplotypes[0,1];
				}
				else
				{
					next LOOPi7_1;
				}
				
				if ($haplotype2counts{$a1}<=2) 
				{
					$sample2gt{$gene}{$sample} = "./.:0";
				}
				elsif ($haplotype2counts{$a2}==0) 
				{
					$sample2gt{$gene}{$sample} = "$a1/$a1:$haplotype2counts{$a1}";
					$haplotypefrq{$gene}{$a1} += 2;
				}
				elsif ($haplotype2counts{$a2}/$haplotype2counts{$a1}<$haplotype_reads_ratio) 
				{
					$sample2gt{$gene}{$sample} = "$a1/$a1:$haplotype2counts{$a1}";
					$haplotypefrq{$gene}{$a1} += 2;
				}
				else
				{
					$sample2gt{$gene}{$sample} = "$a1/$a2:$haplotype2counts{$a1},$haplotype2counts{$a2}";
					$haplotypefrq{$gene}{$a1} += 1;
					$haplotypefrq{$gene}{$a2} += 1;
				}
				#print "$gene\t$sample\t$sample2gt{$gene}{$sample}\n";

				
			}
			

		}
		my $fastafile = "$haplotype2fastadir/$gene.haplotypes.fa";

		open OUT, ">$fastafile";
		open IN, "$haplotype2seqfile";
		my %haplotype2seq = ();
		while (<IN>) 
		{
			chomp;
			my ($a, $c, $s) = split "\t";
			if (exists $haplotypefrq{$gene}{$a}) 
			{
				print OUT ">$a Reads:$c Samples: $haplotypefrq{$gene}{$a}\n$s\n";
				$haplotype2seq{$a} = $s;
			}
		}
		close IN;
		close OUT;
		#if (defined $opts{"m"}) 
		#{
		#	my $msaout = $opts{"m"};
		#	system ("$clustalw_cmd -INFILE=$fastafile -TYPE=DNA -OUTORDER=INPUT -OUTPUT=$msaout -QUIET >& /dev/null");
		#}
		if ($runmsa=~/\w/) 
		{
			my $t= $runmsa;
			$t=~s/xxxxx/$fastafile/g;
			system ($t); #print "\n";
		}
		#my $clustaout= `clustalw -INFILE=$fastafile -OUTFILE=$outfile 2>&1`;	
	}

	open OUT, ">$outdir/hap_genotype";
	print OUT "Locus\tHaplotypes\t";
	print OUT (join "\t", @samples);
	print OUT "\n";
	LOOPi7_2:foreach my $gene (@genes) 
	{
		unless ($haplotypefrq{$gene}) 
		{
			next LOOPi7_2;
		}
		my %frq = %{$haplotypefrq{$gene}};
		my @haplotypes = reverse sort {$frq{$a} <=> $frq{$b}} keys %frq;
		my $total=0;
		foreach my $haplotype (@haplotypes) 
		{
			$total += $frq{$haplotype};
		}
		my $haplotypestr = "";
		if ($total==0) 
		{
			next LOOPi7_2;
		}
		foreach my $haplotype (@haplotypes) 
		{
			my $frq =  sprintf("%.2f", $frq{$haplotype}/$total);
			$haplotypestr .= "$haplotype($frq);";
		}

		print OUT "$gene\t$haplotypestr";
		foreach my $sample (@samples) 
		{
			if (exists $sample2gt{$gene}{$sample}) 
			{
				print OUT "\t", $sample2gt{$gene}{$sample}; 
			}
			else
			{
				print OUT "\t", "./.:0"; 
			}
		}
		print OUT "\n";
	}
	close OUT;

	print "Identifying top 2 haplotypes finished.\n\n";
}
else
{
	print "##Step 7. Identifying top 2 haplotypes skipped!##\n\n";
}

####################################
## call SNPs
####################################
if (((! $opts{"i"})) || (($opts{"i"}) && (! ($opts{"i"}=~/8/))))
{
	print "Step 8. Calling SNPs.\n";
	##check if a filter_hap.log file exists
	my %rejected_haps = ();
	my $filter_fastq = 0;
	if (-e "$statdir/filter_hap.log") 
	{
		$filter_fastq =1;
		open FF, "$statdir/filter_hap.log";
		FFLOOP:while (<FF>) 
		{
			next FFLOOP unless (/\w/);
			my @data = split "\t";
			if ($data[1] eq "rejected")
			{
				push @{$rejected_haps{$data[0]}}, $data[3];
			}
		}
	}

	system("rm $samplegene2fastadir/*");
	unless ($refseqfile) 
	{
		print "Step 8 skipped: cannot perform SNP calling, as reference sequence is required for SNP calling.\n";
		exit;
	}


	my $tmpdir = "$outdir/tmp";
	unless (-d "$tmpdir") 
	{
		mkdir $tmpdir;
	}
	
	my $bwalog = "$samplegene2fastadir/bwa.log";
	system ("touch $bwalog");
	my $samtoollog = "$samplegene2fastadir/samtool.log";
	system ("touch $samtoollog");
	my $runvcf_concat = "$vcfconcatcmd ";
	LOOPi8:foreach my $gene (@genes)
	{
		unless (defined $refseqs->{$gene}) 
		{
			print "Warning: $gene does not have a reference sequence\n";
			next LOOPi8;
		}
		my $refseqstr = $refseqs->{$gene};
		open OUT, ">$samplegene2fastadir/ref_${gene}.fasta";
		print OUT ">$gene\n$refseqstr\n";
		close OUT;
		system ("bwa index $samplegene2fastadir/ref_${gene}.fasta &>> $bwalog");

		#my $pileupcommand = "$samtools_cmd mpileup -uf $samplegene2fastadir/ref_${gene} ";
		LOOPi82:foreach my $sample (@samples) 
		{
			my $rg = "\@RG\\tID:$sample\\tSM:$sample";
			my $fastqfile = "";
			if ((-e "$splitdir/${sample}_tig_$gene") && ((-s "$splitdir/${sample}_tig_$gene")>0) )
			{
				$fastqfile = "$splitdir/${sample}_tig_$gene";
			}
			elsif((-e "$splitdir/${sample}_cat_$gene")&& ((-s "$splitdir/${sample}_cat_$gene")>0))
			{
				$fastqfile = "$splitdir/${sample}_cat_$gene"; 
			}
			else
			{
				next LOOPi82;
			}
			#clean fastq file here
			if ($filter_fastq) 
			{
				if (exists $rejected_haps{$gene}) 
				{
					my @rejectes_haps = @{$rejected_haps{$gene}};

					my $newfastqfile = "$splitdir/tmp.fq";
					if ($fastqfile=~/_cat_/) 
					{
						$newfastqfile = "$splitdir/tmp_cat_.fq";
					}
					
					open FQ, "$fastqfile";
					open FQOUT, ">$newfastqfile";
					my ($line1, $line2, $line3, $line4);
					my $filterreads_count =0;
					my $keptreads_count =0;
					while ($line1=<FQ>) 
					{
						$line2=<FQ>;
						$line3=<FQ>;
						$line4=<FQ>;
						my $reject=0;
						CKLOOP:foreach  (@rejectes_haps) 
						{
							if ($line2=~/^$_/) 
							{
								$reject=1;
								last CKLOOP;
							}
						}
						if ($reject==0) 
						{
							print FQOUT "$line1$line2$line3$line4";
							$keptreads_count ++;
						}
						else
						{
							$filterreads_count ++;
						}

					}
					close FQ;
					close FQOUT;
					$fastqfile = $newfastqfile;
					print FASTQCLEANLOG "$gene\t$sample\t$keptreads_count\t$filterreads_count\n";
				}

			}
			
			#concatenate genes need to be separated into pair before running this command
			if($fastqfile =~/_cat_/) 
			{
				my $f1 = "$splitdir/tmp1.fq";
				my $f2 = "$splitdir/tmp2.fq";
				split_concatenated_fastq($fastqfile, $f1, $f2);
				system ("$bwa_cmd mem -R '\@RG\\tID:$sample\\tSM:$sample' $samplegene2fastadir/ref_${gene}.fasta $f1 $f2 1> $samplegene2fastadir/mygenesample.sam 2>> $bwalog");
			}
			else
			{
				system ("$bwa_cmd mem -R '\@RG\\tID:$sample\\tSM:$sample' $samplegene2fastadir/ref_${gene}.fasta $fastqfile 1> $samplegene2fastadir/mygenesample.sam 2>> $bwalog");
			}

			if (-e "$samplegene2fastadir/${gene}.sam") 
			{
				system ("java -jar ${picard_path} MergeSamFiles INPUT=$samplegene2fastadir/${gene}.sam INPUT=$samplegene2fastadir/mygenesample.sam SORT_ORDER=unsorted OUTPUT=$samplegene2fastadir/t.sam QUIET=true VERBOSITY=ERROR");
				system ("cp $samplegene2fastadir/t.sam $samplegene2fastadir/${gene}.sam");

			}
			else
			{
				system ("cp $samplegene2fastadir/mygenesample.sam $samplegene2fastadir/${gene}.sam");
			}
		}
		system ("java -jar ${picard_path} SamFormatConverter INPUT=$samplegene2fastadir/${gene}.sam OUTPUT=$samplegene2fastadir/${gene}.bam QUIET=true VERBOSITY=ERROR ");
		system ("java -jar ${picard_path} SortSam INPUT=$samplegene2fastadir/${gene}.bam OUTPUT=$samplegene2fastadir/${gene}.sorted.bam SORT_ORDER=coordinate QUIET=true VERBOSITY=ERROR");
		system ("java -jar ${picard_path} CleanSam INPUT=$samplegene2fastadir/${gene}.sorted.bam OUTPUT=$samplegene2fastadir/${gene}.clean.bam QUIET=true VERBOSITY=ERROR");
		system ("java -jar ${picard_path} BuildBamIndex INPUT=$samplegene2fastadir/${gene}.clean.bam QUIET=true VERBOSITY=ERROR");
		system ("$samtools_cmd faidx $samplegene2fastadir/ref_${gene}.fasta");
		system ("java -jar ${picard_path} CreateSequenceDictionary R=$samplegene2fastadir/ref_${gene}.fasta O=$samplegene2fastadir/ref_${gene}.dict QUIET=true VERBOSITY=ERROR");
		system ("${gatk_path}gatk HaplotypeCaller -R $samplegene2fastadir/ref_${gene}.fasta -O $vcfdir/${gene}.vcf -I $samplegene2fastadir/${gene}.clean.bam");
	}
	print "Calling SNP finished.!\n\n";
}
else
{
	print "##Step 8. Calling SNP skipped!##\n\n";
}

sub printhelp
{
	print "Usage: analyze_amplicon.pl -s mySampleFile -k myKeyFile\n";
	print "Options:\n";
	print "-o: outdir\n";
	print "-s: Sample file name\n";
	print "-k: Key file name\n";
	print "-i: skip steps. (steps are specified below.) \n";
	print "-d: parameter for trimming reads, using 2 numbers separated by ':'. Default: 20:100. They are: quality score:minimum_read_length\n";
	print "-w: running modes for paired end reads. 1(default).contig the overlapped reads, discard the reads that do not form contig; 2.concatenate the two reads using reads length defined by -l option, the reverse reads is reverse complemented; 3. mixed mode, contig the reads, and concate reads that do not contig  \n";
	print "-l: switch to concatenate reads mode. Truncate read to specified length then concatenate\n";
	print "-v: parameter for contiging, using 3 numbers separated by ':'. Default: 20:250:0.05. They are: minimum_overlap:max_overlap_length:mismatch_rate\n";
	print "-b: the length of F primer used for identify a locus in the sequence reads, default: the shortest primer length\n";
	print "-c: minimum copies for an haplotype across samples, default:2\n";
	print "-r: reference sequence\n";
	print "-x: filter haplotype by alignment to reference. e.g. -x 50:85:90, alignment with 50nt 5' sequence, 85% of the 50nt can be aligned, 90% identical in in aligned region.\n";
	print "-a: minimum ratio of read count between the two haplotypes  (low/high default: 0.2)\n";
	print "-m: Specify the multiple sequence alignment software, format, and  length per line, separated by : (e.g. clustalw:[CLUSTAL, GCG, GDE, PHYLIP, PIR, NEXUS, FASTA], or clustalo:[fasta, clustal,msf,phylip,selex,stockholm,vienna]:10000). Default: not generated.\n";
	print "-j: number of simultaneous jobs\n";
	print "-t: computing threads\n";
	print "-h: print help\n\n";
	print "Steps: \n";
	print "1. Trimming adapter and low quality reads.\n";
	print "2. Contiging/Concatnating read pairs.\n";
	print "3. Splitting sequences files by locus using primer F.\n";
	print "4. Collapsing identical reads in each sample.\n";
	print "5. Collapsing reads across all samples.\n";
	print "6. Making haplotype-by-sample matrix for each locus.\n";
	print "7. Identifying top 2 haplotypes from each sample.\n";
	print "8. Call SNPs.\n";
}

sub revcom
{
        my $seq = shift;
		my $qual = shift;
		$seq = uc $seq;
        my $seqout = "";
		my $qlout = "";
        my @nt = reverse ($seq=~/(\S)/g);
		my @ql = reverse ($qual=~/(\S)/g);
        foreach my $n (@nt)
        {
			my $q = shift @ql;

			if ($compnt{$n}) 
			{
				$seqout .=$compnt{$n};
				$qlout  .= $q;
			}
             
        }
        return ($seqout, $qlout);
}

sub fasta2hash
{
	my $file = shift;
	my %refseqs = ();
	open (IN, "$file") || die "cannot open $file";
	my $line;
	my $seqstr="";
	my $seqid = "";
	LOOPi8_1:while($line=<IN>)
	{
		if($line=~/^>(\S+)/)
		{
			my $newid = $1;
			if ($seqstr=~/\w/) 
			{
				$refseqs{$seqid} = $seqstr;
			}
			$seqstr = "";
			$seqid = $newid;
		}
		elsif ($line=~/\w/)  
		{
			$line=~s/\s//g;
			$seqstr.= $line;
		}
	}
	if ($seqstr=~/\w/) 
	{
		$refseqs{$seqid} = $seqstr;
	}
	return \%refseqs;
}

sub concatfile
{
    system qq( cat "$_" >> "$_[1]" ) for @{$_[0]};
}

sub clean_hap_file
{
	my ($gene, $mergehapfile, $refseqstr, $hapfilter_aln_len, $hapfilter_aln_len_pct, $hapfilter_aln_pct)= @_;
	## make blast query db and run blast
	open BQ, ">$mergedir/tmpquery";
	open HPIN, "$mergehapfile";
	while (<HPIN>) 
	{
		if (/^>/)
		{
			print BQ $_;
			$_ = <HPIN>;
			print BQ substr($_, 0, $hapfilter_aln_len), "\n";
		}
	}
	close HPIN;
	close BQ;

	open BQ, ">$mergedir/tmpdb";
	print BQ ">db\n$refseqstr\n";
	close BQ;

	system ("$blastdb_cmd -in $mergedir/tmpdb -dbtype nucl");
	system ("$blastn_cmd -num_threads $threads -db $mergedir/tmpdb -query $mergedir/tmpquery -out $mergedir/tmpblast -evalue 0.1  -task blastn-short -dust no -outfmt \"6 qseqid length pident\"");

	open BQ, "$mergedir/tmpblast";
	my %query_accepted;
	my %query2score;
	while (<BQ>) 
	{
		if (/\w/) 
		{
			my @t = split "\t";
			$query2score{$t[0]} = "$t[1]\t$t[2]";
			if (((100*$t[1]/$hapfilter_aln_len) >= $hapfilter_aln_len_pct)  && ($t[2] >= $hapfilter_aln_pct) )
			{
				$query_accepted{$t[0]} = "";
			}
		}
	}
	close BQ;

	open HPIN, "$mergehapfile";
	open BQ, ">$mergedir/tmpfasta";
	while (<HPIN>) 
	{
		if (/^>(\S+)/)
		{
			my $qid = $1;
			my $seqstr = <HPIN>;
			my $score = "no_hit\t";
			if (exists $query2score{$qid}) 
			{
				$score = $query2score{$qid};
			}
			if (exists $query_accepted{$qid}) 
			{
				print BQ ">$qid\n$seqstr\n";
				print HAPCLEANLOG $gene, "\taccepted\t", $qid, "\t", substr($seqstr, 0, $hapfilter_aln_len), "\t", $score, "\n";
			}
			else
			{
				print HAPCLEANLOG $gene, "\trejected\t", $qid, "\t", substr($seqstr, 0, $hapfilter_aln_len), "\t", $score, "\n";
			}
		}
	}
	close HPIN;
	close BQ;

	system ("mv $mergehapfile ${mergehapfile}.replaced");
	system ("mv $mergedir/tmpfasta $mergehapfile");
	
}

sub split_concatenated_fastq
{
	my $infile = shift;
	my $outfile1 = shift;
	my $outfile2 = shift;
	open SIN, $infile;
	open SOUT1, ">$outfile1";
	open SOUT2, ">$outfile2";

	my ($ln1, $ln2, $ln3, $ln4);
	my ($aln1, $aln2, $aln3, $aln4);
	my ($bln1, $bln2, $bln3, $bln4);

	while ($ln1 = <SIN>) 
	{
		$ln2 = <SIN>;
		$ln3 = <SIN>;
		$ln4 = <SIN>;
		chomp $ln1; chomp $ln2;chomp $ln3;chomp $ln4;
		my $len = length ($ln2)/2;
		unless ($len == int($len)) 
		{
			print "Warning: cannot split the sequence string when run split_concatenated_fastq\n";
			exit;
		}
		$aln2 = substr($ln2, 0, $len); $aln4 = substr($ln4, 0, $len);
		$bln2 = substr($ln2, $len, $len); $bln4 = substr($ln4, $len, $len);
		($bln2, $bln4) = revcom($bln2, $bln4);
		$aln1 = $ln1; $aln3 = $ln3;
		$bln1 = $ln1; $bln1=~s/ 1/ 2/; $bln3 = $ln3;
		print SOUT1 ">$aln1\n$aln2\n$aln3\n$aln4\n";
		print SOUT2 ">$bln1\n$bln2\n$bln3\n$bln4\n";
	}
	close SIN;
	close SOUT1;
	close SOUT2;
}

