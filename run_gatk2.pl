#!/usr/bin/perl
use strict;
use warnings;

my $datadir = $ARGV[0];
my $outdir = $ARGV[1];
my $grape_genome_db = $ARGV[2]; #bwa indexed database

$datadir=~s/\/$//;

our $picard_path = "/Programs/picard-tools/picard.jar";
our $fastx_trimmer_cmd = "gunzip -c xxxxx |fastx_trimmer -Q33 -f 1 -l 45 -o yyyyy.r1.fastq";
our $gatk_path = "/Programs/GATK/";

my $threads = 8;

unless (-d $datadir) 
{
	print "Data directory $datadir does not exist!\n";
	exit;
}

if (-d $outdir) 
{
	print "$outdir exists! Please delete before run\n";
	exit;
}
my $trimdir = "$outdir/trim";
mkdir $outdir;
mkdir $trimdir;

my $bwalog = "$outdir/bwalog";
my $trimlog = "$trimdir/trimlog";
open TRIM, ">$trimlog";
print TRIM "sample\ttotal\tsurvived\t%\n";

opendir(my $dh, $datadir)|| die "can't opendir $datadir: $!";
my @allfiles = grep { /fastq\.gz/i } readdir($dh);
closedir $dh;

LOOP1:foreach my $file (@allfiles)
{
	##parse sample name
	my $sample = $file;
	$sample=~s/\.fastq\.gz//;
	
	my $fastx_trimmer_command = $fastx_trimmer_cmd;
	$fastx_trimmer_command=~s/xxxxx/$datadir\/$file/;
	$fastx_trimmer_command=~s/yyyyy/$trimdir\/$sample/;
		print "Trimming $sample ... ";
		system ($fastx_trimmer_command);


		 
	my $bwa_cmd = "bwa mem -t $threads -R '\@RG\\tID:$sample\\tSM:$sample' $grape_genome_db $trimdir/$sample.r1.fastq  1> $outdir/${sample}.sam 2>> $bwalog";
	print "$bwa_cmd\n";
	system ($bwa_cmd);
	if (-e "$outdir/merged.sam") 
	{
		system ("java -jar ${picard_path} MergeSamFiles INPUT=$outdir/merged.sam INPUT=$outdir/${sample}.sam SORT_ORDER=unsorted OUTPUT=$outdir/t.sam QUIET=true VERBOSITY=ERROR");
		system ("cp $outdir/t.sam $outdir/merged.sam");

	}
	else
	{
		system ("cp $outdir/${sample}.sam $outdir/merged.sam");
	}


}
### merged sam created now run gatk
system ("java -jar ${picard_path} SamFormatConverter INPUT=$outdir/merged.sam OUTPUT=$outdir/merged.bam QUIET=true VERBOSITY=ERROR ");
system ("java -jar ${picard_path} SortSam INPUT=$outdir/merged.bam OUTPUT=$outdir/merged.sorted.bam SORT_ORDER=coordinate QUIET=true VERBOSITY=ERROR");
system ("java -jar ${picard_path} CleanSam INPUT=$outdir/merged.sorted.bam OUTPUT=$outdir/merged.clean.bam QUIET=true VERBOSITY=ERROR");
system ("java -jar ${picard_path} BuildBamIndex INPUT=$outdir/merged.clean.bam QUIET=true VERBOSITY=ERROR");
system ("java -jar ${gatk_path}GenomeAnalysisTK.jar -R $grape_genome_db -T UnifiedGenotyper -dt NONE -o $outdir/merged.vcf -I $outdir/merged.clean.bam");
