#!/usr/bin/perl -w
use strict;
use warnings;

my $fastaFile = $ARGV[0];
my $blastFile = $ARGV[1];
my %seenInBlast;
open(my $HANDLE, "<", $blastFile) or die "Could not open $blastFile: $!\n";
while (!eof($HANDLE)) {
    my $line = readline($HANDLE);
    #print $line;

    chomp $line;
    my @splitLine = split("\t", $line);
    $seenInBlast{$splitLine[0]} = 1;    
}
close ($HANDLE);

open(my $FHANDLE, "<", $fastaFile) or die "Could not open $fastaFile $!\n";
my $printNextBool = 0;
while (!eof($FHANDLE)) {
    my $line = readline($FHANDLE);
    #print $line;
    chomp $line;
    if ($line =~ /^>.*/) {
        $line =~ s/>//;
        if (!exists($seenInBlast{$line})){
            $printNextBool = 1;
            print ">".$line."\n";
        }
    }elsif($printNextBool == 1){
        print $line."\n";
        $printNextBool = 0;
    }   
    
}
close ($FHANDLE);
