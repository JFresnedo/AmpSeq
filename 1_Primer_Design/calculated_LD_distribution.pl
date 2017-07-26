#! /usr/bin/perl 
use strict;
#use Math::BigFloat; 
open IN1, $ARGV[0];
open IN2, $ARGV[1];
open OUT2, ">$ARGV[0]\_region.vcf";
my %maker=();
my $chr="";
my $pos="";
my $anchor="";
my @posset;
my $first="";
my $second="";
my $last="";
my @test="";
my $ref="";
my $homoAB=0;
my $heteraB=0;
my $heterAb=0;
my $homoab=0;
my $maxnumber=500;
my @jieceng=();
my $anchorchr;
my $row=0;
my $phase;
my $qq_input;
for (my $p=1;$p<=$maxnumber;$p++){
    $jieceng[0]=1;
    $jieceng[$p]=$jieceng[$p-1]+log($p);
   }


while (<IN1>){
   chomp;
   $row++;
    if ($row<3){
   $_=~s/\s+//g;
   ($chr,$pos)=split /\_/,$_;
   $anchorchr=$chr;
   $chr=~s/S//g;
   push @posset,$pos;
   $maker{$_}="";
    }
   else {
   $anchor=$_;
   $anchor=~s/\s+//g;
      }
      
}

my @vcf=<IN2>;
$first=shift @posset;
$second=shift @posset;
$last=pop @posset;
#print "$first\t$last\n";
for (my $i=0;$i<$#vcf;$i++){
     chomp $vcf[$i];
     my @vcfelements=split /\t/,$vcf[$i];
         if (exists $maker{$vcfelements[2]}){
           print OUT2 "$vcf[$i]\n";
          }
         elsif($vcfelements[0]==$chr && $vcfelements[1]>$first && $vcfelements[1]<$second){

             print OUT2 "$vcf[$i]\n";
    }

}
open INFILE, "$ARGV[0]\_region.vcf";
while (<INFILE>){
     chomp; 
        if ($_=~/.*$anchor.*/i){
       $ref=$_;
       #print "$ref\n";
   }
}
   open IN4, "$ARGV[0]\_region.vcf";
   open OUTSTATISTIC, ">$ARGV[0]\_statistic.out";
   $qq_input="$ARGV[0]\_statistic\.out";
  # print "$qq_input\n";
   print OUTSTATISTIC "Test\_SNP\tAnchor\_SNP\tP\_value\tAB\tAb\taB\tab\tPhase\n";
   @test=<IN4>;
      for (my $j=0;$j<=$#test;$j++){
          chomp $test[$j];
          my @elemtest=split /\t/, $test[$j];
           # print "$ref\n";
            my @refelms=split /\t/, $ref;
               for (my $k=9;$k<=$#elemtest;$k++){
                if ($elemtest[$k]=~/0\/0/ && $refelms[$k]=~/0\/0/){ $homoAB ++;}
                elsif ($elemtest[$k]=~/0\/1/ && $refelms[$k]=~/0\/0/) {$heteraB++; }
                elsif ($elemtest[$k]=~/1\/1/ && $refelms[$k]=~/0\/0/) {$heteraB++;}
                elsif ($elemtest[$k]=~/0\/0/ && $refelms[$k]=~/0\/1/) {$heterAb++;}
                elsif ($elemtest[$k]=~/0\/0/ && $refelms[$k]=~/1\/1/) {$heterAb++;}
                elsif ($elemtest[$k]=~/0\/1/ && $refelms[$k]=~/0\/1/) {$homoab++;}
                elsif ($elemtest[$k]=~/0\/1/ && $refelms[$k]=~/1\/1/) {$homoab++;}
                elsif ($elemtest[$k]=~/1\/1/ && $refelms[$k]=~/0\/1/) {$homoab++;}  
                elsif ($elemtest[$k]=~/1\/1/ && $refelms[$k]=~/1\/1/) {$homoab++;}
                else{} 

    }
    if (($homoAB > $heterAb) && ($homoab > $heteraB)){$phase=1;}
    elsif(($homoAB < $heterAb) && ($homoab < $heteraB)){$phase=2;}
    else{$phase=0;}
    my  $twotailed_value = &fisher($homoAB,$heterAb,$heteraB,$homoab);    
    print OUTSTATISTIC  "$elemtest[2]\t$refelms[2]\t$twotailed_value\t$homoAB\t$heterAb\t$heteraB\t$homoab\t$phase\n";
    $homoAB=0;$heteraB=0;$heterAb=0;$homoab=0; $phase=0;
}

sub fisher {
     my ($a,$b,$c,$d)=@_;
     my $iniatiala=$a;
     my $iniatialb=$b;
     my $iniatialc=$c;
     my $iniatiald=$d;
     my $rawP=&getP($a,$b,$c,$d); 
     #   my $min= (if ($c<$d) ? $c: $d);
         # for (my $k=0;$k<$min;$k++){
             # my $temp=
                 # }
       return $rawP;
     }
sub getrightP{
  my ($n1,$n2,$n3,$n4)=@_;
  my $j = 0;
  my $min=0;
  my $total=$n1+$n2+$n3+$n4;
  my $p1 += &getP($n1,$n2,$n3,$n4);
  if ($n3<$n2) {$min=$n3;}
  else {$min=$n2;}
  for ($j = 0; $j < $min; $j++){
   $p1 += &getP(++$n1,--$n2,--$n3,++$n4);
  }
  return $p1;

}    
sub getP{
 my ($v1,$v2,$v3,$v4)=@_;
 my $sum=$v1+$v2+$v3+$v4;
#  Math::BigFloat->precision(-50);
 my $p =(exp(($jieceng[$v1+$v2]+$jieceng[$v1+$v3]+$jieceng[$v3+$v4]+$jieceng[$v2+$v4])-($jieceng[$v1]+$jieceng[$v2]+$jieceng[$v3]+$jieceng[$v4]+$jieceng[$sum]))); 
 return $p;
    }
