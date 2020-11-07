#!/usr/bin/env perl
use warnings;
use strict;

die "Remove primers and adapters\nUsage: perl $0 [list] [in1] [in2] [o1] [o2]\n" if scalar @ARGV < 3;
our %dict= ('A'=>'T','T'=>'A','C'=>'G','G'=>'C');
my ($list,$i1,$i2,$o1,$o2) = @ARGV;

my ($fwd,$rev,$countF,$countR,$countU) = (0,0,0,0,0);
open L,"<$list" or die "cannot open $list.$!";
while(<L>){
  next if $_=~/^#/;
  chomp;
  my $len = length($_);
  for(my $i=0;$i<=$len-29;$i++){
    my $kmer = substr $_, $i, 29;
    $fwd .= ($fwd)?"|$kmer":"$kmer";
    $rev .= ($rev)?"|".&rev($kmer):&rev($kmer);
  }
}
print STDERR "Looking following seqs to remove:
    FWD: $fwd\n    REV: $rev\n";

open I1,"gzip -dc $i1|" or die $!;
open I2,"gzip -dc $i2|" or die $!;
open O1,"|gzip >$o1" or die $!;
open O2,"|gzip >$o2" or die $!;

while(<I1>){
  my $h1 = $_;   my $s1 = <I1>; <I1>; my $q1 = <I1>;
  my $h2 = <I2>; my $s2 = <I2>; <I2>; my $q2 = <I2>;
  my $hitF = ("$s1$s2" =~ /$fwd/g);
  my $hitR = ("$s1$s2" =~ /$rev/g);
  if($hitF + $hitR>0){
    my $debug =1 ;
  }
  if($hitF > 0){
    $countF ++;
  }elsif($hitR >0){
    $countR ++;
  }else{
    $countU ++;
    print O1 "$h1"."$s1+\n$q1";
    print O2 "$h2"."$s2+\n$q2";
  }
}

print STDERR "Found FWD: $countF ; REV: $countR; Written read-pairs: $countU\n";

sub rev{
  my $seq = shift;
  my $rev = "";
  for(my $i=length($seq)-1;$i>=0;$i--){
    $rev .= $dict{substr($seq,$i,1)};
  }
  return($rev)
}
