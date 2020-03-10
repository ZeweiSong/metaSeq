#!/usr/bin/env perl
use warnings;
use strict;

die "Remove primers and adapters\nUsage: perl $0 [list] [in1] [in2] [o1] [o2]\n" if scalar @ARGV < 3;
our %dict= ('A'=>'T','T'=>'A','C'=>'G','G'=>'C');
my ($list,$i1,$i2,$o1,$o2) = @ARGV;

my $Ad = "";
open L,"<$list" or die "cannot open $list.$!";
while(<L>){
  next if $_=~/^#/;
  chomp;
  if($Ad eq ""){
    $Ad = "$_|".&rev($_);
  }else{
    $Ad .= "|$_|".&rev($_);
  }
}
print STDERR "Looking following seqs to remove: $Ad\n";

open I1,"<$i1" or die $!;
open I2,"<$i2" or die $!;
open O1,">$o1" or die $!;
open O2,">$o2" or die $!;

while(<I1>){
  my $h1 = $_;   my $s1 = <I1>; <I1>; my $q1 = <I1>;
  my $h2 = <I2>; my $s2 = <I2>; <I2>; my $q2 = <I2>;
  unless($s1=~/$Ad/||$s2=~/$Ad/){
    print O1 "$h1"."$s1+\n$q1";
    print O2 "$h2"."$s2+\n$q2";
  }
}

sub rev{
  my $seq = shift;
  my $rev = "";
  for(my $i=length($seq)-1;$i>=0;$i--){
    $rev .= $dict{substr($seq,$i,1)};
  }
  return($rev)
}
