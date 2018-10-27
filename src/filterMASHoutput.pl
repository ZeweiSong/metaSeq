#!/usr/bin/env perl
use strict;

my ($in,$cut,$out) = @ARGV;

open IN, "<$in" or die $!;
open OUT, ">$out" or die $!;
while(<IN>){
  my @line = split /\t/;
  if($line[0] ne $line[1] && $line[2] >= 0.02 && $line[2] <= $cut){
    my @L0 = split("/",$line[0]);
    my @L1 = split("/",$line[1]);
    print OUT "$L0[-1]\t$L1[-1]\t$line[2]\n";
  }
}
close IN;
close OUT;
exit;
