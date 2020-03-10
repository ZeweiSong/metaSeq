#!/usr/bin/env perl
# (c) 2019 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Annotate and stat reads mapped to references.
# Author:            Chao | fangchao@genomics.cn
# Version:           V0.1
# Last modified:    31 May 2019 (since 31 May 2019)
# ===================================================================
# see detail below
use strict;

unless (@ARGV){
  print "Usage: perl anno.pl <ID annotation file> > <blast6 or sam file>\n";
  exit;
}

my($anno,$align) = @ARGV;

open ANN, "<$anno" or die $!;
open ALN, "<$align" or die $!;

my %HASH;
if($anno =~ /silva.*.txt$/){
  while(<ANN>){
    chomp; my @a= split(/\t/,$_);
    $HASH{"$a[0].$a[1].$a[2]"} = "$a[4]\t$a[3]";
  }
}elsif($anno =~ /silva.*tax$/){
  while(<ANN>){
    chomp; my @a= split(/\t/,$_);
    $HASH{$a[0]} = "$a[1]\t$a[2]\t$a[3]\t$a[5]";
  }
}else{
  while(<ANN>){
    chomp; my @a= split(/\t/,$_);
    $HASH{$a[0]} = $_;
  }
}
close ANN;

if($align =~/.(blast6|m6)$/){
  while(<ALN>){
  	chomp;
  	my @a= split(/\t/,$_);
  	my $ann = $HASH{$a[1]};
  	print "$_\t$ann\n";
  }
}else{
  while(<ALN>){
    next if $_ =~ /^\@/;
  	chomp;
  	my @a= split(/\t/,$_);
  	my $ann = $HASH{$a[2]};
  	print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$ann\n";
  }
}

close ALN;
exit;
