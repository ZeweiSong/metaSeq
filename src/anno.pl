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
use Getopt::Long;

unless (@ARGV){
  print "Usage: perl $0 -l level <ID annotation file> > <blast6 or sam file>\n";
  print "              -l  level to merge. default [0]:\n";
  print "                    6:species  5:genus  4:class  0:all\n";

  exit;
}

my ($lv);
GetOptions(
  "l=s" => \$lv,
);
$lv||=0;
my($anno,$align) = @ARGV;

open ANN, "<$anno" or die $!;
open ALN, "<$align" or die $!;

my %HASH;
if($anno =~ /silva.*.txt$/){
  while(<ANN>){
    chomp; my @a= split(/\t/,$_);
    $HASH{"$a[0].$a[1].$a[2]"}{A} = "$a[4]\t$a[3]";
    if($lv == 0){
      $HASH{"$a[0].$a[1].$a[2]"}{T} = $a[4];
    }else{
      my @ranks = split /;/, $a[3];
      $HASH{$a[0]}{T} = $ranks[$lv];
    }
  }
}elsif($anno =~ /.*tax$/){
  while(<ANN>){
    chomp; my @a= split(/\t/,$_);
    $HASH{$a[0]}{A} = "$a[1]\t$a[2]\t$a[3]\t$a[5]";
    if($lv == 0){
      $HASH{$a[0]}{T} = $a[2];
    }else{
      my @ranks = split /;/, $a[5];
      $HASH{$a[0]}{T} = $ranks[$lv];
    }
  }
}else{
  while(<ANN>){
    chomp; my @a= split(/\t/,$_);
    $HASH{$a[0]}{A} = $_;
    $HASH{$a[0]}{T} = $_;
  }
}
close ANN;

my (%TMPANN,$preQuery);
if($align =~/.(blast6|m6)$/){
  while(<ALN>){
  	chomp;
  	my @a= split(/\t/,$_);
    if($a[0] eq "Penicillium_expansum" && ($a[1] eq "AB028137.1.1770" or $HASH{$a[1]}{T} eq "Penicillium expansum")){
      my $debug = 1;
    }
    %TMPANN = () if $preQuery ne $a[0];
    next if exists $TMPANN{$HASH{$a[1]}{T}};
  	my $ann = $HASH{$a[1]}{A};
  	print "$_\t$ann\n";
    $TMPANN{$HASH{$a[1]}{T}} ++;
    $preQuery = $a[0];
  }
}else{
  while(<ALN>){
    next if $_ =~ /^\@/;
  	chomp;
  	my @a= split(/\t/,$_);
    %TMPANN = () if $preQuery ne $a[0];
    next if exists $TMPANN{$HASH{$a[2]}{T}};
  	my $ann = $HASH{$a[2]}{A};
  	print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$ann\n";
    $TMPANN{$HASH{$a[2]}{T}} ++;
    $preQuery = $a[0];
  }
}

close ALN;
exit;
