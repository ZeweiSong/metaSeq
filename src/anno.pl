#!/usr/bin/perl
use strict;

unless (@ARGV){
  print "Usage: perl anno.pl <ID annotation> <blast6 or sam>\n";
  exit;
}

my($anno,$align) = @ARGV;

open ANN, "<$anno" or die $!;
open ALN, "<$align" or die $!;

my %HASH;
while(<ANN>){
	chomp;
	my @a= split(/\t/,$_);
	$HASH{$a[1]} = $_;
}
close ANN;

if($align =~/.blast6$/){
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
