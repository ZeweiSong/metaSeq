#!/usr/bin/perl
use strict;

my($anno,$b6) = @ARGV;

open ANN, "<$anno" or die $!;
open B6, "<$b6" or die $!;

my %HASH;
while(<ANN>){
	chomp;
	my @a= split(/\t/,$_);
	$HASH{$a[1]} = $_;
}
close ANN;

while(<B6>){
	chomp;
	my @a= split(/\t/,$_);
	my $ann = $HASH{$a[1]};
	print "$_\t$ann\n";
}
close B6;
exit;
