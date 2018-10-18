#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;

die &usage if @ARGV < 2;
sub usage {
	print <<USAGE;
This script is to detect the position of bead's barcodes when you don't know it.
usage:
	perl $0 <barcode list file> <target sequence file>
USAGE
	exit;
}
my ($bf,$sf) = @ARGV;
my $baseFilename = basename($sf);
my (%BARCODE,%POS);
open BF, "< $bf" or die $!;
while(<BF>){
  chomp;
  my @read = split;
  $BARCODE{$read[0]} = $read[1];
}
close BF;

my $support = 0;
my $attempts= 0;
my $readLen;
open SF, openMethod($sf) or die $!;
while(<SF> && $support < 1000){
  chomp(my $id = $_);
  chomp(my $seq = <SF>);
  <SF>;<SF>;
  $readLen = length($seq);
  for(my $i=1;$i<=$readLen-9;$i++){
    my $kmer = substr($seq,$i-1,10);
    if(defined $BARCODE{$kmer}){
      $POS{$i} ++;
      $support = ($support>$POS{$i})?$support:$POS{$i};
    }
  }
	$attempts ++;
}
close SF;
print "Check file: $baseFilename. Read length range: [1,$readLen]. Tried $attempts times.\n";
my %POS2;
foreach my $code (keys %POS){
  $POS2{$POS{$code}} = $code;
}
my @sup =(sort {$b<=>$a}  keys %POS2);

printf("Top4 position(supports%%):  %3d (%.2f%%)\t%3d (%.2f%%)\t%3d (%.2f%%)\t*%3d (%.2f%%)\n",
$POS2{$sup[0]},$sup[0]/$attempts * 100,
$POS2{$sup[1]},$sup[1]/$attempts * 100,
$POS2{$sup[2]},$sup[2]/$attempts * 100,
$POS2{$sup[3]},$sup[3]/$attempts * 100);

exit;


###
sub openMethod {$_ = shift; return(($_ =~ /\.gz$/)?"pigz -dc $_|":"$_")}
