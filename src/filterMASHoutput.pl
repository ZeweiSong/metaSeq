#!/usr/bin/env perl
# Change barcode id and/or filter MASH distance by threshold.
#-----------------------------------------------------------------------------
# Author : Chao Fang
# Email  : fangchao@genomics.cn
# Create : Nov 2018
#-----------------------------------------------------------------------------
# see usage below
use strict;
use Getopt::Long;

sub usage {
  my $msg = shift;
print <<USAGE;
$msg
usage:
  $0 -i input -o output [-c cutoff] [-m id_map] [-M new_map]
    -i  distance file with three column (src dist weight)
    -o  output filename
    -c  if specified, only weight below cutoff will be kept.
    -m  if specified, id will changed by mapping in the id_map file
    -M  if specified, output id will renumbered and the mapping will sotored.
        Skipped when -m specified.
    -p  [b|p|u]
          b: bc format, default
          p: flag when distance contained pairwise seq id.
          u: universal format
    -v  verbose mode
    -h  show help info
USAGE
}
&usage && exit unless @ARGV;
my ($inf,$cut,$map,$MAP,$out,$pair,$verbose,$help);
GetOptions(
  "i=s" => \$inf,
  "o=s" => \$out,
  "c=f" => \$cut,
  "m=s" => \$map,
  "M=s" => \$MAP,
  "p=s" => \$pair,
  "v" => \$verbose,
  "h|help|?" => \$help,
);
&usage && exit if $help;
#&usage("[fatal] Essential input is missing") && exit unless defined $inf;

my($nid,%MAP)=(0,());
$pair ||= "b";
open INF, ($inf)?"<$inf":"-" or die $!;
open OUT, ($out)?">$out":">STDOUT" or die $!;

if($map){
  &verbose("[log] map enabled. Reading $map ... ");
  open MAP,"<$map" or die $!;
  while(<MAP>){
    chomp;
    my @L = split(/\t/,$_);
    $MAP{$L[0]} = $L[2];
  }
  close MAP;
  &verbose("done!\n");
}elsif($MAP){
  &verbose("[log] new_map enabled. The mapping will be writen in $MAP\n");
  open MAP,">$MAP" or die $!;
}
if($cut){
  &verbose("[log] filter threshold : $cut \n");
}else{
  $cut = 1;
  &verbose("[log] filter threshold undefined. Set to max ($cut)\n");
}
&verbose("[log] Start Writing ... ");
while(<INF>){
  my @line = split /\t/;
  if($line[0] ne $line[1] && $line[2] > 0 && $line[2] <= $cut){
    my($src,$dist,$s_p,$d_p) = ();
    if($pair eq "p"){
      $line[0] =~ /(\d\d\d\d_\d\d\d\d_\d\d\d\d.\d)/;
      $src = $1;
      $line[1] =~ /(\d\d\d\d_\d\d\d\d_\d\d\d\d.\d)/;
      $dist = $1;
    }elsif($pair eq "b"){
      $line[0] =~ /(\d\d\d\d_\d\d\d\d_\d\d\d\d)/;
      $src = $1;
      $line[1] =~ /(\d\d\d\d_\d\d\d\d_\d\d\d\d)/;
      $dist = $1;
    }elsif($pair eq "u"){
      $src = $line[0];
      $dist = $line[1];
    }

    if($map){
      print OUT "$MAP{$src}\t$MAP{$dist}\t$line[2]\n";
    }elsif($MAP){
      &generateMAP($src);
      &generateMAP($dist);
      print OUT "$MAP{$src}\t$MAP{$dist}\t$line[2]\n";
    }else{
      print OUT "$src\t$dist\t$line[2]\n";
    }
  }
  if($. % 1000000 == 0){
    &verbose(sprintf("[log] Processing %10d\n", $.));
  }
}
close INF;
close OUT;
close MAP if $MAP;

&verbose("done!\n");

exit;

sub verbose{
  my $msg = shift;
  print STDERR $msg if $verbose;
}

sub generateMAP{
  my $id = shift;
  unless(defined $MAP{$id}){
    $MAP{$id} = $nid;
    print MAP "$nid\t$id\n";
    $nid ++;
  }
}
