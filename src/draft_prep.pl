#!/usr/bin/env perl
# (c) 2016 - 2019 Chao
# ===================================================================
# Desc:             Prepare draft mapping from SPAdes scaffolds
# Author:           fangchao@genomics.cn
# Dev cycle:        16 Dec 2018 - 16 Dec 2018
# Version:          Alpha_0.1
# ===================================================================
#
use Getopt::Long;

sub usage {
  my $msg = shift;
print <<USAGE;
$msg
usage:
  $0 -i spadesDir -m minLength -o outputDir
    -i  Directory containing SPAdes results
    -m  minimal length of sequence selected as index
    -o  output direcotry
    -v  verbose
    -h  show help info
USAGE
}
&usage && exit unless @ARGV;

my ($inDir,$minL,$outDir,$verbose,$help);
GetOptions(
  "i=s" => \$inDir,
  "m=i" => \$minL,
  "o=s" => \$outDir,
  "v" => \$verbose,
  "h|help|?" => \$help,
);
&usage && exit if $help;
&usage("[fatal] Essential input dir is missing") && exit unless defined $inDir;

open INF, "<$inDir/assembly_graph_with_scaffolds.gfa" or die $!;

open REF, ($outDir)?">$outDir/edges.fasta":">./edges.fasta" or die $!;
open LIK, ($outDir)?">$outDir/edges.fasta":">./edges.fasta" or die $!;
open NOD, ($outDir)?">$outDir/edges.fasta":">./edges.fasta" or die $!;

# Main start
&verbose("[log] Start ... ");
## module1: select edges with length > 1000bp
my(%SEQ,%LINK,%NODE,$log);
while(<INF>){
  chomp;
  my @cells = split("\t",$_);
  if($cells[0] eq "S"){
    $SEQ{$cells[1]} = $cells[2];
    my $len = length($cells[2]);
    if($len >= $minL){
      $log ++;
      print REF "\@$cells[1] length=$len\n$cells[2]\n"
    }
  }elsif($cells[0] eq "L"){
    
  }
}

close INF;
close REF;
# Main end

&verbose("$log sequences generated. done!\n");

exit;

sub verbose{
  my $msg = shift;
  print STDERR $msg if $verbose;
}
