#!/usr/bin/env perl
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       stat the overlap depth between each 2 beads.
# Author:            Chao | fangchao@genomics.cn
# Version:           V0.1
# Last modified:    03 Jan 2019 (since 03 Jan 2019)
# ===================================================================
# see detail below
use strict;
use Getopt::Long;

sub usage {
  my $msg = shift;
print <<USAGE;
$msg
usage:
  $0 -i input -o output
    -i  unchime file with 10 column
    -m  mode [uc|tile]
    -o  output filename
    -v  verbose
    -h  show help info
USAGE
}
&usage("Stat the overlap depth between each 2 beads.") && exit unless @ARGV;

my ($inf,$mode,$ident,$match,$out,$verbose,$help);
GetOptions(
  "i=s" => \$inf,
  "m=s" => \$mode,
  "ident=i" => \$ident,
  "match=i" => \$match,
  "o=s" => \$out,
  "v" => \$verbose,
  "h|help|?" => \$help,
);
&usage && exit if $help;
#&usage("[fatal] Essential input is missing");

$ident||= 98;
$match||= 30;
#open INF, ($inf)?"<$inf":"<-" or die $!;
open OUT, ($out)?">$out":">-" or die $!;

# Main start
&run_uc    if $mode eq "uc";
&run_merge if $mode eq "merge";
# Main end

&verbose("[log] All done!\n");

exit;

sub verbose{
  my $msg = shift;
  print STDERR $msg if $verbose;
}

sub BC2Num {
  my $id = shift;
  my $num= shift;
  $num ||= '1';
  $id =~ /[@\/](\d+)_(\d+)_(\d+)\//;
  $id = "$1\_$2\_$3";
  my $status = ($id =~ /0000/)?"F":"T"; # FAIL or TRUE BARCODE
  return($id,$status);
}

sub run_uc {
  &verbose("[log] Mode [UC] start ... \n");
  my(%HASH,%CLUST,%STAT,$clusterCount,$linkCount);
  #Here is the shorts for each column. Find more detail in VSEARCH manual.
  ##   0         1       2            3            4  5  6      7     8        9##
  ##Type  #Cluster  length  similarity%  orientation  *  *  CIGAR query centroid##
  open INF, ($inf)?"<$inf":"<-" or die $!;
  while(<INF>){
    if($. % 1000000 == 0 ){
      &verbose("[log] lines: $. | links: $linkCount | cluster: $clusterCount\n");
    }
    $clusterCount ++ if $_ =~ /^S/;
    %CLUST = () && next if $_ =~ /^[CS]/; #Before the start of a new cluster
    chomp;
    my @info = split(/\t/,$_);
    my $iMatch = ($info[7] =~ /(\d+)M/)?$1:$info[2];
    $STAT{'O'} ++ && next if $iMatch < $match || $info[3] < $ident;  #Skip when identity or length less than cutoff;

    my ($BC0,$BC1,$BC8,$BC9,$s0,$s8,$s9) = (); #init
    unless (%CLUST){
      ($BC9,$s9) = &BC2Num($info[9]);
      $CLUST{$BC9} = 1 if $s9 eq "T";
    }
    ($BC8,$s8) = &BC2Num($info[8]);
    $STAT{'F'} ++ && next if $s8 eq "F";

    foreach my $BC9 (keys %CLUST){
      $STAT{'S'} ++ && next if $BC8 eq $BC9 || $BC9 eq "";
      $BC0 = ($BC8 lt $BC9)?$BC8:$BC9;
      $BC1 = ($BC8 lt $BC9)?$BC9:$BC8;
      $linkCount ++ unless defined $HASH{$BC0}{$BC1}{'count'};
      #$HASH{$BC0}{$BC1}{'status'} = ($BC0 eq $BC1)?"S":"T" unless defined $HASH{$BC0}{$BC1}{'status'};
      $HASH{$BC0}{$BC1}{'count'} ++;
      $HASH{$BC0}{$BC1}{'length'} += $info[2];
      $STAT{'T'} ++;
    }
    $CLUST{$BC8} = 1;
  }
  close INF;
  &verbose("[log] Finish read. stat ... \n");

  foreach my $BC0 (sort keys %HASH){
    foreach my $BC1 (sort keys %{$HASH{$BC0}}){
      my $status  = "T";#$HASH{$BC0}{$BC1}{'status'};
      my $count   = $HASH{$BC0}{$BC1}{'count'};
      my $length  = $HASH{$BC0}{$BC1}{'length'} / $count;
      my $output = sprintf("%s\t%s\t%s\t%5d\t%7.2f\n",
      $status, $BC0, $BC1, $count, $length);
      printf OUT $output;
    }
  }
  close OUT;
  &verbose(sprintf("[log] DUPLICATES SUMMARY:\n\nclusters\tedges\n%d\t%d\n",
  $clusterCount,$linkCount));
  &verbose(sprintf("[log] READS PAIR SUMMARY:\n\nFAIL\tSELF\tKEPT\n%d\t%d\t%d\n\n",
  $STAT{F},$STAT{S},$STAT{T}));
  &verbose("[log] Mode [UC] done ... \n");
}

sub run_merge {
  &verbose("[log] Mode [merge] start ... \n");
  my @files = split (",",$inf);
  &verbose("[err] [merge] needs two files. Pls check\n") & die $! if @files != 2;
  open IN1, "<$files[0]" or die $!;
  open IN2, "<$files[1]" or die $!;
  #Here is the shorts for each column.
  ##   0         1         2      3       4##
  ##Type  Barcode1  Barcode2  count  length##
  my ($run1, $run2, $count1, $count2) = (1,1,0,0);
  chomp(my $read1 = <IN1>);
  my @inf1 = split (/\t/,$read1);
  chomp(my $read2 = <IN2>);
  my @inf2 = split (/\t/,$read2);

  while($run1 && $run2){

    while($run1 && "$inf1[1] $inf1[2]" lt "$inf2[1] $inf2[2]"){
      $count1 ++;
      print OUT "$read1\n";
      chomp($read1 = <IN1>);
      $run1 = 0 && close IN1 if $read1 eq <EOF>;
      @inf1 = split (/\t/,$read1);
    }

    while("$inf1[1] $inf1[2]" eq "$inf2[1] $inf2[2]"){
      my $length = ($inf1[3] * $inf1[4] + $inf2[3] * $inf2[4]) / ($inf1[3] + $inf2[3]);
      my $output = sprintf("%s\t%s\t%s\t%5d\t%7.2f\n",
      $inf1[0], $inf1[1], $inf1[2], $inf1[3] + $inf2[3], $length);

      $count1 ++; $count2 ++;
      print OUT $output;

      chomp($read1 = <IN1>);
      @inf1 = split (/\t/,$read1);
      chomp($read2 = <IN2>);
      @inf2 = split (/\t/,$read2);
    }

    while($run2 && "$inf1[1] $inf1[2]" gt "$inf2[1] $inf2[2]"){
      $count2 ++;
      print OUT "$read2\n";
      chomp($read2 = <IN2>);
      $run2 = 0 && close IN2 if $read2 eq <EOF>;
      @inf2 = split (/\t/,$read2);
    }
    if($count1 % 1000000 == 0 || $count2 % 1000000 == 0){
      &verbose(sprintf("[log] Processing | file1: %9d | file2: %9d\n", $count1, $count2));
    }
  }
  &verbose("[log] file1 closed.\n") unless $run1;
  &verbose("[log] file2 closed.\n") unless $run2;
  while($run1){
    $count1 ++;
    print OUT "$read1\n";
    chomp($read1 = <IN1>);
    $run1 = 0 && close IN2 if $read1 eq <EOF>;
    @inf1 = split (/\t/,$read1);
  }
  while($run2){
    $count2 ++;
    print OUT "$read2\n";
    chomp($read2 = <IN2>);
    $run2 = 0 && close IN2 if $read2 eq <EOF>;
    @inf2 = split (/\t/,$read2);
  }

  close OUT;
  &verbose("[log] Mode [merge] done ... \n");
}
