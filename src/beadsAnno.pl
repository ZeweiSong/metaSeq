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

sub usage {
  my $msg = shift;
print <<USAGE;
$msg
usage:
  $0 -i input -o output
    -m  mode [uc|tile|bc|bin|otu2b|f2b|cbr]
    -i  unchime file with 10 column
    -l  rank (for f2b)
    -o  output filename
    -b  output beadfile name
    -v  verbose
    -h  show help info
USAGE
}
&usage("Stat the overlap depth between each 2 beads.") && exit unless @ARGV;

our ($inf,$rank,$mode,$ident,$match,$out,$obd,$verbose,$debug,$help);
GetOptions(
  "i=s" => \$inf,
  "m=s" => \$mode,
  "l=s" => \$rank,
  "ident=i" => \$ident,
  "match=i" => \$match,
  "o=s" => \$out,
  "b=s" => \$obd,
  "v" => \$verbose,
  "d" => \$debug,
  "h|help|?" => \$help,
);
&usage && exit if $help;
#&usage("[fatal] Essential input is missing");

$ident||= 98;
$match||= 30;
$rank||=7;
#open INF, ($inf)?"<$inf":"<-" or die $!;
open OUT, ($out)?">$out":">-" or die $!;

# Main start
my (%CLUST);
&run_uc    if $mode eq "uc";
&run_merge if $mode eq "merge";
&run_BCent if $mode eq "bc";
&run_otu2b if $mode eq "otu2b";
&run_F2b   if $mode eq "f2b";
&run_B2b   if $mode eq "b2b";
&run_cbr   if $mode eq "cbr"; #check beads region
&run_otupair if $mode eq "otupair";
&run_otuGC if $mode eq "gc";
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

sub run_BCent {
  &verbose("[log] Mode [Bead Centroid] start ... \n");
  my (%STAT,$RCCount,$BCCount);
  #Here is the shorts for each column. Find more detail in VSEARCH manual.
  ##   0         1       2            3            4  5  6      7     8        9##
  ##Type  #Cluster  length  similarity%  orientation  *  *  CIGAR query centroid##
  open INF, ($inf)?"<$inf":"<-" or die $!;
  $RCCount = 0;
  $BCCount = 1;
  while(<INF>){
    next if $_ =~ /^C/;
    if($. % 1000000 == 0 ){
      &verbose("[log] lines: $. | RC: $RCCount | BC: $BCCount\n");
    }
    if($_ =~ /^S/){
      &writePreCluster;
      %CLUST = () ;
      $RCCount ++ ;
    }
    chomp;
    my @info = split(/\t/,$_);
    #my $iMatch = ($info[7] =~ /(\d+)M/)?$1:$info[2];
    #$STAT{'O'} ++ && next if $iMatch < $match || $info[3] < $ident;  #Skip when identity or length less than cutoff;

    my ($BC0,$BC1,$BC8,$BC9,$s0,$s8,$s9) = (); #init
    ($BC8,$s8) = &BC2Num($info[8]);
    if($CLUST{$BC8}{'S'}){
      $CLUST{$BC8}{'H'}{$CLUST{$BC8}{'C'}} = sprintf("H\t%8d\t%3d\t%5s\t%s\t%s\t%s\t%s\t%s\t%s\n",
      $BCCount-1,$info[2],$info[3],$info[4],$info[5],$info[6],$info[7],$info[8],$CLUST{$BC8}{'Center'});
      $CLUST{$BC8}{'C'} ++;
    }else{
      $CLUST{$BC8}{'S'} = sprintf("S\t%8d\t%3d\t%5s\t%s\t%s\t%s\t%s\t%s",
      $BCCount-1,$info[2],$info[3],$info[4],$info[5],$info[6],$info[7],$info[8]);

      $CLUST{$BC8}{'C'} = 1;
      $CLUST{$BC8}{'Center'} = $info[8];
      $BCCount ++;
    }
  }

  close INF;
  close OUT;

  &verbose("[log] Mode [Bead Centroid] done ... \n");
}

sub writePreCluster{
  return(1) unless %CLUST;
  foreach my $BC (sort keys %CLUST){
    print OUT "$CLUST{$BC}{'S'}\t$CLUST{$BC}{'C'}\n";
    print "$BC\t$CLUST{$BC}{'C'}\n";
    foreach my $HT (sort {$a<=>$b} keys %{$CLUST{$BC}{'H'}}){
      print OUT "$CLUST{$BC}{'H'}{$HT}";
    }
  }
}


################################################################################
# summary SSU and LSU annotation for each otu
################################################################################
sub usage4otu2b {
  my $msg = shift;
  print <<USAGE;
$msg
usage [v0.1] :
  $0 otu2b -L <LSU annotation file> -S <SSU annotation file> -o output
    -L  LSU annotation file
    -S  SSU annotation file
    -o  output
    -v  verbose
    -h  show help info
USAGE
}

sub run_otu2b {
  &verbose("[log] Mode [O2B] start ... \n");
  my @files = split (",",$inf);
  &verbose("[err] [O2B] needs two files. Pls check\n") & die $! if @files != 2;
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
  &verbose("[log] Mode [O2B] done ... \n");
}

################################################################################
# For bacteria, summary silva SSU and LSU annotation for each otu
################################################################################
sub usage4B2b {
  my $msg = shift;
  print <<USAGE;
$msg
usage [v0.1] :
  $0 b2b -m f2b -i SSU,LSU -o output
    -i  LSU,SSU annotation file
    -l  rank
    -o  output
    -v  verbose
    -h  show help info
USAGE
}
sub run_B2b {
  &verbose("[log] Mode [B2B] start ... \n");
  my @files = split (",",$inf);
  &verbose("[err] [B2B] needs three files. Pls check\n") & die $! if @files != 2;
  open my $IN1, "<$files[0]" or die $!;
  open my $IN2, "<$files[1]" or die $!;

  my ($count,%LSU,%SSU,%ITS,%PLSU,%PSSU,%PITS) = (0,);
  my %EMPTY = ("MAIN" => {"NONE" => {"A" => [(0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0)]}});
  chomp(my $read1 = <$IN1>);
  my @inf1 = split (/\t/,$read1);
  $inf1[0] =~ /(\d+)$/; my $no1 = $1;
  $LSU{NEXT}{A} = \@inf1; $LSU{NEXT}{N} = $no1; $LSU{FILEHANDLE} = $IN1;
  chomp(my $read2 = <$IN2>);
  my @inf2 = split (/\t/,$read2);
  $inf2[0] =~ /(\d+)$/; my $no2 = $1;
  $SSU{NEXT}{A} = \@inf2; $SSU{NEXT}{N} = $no2; $SSU{FILEHANDLE} = $IN2;

  my @pick = &checkOrder($no1,$no2,<EOF>);
  while($pick[1] || $pick[2] ){
    # read LSU annotations of 1 OTU if allowed
    if($pick[1]){
      %LSU  = &read1otu(\%LSU);
      %PLSU = %LSU;
      $no1  = $LSU{NEXT}{N};
    }else{
      %PLSU = %EMPTY;
    }
    # read SSU annotations of 1 OTU if allowed
    if($pick[2]){
      %SSU  = &read1otu(\%SSU);
      %PSSU = %SSU;
      $no2  = $SSU{NEXT}{N};
    }else{
      %PSSU = %EMPTY;
    }

    #summary annotation of this OTU
    &summaryAnno(\%PLSU,\%PSSU,\%EMPTY);
    $count ++;
    #verbose
    if($count % 1000 == 0 ){
      &verbose(sprintf("[log] Processing %6d th LOTU\r",$count));
    }
    #prepare for next loop
    @pick = &checkOrder($no1,$no2,<EOF>);
    %EMPTY = ("MAIN" => {"NONE" => {"A" => [(0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0)]}});
  }
  &verbose(sprintf("[log] Processing %6d th OTU [done]\n",$count));
  close OUT;
  &verbose("[log] Mode [B2B] done ... \n");
}


################################################################################
# For fungi, summary SSU, LSU and UNITE annotation for each otu
################################################################################
sub usage4F2b {
  my $msg = shift;
  print <<USAGE;
$msg
usage [v0.1] :
  $0 f2b -m f2b -i SSU,LSU,UNITE -o output
    -i  LSU,SSU,ITS annotation file
    -l  rank
    -U  UNITE annotation file
    -o  output
    -v  verbose
    -h  show help info
USAGE
}

sub run_F2b {
  &verbose("[log] Mode [F2B] start ... \n");
  my @files = split (",",$inf);
  &verbose("[err] [F2B] needs three files. Pls check\n") & die $! if @files != 3;
  open my $IN1, "<$files[0]" or die $!;
  open my $IN2, "<$files[1]" or die $!;
  open my $IN3, "<$files[2]" or die $!;

  my ($count,%LSU,%SSU,%ITS,%PLSU,%PSSU,%PITS) = (0,);
  my %EMPTY = ("MAIN" => {"NONE" => {"A" => [(0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0)]}});
  chomp(my $read1 = <$IN1>);
  my @inf1 = split (/\t/,$read1);
  $inf1[0] =~ /(\d+)$/; my $no1 = $1;
  $LSU{NEXT}{A} = \@inf1; $LSU{NEXT}{N} = $no1; $LSU{FILEHANDLE} = $IN1;
  chomp(my $read2 = <$IN2>);
  my @inf2 = split (/\t/,$read2);
  $inf2[0] =~ /(\d+)$/; my $no2 = $1;
  $SSU{NEXT}{A} = \@inf2; $SSU{NEXT}{N} = $no2; $SSU{FILEHANDLE} = $IN2;
  chomp(my $read3 = <$IN3>);
  my @inf3 = split (/\t/,$read3);
  $inf3[0] =~ /(\d+)$/; my $no3 = $1;
  $ITS{NEXT}{A} = \@inf3; $ITS{NEXT}{N} = $no3; $ITS{FILEHANDLE} = $IN3;

  my @pick = &checkOrder($no1,$no2,$no3);
  while($pick[1] || $pick[2] || $pick[3]){
    # read LSU annotations of 1 OTU if allowed
    if($pick[1]){
      %LSU  = &read1otu(\%LSU);
      %PLSU = %LSU;
      $no1  = $LSU{NEXT}{N};
    }else{
      %PLSU = %EMPTY;
    }
    # read SSU annotations of 1 OTU if allowed
    if($pick[2]){
      %SSU  = &read1otu(\%SSU);
      %PSSU = %SSU;
      $no2  = $SSU{NEXT}{N};
    }else{
      %PSSU = %EMPTY;
    }
    # read ITS annotations of 1 OTU if allowed
    if($pick[3]){
      %ITS  = &read1otu(\%ITS);
      %PITS = %ITS;
      $no3  = $ITS{NEXT}{N};
    }else{
      %PITS = %EMPTY;
    }
    #summary annotation of this OTU
    &summaryAnno(\%PLSU,\%PSSU,\%PITS);
    $count ++;
    #verbose
    if($count % 1000 == 0 ){
      &verbose(sprintf("[log] Processing %6d th LOTU\r",$count));
    }
    #prepare for next loop
    @pick = &checkOrder($no1,$no2,$no3);
    %EMPTY = ("MAIN" => {"NONE" => {"A" => [(0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0)]}});
  }
  &verbose(sprintf("[log] Processing %6d th OTU [done]\n",$count));
  close OUT;
  &verbose("[log] Mode [F2B] done ... \n");
}

sub checkOrder{
  my ($id1,$id2,$id3) =@_;
  $id1=($id1 eq <EOF>)?"~":$id1;
  $id2=($id2 eq <EOF>)?"~":$id2;
  $id3=($id3 eq <EOF>)?"~":$id3;
  my %HS = ( 1 => $id1, 2 => $id2, 3 => $id3);
  my @od = sort { $HS{$a} cmp $HS{$b} } keys %HS;
  my @pick = (0,0,0,0,$HS{$od[0]});
    # $pick[0] => number of picked files to read
    # $pick[1] => pick file 1 or not
    # $pick[2] => pick file 2 or not
    # $pick[3] => pick file 3 or not
    # $pick[4] => order of the first file
  if($HS{$od[0]} eq "~"){
    #Do nothing;
  }elsif($HS{$od[0]} lt $HS{$od[1]}){
    $pick[0] =1;
    $pick[$od[0]] =1;
  }elsif($HS{$od[1]} lt $HS{$od[2]}){
    $pick[0] =2;
    $pick[$od[0]] =1;
    $pick[$od[1]] =1;
  }else{
    $pick[0] =3;
    $pick[$od[0]] =1;
    $pick[$od[1]] =1;
    $pick[$od[2]] =1;
  }
  return(@pick);
}

sub read1otu {
  my $HS = shift;
  #load first line:
  my ($no1,$IN,@inf1) = ($$HS{NEXT}{N},$$HS{FILEHANDLE},@{$$HS{NEXT}{A}});
  #rest hash
  %$HS = ();
  my ($no,@inf,$order) = ($no1,@inf1,0);
  while($no eq $no1){
    #treat previous record
    my @tags = (exists $$HS{MAIN})?("ALL"):("ALL","MAIN");
    my @ranks = split /;/, $inf[17];
    $order ++;
    foreach my $tag (@tags){
      $$HS{$tag}{$inf[15]}{RefID}    ||= $inf[1];
      $$HS{$tag}{$inf[15]}{taxID}    ||= $inf[15];
      $$HS{$tag}{$inf[15]}{taxName}  ||= $inf[13];
      $$HS{$tag}{$inf[15]}{bitScore} ||= $inf[11];
      $$HS{$tag}{$inf[15]}{identity} ||= $inf[2];
      $$HS{$tag}{$inf[15]}{rank}     ||= $ranks[$rank];
      $$HS{$tag}{$inf[15]}{rank_1}   ||= ($ranks[$rank-1])?$ranks[$rank-1]:$ranks[$rank];
      $$HS{$tag}{$inf[15]}{order}    ||= $order;
      $$HS{$tag}{$inf[15]}{A}        ||= [@inf];
    }

    #read new record
    chomp(my $read = <$IN>);
    if($read eq <EOF>){
      @inf = (<EOF>);
      $no = <EOF>;
      close $IN;
      last;
    }else{
      @inf = split (/\t/,$read);
      $inf[0] =~ /(\d+)$/; $no = $1;
    }

  }
  $$HS{NEXT}{A} = [@inf];
  $$HS{NEXT}{N} = $no;
  $$HS{FILEHANDLE} = $IN;
  return(%$HS);
}

sub summaryAnno{
  my ($L,$S,$I) = @_;
  my %MHASH;
  my ($kl,$ks,$ki) = (keys %{$$L{MAIN}}, keys %{$$S{MAIN}}, keys %{$$I{MAIN}});
  #check point
  # if($$L{MAIN}{$kl}{A}[0] eq "LOTU_000080"){
  #   print "checkpoint. Pause\n";
  # }
  my %TMP;
  my @derep = (keys %{$$L{ALL}}, keys %{$$S{ALL}}, keys %{$$I{ALL}});
  # dereplicate
  @derep = grep { ++$TMP{$_} < 2 } @derep;
  my ($pickTax,$maxScore,$maxOrder,$pickGenus,$maxGscore,$maxGorder,%NAMEID) = ("",0,9999,"",0,9999,);
  my $identityAchieveGenusThreshold = 0;
  foreach my $tax (@derep){
    #log positions info from each database:
    next unless $tax;
    my %regions;
    #SSU:
    foreach my $unit("SSU","ITS","LSU"){
      my $U = ($unit eq "SSU")?\%$S:($unit eq "ITS")?\%$I:\%$L;
      if($$U{ALL}{$tax}{A}){
        $MHASH{$tax}{$unit}{A}   = $$U{ALL}{$tax}{A};
        $MHASH{$tax}{$unit}{ord} = $$U{ALL}{$tax}{order};
        $MHASH{$tax}{$unit}{std} = ($$U{ALL}{$tax}{A}[8]<$$U{ALL}{$tax}{A}[9])?"+":"-";
        $MHASH{$tax}{$unit}{sf}  = ($$U{ALL}{$tax}{A}[8]<$$U{ALL}{$tax}{A}[9])?$$U{ALL}{$tax}{A}[8]:$$U{ALL}{$tax}{A}[9];
        $MHASH{$tax}{$unit}{st}  = ($$U{ALL}{$tax}{A}[8]<$$U{ALL}{$tax}{A}[9])?$$U{ALL}{$tax}{A}[9]:$$U{ALL}{$tax}{A}[8];
        $MHASH{$tax}{$unit}{err} = $$U{ALL}{$tax}{A}[4] + 2* $$U{ALL}{$tax}{A}[5];
      }
    }

    $MHASH{$tax}{qlen} = ($$S{ALL}{$tax}{A}[12])?$$S{ALL}{$tax}{A}[12]:($$L{ALL}{$tax}{A}[12])?$$L{ALL}{$tax}{A}[12]:$$I{ALL}{$tax}{A}[12];
    $MHASH{$tax}{LOTU}  = ($$S{ALL}{$tax}{A}[0])?$$S{ALL}{$tax}{A}[0]:($$L{ALL}{$tax}{A}[0])?$$L{ALL}{$tax}{A}[0]:$$I{ALL}{$tax}{A}[0];
    $MHASH{$tax}{taxN}  = ($$S{ALL}{$tax}{rank})?$$S{ALL}{$tax}{rank}:($$L{ALL}{$tax}{rank})?$$L{ALL}{$tax}{rank}:$$I{ALL}{$tax}{rank};
    $MHASH{$tax}{taxN1}  = ($$S{ALL}{$tax}{rank_1})?$$S{ALL}{$tax}{rank_1}:($$L{ALL}{$tax}{rank_1})?$$L{ALL}{$tax}{rank_1}:$$I{ALL}{$tax}{rank_1};
    $NAMEID{$MHASH{$tax}{taxN1}}{pickTax}  ||= "";
    $NAMEID{$MHASH{$tax}{taxN1}}{minScore} ||= 5000;
    $NAMEID{$MHASH{$tax}{taxN1}}{maxOrder} ||= 999;

    $MHASH{$tax}{mLen}  = $$L{ALL}{$tax}{A}[3] + $$S{ALL}{$tax}{A}[3] + $$I{ALL}{$tax}{A}[3];
    $MHASH{$tax}{mis}   = $$L{ALL}{$tax}{A}[4] + $$S{ALL}{$tax}{A}[4] + $$I{ALL}{$tax}{A}[4];
    $MHASH{$tax}{gap}   = $$L{ALL}{$tax}{A}[5] + $$S{ALL}{$tax}{A}[5] + $$I{ALL}{$tax}{A}[5];
    if($MHASH{$tax}{mLen} == 0){
      die "$MHASH{$tax}{mLen}\n";
    }
    my @orders = sort (($$L{ALL}{$tax}{order})?$$L{ALL}{$tax}{order}:999,($$S{ALL}{$tax}{order})?$$S{ALL}{$tax}{order}:999,($$I{ALL}{$tax}{order})?$$I{ALL}{$tax}{order}:999);
    $MHASH{$tax}{order} = $orders[0];
    $MHASH{$tax}{ident} = 100 - 100 * ($MHASH{$tax}{mis} + 2 * $MHASH{$tax}{gap}) / $MHASH{$tax}{mLen};
    $MHASH{$tax}{bit}   = $$L{ALL}{$tax}{A}[11] + $$S{ALL}{$tax}{A}[11] + $$I{ALL}{$tax}{A}[11];
    $MHASH{$tax}{score} = ($MHASH{$tax}{taxN}=~/metagenome|unidentified|uncultured|unknown/)? $MHASH{$tax}{bit} - 10: $MHASH{$tax}{bit};
    $MHASH{$tax}{qcov}  = 100 * $MHASH{$tax}{mLen} / $MHASH{$tax}{qlen};
    $MHASH{$tax}{scov}  = 100 * $MHASH{$tax}{mLen} / ($$L{ALL}{$tax}{A}[13] + $$S{ALL}{$tax}{A}[13] + $$I{ALL}{$tax}{A}[13]);
    $MHASH{$tax}{scovL} = ($$L{ALL}{$tax}{A}[13])?(100 * $$L{ALL}{$tax}{A}[3] / $$L{ALL}{$tax}{A}[13]):0;
    $MHASH{$tax}{scovS} = ($$S{ALL}{$tax}{A}[13])?(100 * $$S{ALL}{$tax}{A}[3] / $$S{ALL}{$tax}{A}[13]):0;
    $MHASH{$tax}{scovI} = ($$I{ALL}{$tax}{A}[13])?(100 * $$I{ALL}{$tax}{A}[3] / $$I{ALL}{$tax}{A}[13]):0;
    $MHASH{$tax}{anno}  = ($$L{ALL}{$tax}{A}[1])?"$$L{ALL}{$tax}{A}[15]\t$$L{ALL}{$tax}{A}[16]\t$$L{ALL}{$tax}{A}[17]":
    ($$S{ALL}{$tax}{A}[1])?"$$S{ALL}{$tax}{A}[15]\t$$S{ALL}{$tax}{A}[16]\t$$S{ALL}{$tax}{A}[17]":
    "$$I{ALL}{$tax}{A}[15]\t$$I{ALL}{$tax}{A}[16]\t$$I{ALL}{$tax}{A}[17]";
    #pick top genus if identity > 97:
    if($MHASH{$tax}{ident}>97){
      if($identityAchieveGenusThreshold == 0){
        $identityAchieveGenusThreshold = 1;
        ($pickTax,$maxScore,$maxOrder) = ("",0,9999);
      }
      if($MHASH{$tax}{score} == $maxScore){
        $pickGenus = ($MHASH{$tax}{order} < $maxGorder)?$MHASH{$tax}{taxN1}:$pickGenus;
        $maxGscore = ($MHASH{$tax}{order} < $maxGorder)?$MHASH{$tax}{score}:$maxGscore;
        $maxGorder = ($MHASH{$tax}{order} < $maxGorder)?$MHASH{$tax}{order}:$maxGorder;
      }else{
        $pickGenus = ($MHASH{$tax}{score} > $maxGscore)?$MHASH{$tax}{taxN1}:$pickGenus;
        $maxGscore = ($MHASH{$tax}{score} > $maxGscore)?$MHASH{$tax}{score}:$maxGscore;
        $maxGorder = ($MHASH{$tax}{score} > $maxGscore)?$MHASH{$tax}{order}:$maxGorder;
      }
      #now pick top species under top genus:
      my $mis = $MHASH{$tax}{mis} + $MHASH{$tax}{gap};
      $NAMEID{$pickGenus}{pickTax}  = ( $mis < $NAMEID{$pickGenus}{minScore})?$tax:$NAMEID{$pickGenus}{pickTax};
      $NAMEID{$pickGenus}{minScore} = ( $mis < $NAMEID{$pickGenus}{minScore})?$mis:$NAMEID{$pickGenus}{minScore};
      $NAMEID{$pickGenus}{maxOrder} = ( $mis < $NAMEID{$pickGenus}{minScore})?$MHASH{$tax}{order}:$NAMEID{$pickGenus}{maxOrder};
    }elsif($identityAchieveGenusThreshold == 0){
      if($MHASH{$tax}{score} == $maxScore){
        $pickTax = ($MHASH{$tax}{order} < $maxOrder)?$tax:$pickTax;
        $maxScore= ($MHASH{$tax}{order} < $maxOrder)?$MHASH{$tax}{score}:$maxScore;
        $maxOrder= ($MHASH{$tax}{order} < $maxOrder)?$MHASH{$tax}{order}:$maxOrder;
      }else{
        $pickTax = ($MHASH{$tax}{score} > $maxScore)?$tax:$pickTax;
        $maxScore= ($MHASH{$tax}{score} > $maxScore)?$MHASH{$tax}{score}:$maxScore;
        $maxOrder= ($MHASH{$tax}{score} > $maxScore)?$MHASH{$tax}{order}:$maxOrder;
      }
    }
    #&verbose("[DEBUG] $MHASH{$tax}{LOTU}: $MHASH{$tax}{taxN}\t$MHASH{$tax}{ident}\t$MHASH{$tax}{mis}\t$MHASH{$tax}{gap}\t$MHASH{$tax}{bit}\t$MHASH{$tax}{score}\n")
  }
  my $tax = ($identityAchieveGenusThreshold)?$NAMEID{$pickGenus}{pickTax}:$pickTax;

  print OUT sprintf("%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",$MHASH{$tax}{LOTU},$MHASH{$tax}{taxN},
  $MHASH{$tax}{ident}, $MHASH{$tax}{mLen}, $MHASH{$tax}{mis}, $MHASH{$tax}{gap}, $MHASH{$tax}{bit}, $MHASH{$tax}{qlen}, $MHASH{$tax}{qcov},
  $MHASH{$tax}{scov}, $MHASH{$tax}{scovI}, $MHASH{$tax}{scovS}, $MHASH{$tax}{scovL}, $MHASH{$tax}{anno});
#  }
  #return %MHASH if needed
  return(%MHASH);
}

###############################################################################
# Summary SSU, LSU and/or UNITE annotation for each bead
################################################################################
sub usage4cbr {
  my $msg = shift;
  print <<USAGE;
$msg
usage [v0.1] :
  $0 f2b -m cbr -i SSU,LSU,UNITE -o output
    -i  LSU,SSU,ITS annotation file
    -l  rank
    -U  UNITE annotation file
    -o  output
    -v  verbose
    -h  show help info
USAGE
}

sub run_cbr {
  &verbose("[log] Mode [CBR] (Check Bead Region) start ... \n");
  my @files = split (",",$inf);
  &verbose("[err] [CBR] needs two or three files. Pls check\n") & die $! if @files < 2;
  open my $IN1, "<$files[0]" or die $!;
  open my $IN2, "<$files[1]" or die $!;
  open my $IN3, ($files[2])?"<$files[2]":"</dev/null"; # If no ITS input, open an empty filehandle;
  open BEAD, ($obd)?">$obd":($out)?">$out.bead":"> -" or die $!;
  my (%BEADINFO);
  my ($count,$countB,%LSU,%SSU,%ITS,%PLSU,%PSSU,%PITS) = (0,0,);
  my %EMPTY = ("MAIN" => {"NONE" => {"A" => [(0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0)]}});
  chomp(my $read1 = <$IN1>);
  my @inf1 = split (/\t/,$read1);
  my $no1 = $inf1[0];#($inf1[0] =~ /^(\S+\d+)_k\d+_(\d+)_.*_C(\d+|-)_\(/)?"$1.$2.$3":($inf1[0])?$inf1[0]:<EOF>;
  $LSU{NEXT}{A} = \@inf1; $LSU{NEXT}{N} = $no1; $LSU{FILEHANDLE} = $IN1;
  chomp(my $read2 = <$IN2>);
  my @inf2 = split (/\t/,$read2);
  my $no2 = $inf2[0];#($inf2[0] =~ /^(\S+\d+)_k\d+_(\d+)_.*_C(\d+|-)_\(/)?"$1.$2.$3":($inf2[0])?$inf2[0]:<EOF>;
  $SSU{NEXT}{A} = \@inf2; $SSU{NEXT}{N} = $no2; $SSU{FILEHANDLE} = $IN2;
  chomp(my $read3 = <$IN3>);
  my @inf3 = split (/\t/,$read3);
  my $no3 = $inf3[0];#($inf3[0] =~ /^(\S+\d+)_k\d+_(\d+)_.*_C(\d+|-)_\(/)?"$1.$2.$3":($inf3[0])?$inf3[0]:<EOF>;
  $ITS{NEXT}{A} = \@inf3; $ITS{NEXT}{N} = $no3; $ITS{FILEHANDLE} = $IN3;

  my @pick = &checkOrder($no1,$no2,$no3);
  while($pick[1] || $pick[2] || $pick[3]){
    # read LSU annotations of 1 OTU if allowed
    if($pick[1]){
      %LSU  = &read1clip(\%LSU);
      %PLSU = %LSU;
      $no1  = $LSU{NEXT}{N};
    }else{
      %PLSU = %EMPTY;
    }
    # read SSU annotations of 1 OTU if allowed
    if($pick[2]){
      %SSU  = &read1clip(\%SSU);
      %PSSU = %SSU;
      $no2  = $SSU{NEXT}{N};
    }else{
      %PSSU = %EMPTY;
    }
    # read ITS annotations of 1 OTU if allowed
    if($pick[3]){
      %ITS  = &read1clip(\%ITS);
      %PITS = %ITS;
      $no3  = $ITS{NEXT}{N};
    }else{
      %PITS = %EMPTY;
    }
    #summary annotation of this OTU
    if($pick[4] =~ "BI00000022"){
      my $debug = 1;
    }
    %{$BEADINFO{$pick[4]}} = &summaryAnno2(\%PLSU,\%PSSU,\%PITS,$pick[4]);
    $count ++;
    #verbose
    if($count % 1000 == 0 ){
      &verbose(sprintf("[CBR] Processing %6d th clips and %6d beads\r",$count,$countB));
    }

    #summary bead if all clips read:
    my @nextpick = &checkOrder($no1,$no2,$no3);
    my $bid_current = ($pick[4]  =~ /(\S+\d+)\.(\d+)\.(\S+)/)?$1:($pick[4]  =~ /^(.*BI\d+)_k\d+_(\d+)_.*C([0-9]+|-)_/)?$1:$pick[4];
    my $bid_next = ($nextpick[4] =~ /(\S+\d+)\.(\d+)\.(\S+)/)?$1:($nextpick[4] =~ /^(.*BI\d+)_k\d+_(\d+)_.*C([0-9]+|-)_/)?$1:$nextpick[4];
    if ($bid_next ne $bid_current){
      &summaryBead($bid_current,\%BEADINFO);
      %BEADINFO = ();
      $countB ++;
    }


    #prepare for next loop
    %EMPTY = ("MAIN" => {"NONE" => {"A" => [(0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0)]}});
    @pick = @nextpick;
  }
  &verbose(sprintf("[CBR] Processing %6d th clips and %6d beads [done]\n",$count,$countB));
  close OUT;
  &verbose("[CBR] Mode [F2B] done ... \n");
}

sub read1clip {
  my $HS = shift;
  #load first line:
  my ($no1,$IN,@inf1) = ($$HS{NEXT}{N},$$HS{FILEHANDLE},@{$$HS{NEXT}{A}});
  #rest hash
  %$HS = ();
  my ($no,@inf,$order) = ($no1,@inf1,0);
  while($no eq $no1){
    #treat previous record
    my @tags = (exists $$HS{MAIN})?("ALL"):("ALL","MAIN");
    my @ranks = split /;/, $inf[17];
    $order ++;
    foreach my $tag (@tags){
      $$HS{$tag}{$inf[15]}{RefID}    ||= $inf[1];
      $$HS{$tag}{$inf[15]}{taxID}    ||= $inf[15];
      $$HS{$tag}{$inf[15]}{taxName}  ||= $inf[13];
      $$HS{$tag}{$inf[15]}{bitScore} ||= $inf[11];
      $$HS{$tag}{$inf[15]}{identity} ||= $inf[2];
      $$HS{$tag}{$inf[15]}{rank}     ||= $ranks[$rank];
      $$HS{$tag}{$inf[15]}{rank_1}   ||= ($ranks[$rank-1])?$ranks[$rank-1]:$ranks[$rank];
      $$HS{$tag}{$inf[15]}{order}    ||= $order;
      $$HS{$tag}{$inf[15]}{A}        ||= [@inf];
    }

    #read new record
    chomp(my $read = <$IN>);
    if($read eq <EOF>){
      @inf = (<EOF>);
      $no = <EOF>;
      close $IN;
      last;
    }else{
      @inf = split (/\t/,$read);
      $no = $inf[0]; #($inf[0] =~ /^(\S+\d+)_k\d+_(\d+)_.*_C(\d+|-)_\(/)?"$1.$2.$3":($inf[0])?$inf[0]:<EOF>;
    }
  }
  $$HS{NEXT}{A} = [@inf];
  $$HS{NEXT}{N} = $no;
  $$HS{FILEHANDLE} = $IN;
  return(%$HS);
}

sub summaryAnno2{
	my ($L,$S,$I,$CID) = @_;
	my %MHASH;
	my ($kl,$ks,$ki) = (keys %{$$L{MAIN}}, keys %{$$S{MAIN}}, keys %{$$I{MAIN}});
	#check point
	# if($$L{MAIN}{$kl}{A}[0] eq "LOTU_000080"){
	#   print "checkpoint. Pause\n";
	# }
	my %TMP;
	my @derep = (keys %{$$L{ALL}}, keys %{$$S{ALL}}, keys %{$$I{ALL}});
	# dereplicate
	@derep = grep { ++$TMP{$_} < 2 } @derep;
	my ($pickTax,$maxScore,$maxOrder,$pickGenus,$maxGscore,$maxGorder,%NAMEID) = ("",0,9999,"",0,9999,);
	my $identityAchieveGenusThreshold = 0;
	foreach my $tax (@derep){
		#log positions info from each database:
		next unless $tax;
		my %regions;
		#SSU:
		foreach my $unit("SSU","ITS","LSU"){
			my $U = ($unit eq "SSU")?\%$S:($unit eq "ITS")?\%$I:\%$L;
			if($$U{ALL}{$tax}{A}){
				$MHASH{$tax}{$unit}{A}   = $$U{ALL}{$tax}{A};
				$MHASH{$tax}{$unit}{ord} = $$U{ALL}{$tax}{order};
				$MHASH{$tax}{$unit}{std} = ($$U{ALL}{$tax}{A}[8]<$$U{ALL}{$tax}{A}[9])?"+":"-";
				$MHASH{$tax}{$unit}{sf}  = ($$U{ALL}{$tax}{A}[8]<$$U{ALL}{$tax}{A}[9])?$$U{ALL}{$tax}{A}[8]:$$U{ALL}{$tax}{A}[9];
				$MHASH{$tax}{$unit}{st}  = ($$U{ALL}{$tax}{A}[8]<$$U{ALL}{$tax}{A}[9])?$$U{ALL}{$tax}{A}[9]:$$U{ALL}{$tax}{A}[8];
				$MHASH{$tax}{$unit}{err} = $$U{ALL}{$tax}{A}[4] + 2* $$U{ALL}{$tax}{A}[5];
			}
		}

		$MHASH{$tax}{qlen} = ($$S{ALL}{$tax}{A}[12])?$$S{ALL}{$tax}{A}[12]:($$L{ALL}{$tax}{A}[12])?$$L{ALL}{$tax}{A}[12]:$$I{ALL}{$tax}{A}[12];
		$MHASH{$tax}{LOTU}  = ($$S{ALL}{$tax}{A}[0])?$$S{ALL}{$tax}{A}[0]:($$L{ALL}{$tax}{A}[0])?$$L{ALL}{$tax}{A}[0]:$$I{ALL}{$tax}{A}[0];
		$MHASH{$tax}{taxN}  = ($$S{ALL}{$tax}{rank})?$$S{ALL}{$tax}{rank}:($$L{ALL}{$tax}{rank})?$$L{ALL}{$tax}{rank}:$$I{ALL}{$tax}{rank};
		$MHASH{$tax}{taxN1}  = ($$S{ALL}{$tax}{rank_1})?$$S{ALL}{$tax}{rank_1}:($$L{ALL}{$tax}{rank_1})?$$L{ALL}{$tax}{rank_1}:$$I{ALL}{$tax}{rank_1};
		$NAMEID{$MHASH{$tax}{taxN1}}{pickTax}  ||= "";
		$NAMEID{$MHASH{$tax}{taxN1}}{minScore} ||= 5000;
		$NAMEID{$MHASH{$tax}{taxN1}}{maxOrder} ||= 999;

		$MHASH{$tax}{mLen}  = $$L{ALL}{$tax}{A}[3] + $$S{ALL}{$tax}{A}[3] + $$I{ALL}{$tax}{A}[3];
		$MHASH{$tax}{mis}   = $$L{ALL}{$tax}{A}[4] + $$S{ALL}{$tax}{A}[4] + $$I{ALL}{$tax}{A}[4];
		$MHASH{$tax}{gap}   = $$L{ALL}{$tax}{A}[5] + $$S{ALL}{$tax}{A}[5] + $$I{ALL}{$tax}{A}[5];
		if($MHASH{$tax}{mLen} == 0){
			die "$MHASH{$tax}{mLen}\n";
		}
		my @orders = sort (($$L{ALL}{$tax}{order})?$$L{ALL}{$tax}{order}:999,($$S{ALL}{$tax}{order})?$$S{ALL}{$tax}{order}:999,($$I{ALL}{$tax}{order})?$$I{ALL}{$tax}{order}:999);
		$MHASH{$tax}{order} = $orders[0];
		$MHASH{$tax}{ident} = 100 - 100 * ($MHASH{$tax}{mis} + 2 * $MHASH{$tax}{gap}) / $MHASH{$tax}{mLen};
		$MHASH{$tax}{bit}   = $$L{ALL}{$tax}{A}[11] + $$S{ALL}{$tax}{A}[11] + $$I{ALL}{$tax}{A}[11];
		$MHASH{$tax}{score} = ($MHASH{$tax}{taxN}=~/metagenome|unidentified|uncultured|unknown/)? $MHASH{$tax}{bit} - 10: $MHASH{$tax}{bit};
		$MHASH{$tax}{qcov}  = 100 * $MHASH{$tax}{mLen} / $MHASH{$tax}{qlen};
		$MHASH{$tax}{scov}  = 100 * $MHASH{$tax}{mLen} / ($$L{ALL}{$tax}{A}[13] + $$S{ALL}{$tax}{A}[13] + $$I{ALL}{$tax}{A}[13]);
		$MHASH{$tax}{scovL} = ($$L{ALL}{$tax}{A}[13])?(100 * $$L{ALL}{$tax}{A}[3] / $$L{ALL}{$tax}{A}[13]):0;
		$MHASH{$tax}{scovS} = ($$S{ALL}{$tax}{A}[13])?(100 * $$S{ALL}{$tax}{A}[3] / $$S{ALL}{$tax}{A}[13]):0;
		$MHASH{$tax}{scovI} = ($$I{ALL}{$tax}{A}[13])?(100 * $$I{ALL}{$tax}{A}[3] / $$I{ALL}{$tax}{A}[13]):0;
		$MHASH{$tax}{anno}  = ($$L{ALL}{$tax}{A}[1])?"$$L{ALL}{$tax}{A}[15]\t$$L{ALL}{$tax}{A}[16]\t$$L{ALL}{$tax}{A}[17]":
		($$S{ALL}{$tax}{A}[1])?"$$S{ALL}{$tax}{A}[15]\t$$S{ALL}{$tax}{A}[16]\t$$S{ALL}{$tax}{A}[17]":
		"$$I{ALL}{$tax}{A}[15]\t$$I{ALL}{$tax}{A}[16]\t$$I{ALL}{$tax}{A}[17]";
		#pick top genus if identity > 97:
		if($MHASH{$tax}{ident}>97){
			if($identityAchieveGenusThreshold == 0){
				$identityAchieveGenusThreshold = 1;
				($pickTax,$maxScore,$maxOrder) = ("",0,9999);
			}
			if($MHASH{$tax}{score} == $maxScore){
				$pickGenus = ($MHASH{$tax}{order} < $maxGorder)?$MHASH{$tax}{taxN1}:$pickGenus;
				$maxGscore = ($MHASH{$tax}{order} < $maxGorder)?$MHASH{$tax}{score}:$maxGscore;
				$maxGorder = ($MHASH{$tax}{order} < $maxGorder)?$MHASH{$tax}{order}:$maxGorder;
			}else{
				$pickGenus = ($MHASH{$tax}{score} > $maxGscore)?$MHASH{$tax}{taxN1}:$pickGenus;
				$maxGscore = ($MHASH{$tax}{score} > $maxGscore)?$MHASH{$tax}{score}:$maxGscore;
				$maxGorder = ($MHASH{$tax}{score} > $maxGscore)?$MHASH{$tax}{order}:$maxGorder;
			}
			#now pick top species under top genus:
			my $mis = $MHASH{$tax}{mis} + $MHASH{$tax}{gap};
			$NAMEID{$pickGenus}{pickTax}  = ( $mis < $NAMEID{$pickGenus}{minScore})?$tax:$NAMEID{$pickGenus}{pickTax};
			$NAMEID{$pickGenus}{minScore} = ( $mis < $NAMEID{$pickGenus}{minScore})?$mis:$NAMEID{$pickGenus}{minScore};
			$NAMEID{$pickGenus}{maxOrder} = ( $mis < $NAMEID{$pickGenus}{minScore})?$MHASH{$tax}{order}:$NAMEID{$pickGenus}{maxOrder};
		}elsif($identityAchieveGenusThreshold == 0){
			if($MHASH{$tax}{score} == $maxScore){
				$pickTax = ($MHASH{$tax}{order} < $maxOrder)?$tax:$pickTax;
				$maxScore= ($MHASH{$tax}{order} < $maxOrder)?$MHASH{$tax}{score}:$maxScore;
				$maxOrder= ($MHASH{$tax}{order} < $maxOrder)?$MHASH{$tax}{order}:$maxOrder;
			}else{
				$pickTax = ($MHASH{$tax}{score} > $maxScore)?$tax:$pickTax;
				$maxScore= ($MHASH{$tax}{score} > $maxScore)?$MHASH{$tax}{score}:$maxScore;
				$maxOrder= ($MHASH{$tax}{score} > $maxScore)?$MHASH{$tax}{order}:$maxOrder;
			}
		}
		#&verbose("[DEBUG] $MHASH{$tax}{LOTU}: $MHASH{$tax}{taxN}\t$MHASH{$tax}{ident}\t$MHASH{$tax}{mis}\t$MHASH{$tax}{gap}\t$MHASH{$tax}{bit}\t$MHASH{$tax}{score}\n")
	}

	#part 2: summary

	my $BID = ($CID =~ /^(.*BI\d+)_k\d+_(\d+)_.*C([0-9]+|-)_/)?$1:$CID;
  my %BI;
	%{$BI{$CID}} = %MHASH;

  my $B = \%BI;
	my (%QP,%BINFO,$fi,%RE,%MAP,%SCORE,%SCOREBAK,%SUPPORTCNO2TAX);
	return () if $BID eq "";
	my @cNos = ($CID);
	foreach my $cNo (@cNos){
		my @s = sort {$$B{$cNo}{$b}{score} <=> $$B{$cNo}{$a}{score}} keys %{$$B{$cNo}};
		my ($s1,$s2) = ($s[0],$s[1]);
		$BINFO{piece} ++;
		$BINFO{bLen}  += $$B{$cNo}{$s1}{qlen};
		# sort regions on clips
		my (%TMPSCORE,%TMPSCOREBAK) = ();
		for(my $i=0;$i<@s;$i++){
			my $tax = $s[$i]; # the ref ID
			my $sp = $$B{$cNo}{$tax}{taxN};
			my $ge = $$B{$cNo}{$tax}{taxN1};
			my $spPenalty = ($sp =~ /metagenome|unidentified|uncultured|unknown| sp\.$/)?10:0;
			my $gePenalty = ($ge =~ /metagenome|unidentified|uncultured|unknown/)?10:0;
			# pick best/alternative score for genus and species(if multiple subspecies were found) level:
			foreach my $unit ("SSU","ITS","LSU"){
				#genus
				# if($cNo eq "BI00004826.1.-" && $sp eq "Bacillus sp. 17376" && $unit eq "LSU"){
				#   print STDERR "checkpoint\n";
				# }
				if($$B{$cNo}{$tax}{$unit}{A}[11] && $$B{$cNo}{$tax}{$unit}{A}[11] > $TMPSCORE{G}{$ge}{$unit}{score}){
					if($$B{$cNo}{$tax}{$unit}{A}[2] > 97){
						$TMPSCORE{G}{$ge}{$unit}{score} = $$B{$cNo}{$tax}{$unit}{A}[11];
						$TMPSCORE{G}{$ge}{$unit}{len} = $$B{$cNo}{$tax}{$unit}{A}[3];
						$TMPSCORE{G}{$ge}{$unit}{mis} = $$B{$cNo}{$tax}{$unit}{A}[4];
						$TMPSCORE{G}{$ge}{$unit}{gap} = $$B{$cNo}{$tax}{$unit}{A}[5];
						$TMPSCORE{G}{$ge}{$unit}{$cNo} = $tax; #needed in later summarise part
					}else{
						$TMPSCOREBAK{G}{$ge}{$unit}{score} = $$B{$cNo}{$tax}{$unit}{A}[11];
						$TMPSCOREBAK{G}{$ge}{$unit}{len} = $$B{$cNo}{$tax}{$unit}{A}[3];
						$TMPSCOREBAK{G}{$ge}{$unit}{mis} = $$B{$cNo}{$tax}{$unit}{A}[4];
						$TMPSCOREBAK{G}{$ge}{$unit}{gap} = $$B{$cNo}{$tax}{$unit}{A}[5];
						$TMPSCOREBAK{G}{$ge}{$unit}{$cNo} = $tax;
					}
				}
				#sp
				if($$B{$cNo}{$tax}{$unit}{A}[11] && $$B{$cNo}{$tax}{$unit}{A}[11] > $TMPSCORE{S}{$ge}{$sp}{$unit}{score}){
					if($$B{$cNo}{$tax}{$unit}{A}[2] > 99){
						$TMPSCORE{S}{$ge}{$sp}{$unit}{score} = $$B{$cNo}{$tax}{$unit}{A}[11];
						$TMPSCORE{S}{$ge}{$sp}{$unit}{$cNo} = $tax; #needed in later summarise part
					}else{
						$TMPSCOREBAK{S}{$ge}{$sp}{$unit}{score} = $$B{$cNo}{$tax}{$unit}{A}[11];
						$TMPSCOREBAK{S}{$ge}{$sp}{$unit}{$cNo} = $tax;
					}
				}
				# record supports of clips to this tax
				if($$B{$cNo}{$tax}{$unit}{A}[11] && $$B{$cNo}{$tax}{$unit}{A}[2] > 99){
					$SUPPORTCNO2TAX{$cNo}{$sp}{S} += $$B{$cNo}{$tax}{$unit}{A}[11];
					$SUPPORTCNO2TAX{$cNo}{$sp}{L} += $$B{$cNo}{$tax}{$unit}{A}[3];
					$SUPPORTCNO2TAX{$cNo}{$sp}{C} = 1;
					$SUPPORTCNO2TAX{$cNo}{$sp}{G} = $ge;
				}
				# For unique tax records:
				$SCORE{T}{$tax}{score} += $$B{$cNo}{$tax}{$unit}{A}[11];
			}
		}
		#sum total score for each tax from all ref database:
		foreach my $ge (sort keys %{$TMPSCORE{G}}){
			foreach my $unit ("SSU","ITS","LSU"){
				#sum bitscore
				$SCORE{G}{$ge}  += $TMPSCORE{G}{$ge}{$unit}{score};
			}
		}
		foreach my $ge (sort keys %{$TMPSCORE{S}}){
			foreach my $unit ("SSU","ITS","LSU"){
				foreach my $sp (sort keys %{$TMPSCORE{S}{$ge}}){
					# if($cNo eq "BI00004826.1.-" && $sp eq "Bacillus sp. 17376" && $unit eq "LSU"){
					#   print STDERR "checkpoint\n";
					# }
					$SCORE{S}{$ge}{$sp}{score}  += $TMPSCORE{S}{$ge}{$sp}{$unit}{score};
					# calculate more info for species level
					if(my $tax = $TMPSCORE{S}{$ge}{$sp}{$unit}{$cNo}){
						$SCORE{S}{$ge}{$sp}{len} += $$B{$cNo}{$tax}{$unit}{A}[3];
						$SCORE{S}{$ge}{$sp}{mis} += $$B{$cNo}{$tax}{$unit}{A}[4];
						$SCORE{S}{$ge}{$sp}{gap} += $$B{$cNo}{$tax}{$unit}{A}[5];
						$SCORE{S}{$ge}{$sp}{$unit}{$cNo} = $tax;
					}
				}
			}
		}
		#bak
		foreach my $ge (sort keys %{$TMPSCOREBAK{G}}){
			foreach my $unit ("SSU","ITS","LSU"){
				#sum bitscore
				$SCOREBAK{G}{$ge} += $TMPSCOREBAK{G}{$ge}{$unit}{SCOREBAK};
			}
		}
		foreach my $ge (sort keys %{$TMPSCOREBAK{S}}){
			foreach my $unit ("SSU","ITS","LSU"){
				foreach my $sp (sort keys %{$TMPSCOREBAK{S}{$ge}}){
					$SCOREBAK{S}{$ge}{$sp}{score}  += $TMPSCOREBAK{S}{$ge}{$sp}{$unit}{score};
					# calculate more info for species level
					if(my $tax = $TMPSCOREBAK{S}{$ge}{$sp}{$unit}{$cNo}){
						$SCOREBAK{S}{$ge}{$sp}{len} += $$B{$cNo}{$tax}{$unit}{A}[3];
						$SCOREBAK{S}{$ge}{$sp}{mis} += $$B{$cNo}{$tax}{$unit}{A}[4];
						$SCOREBAK{S}{$ge}{$sp}{gap} += $$B{$cNo}{$tax}{$unit}{A}[5];
						$SCOREBAK{S}{$ge}{$sp}{$unit}{$cNo} = $tax;
					}
				}
			}
		}
	}
	# summary
	############################################################################
	# pick top tax among top species which belonging to the top genus:
	my (@genus,@topGenus,@species,@topSpecies,$SCO,$ge_1st,$sp_1st,@ge_2nds,@sp_2nds);
	@genus = sort {$SCORE{G}{$b} <=> $SCORE{G}{$a}} keys %{$SCORE{G}};
	if( @genus && $SCORE{G}{$genus[0]} ){
		@topGenus = grep {$SCORE{G}{$_} == $SCORE{G}{$genus[0]} } @genus;
	}else{
		@genus = sort (sort {$SCOREBAK{G}{$b} <=> $SCOREBAK{G}{$a}} keys %{$SCOREBAK{G}});
		@topGenus = grep {$SCOREBAK{G}{$_} == $SCOREBAK{G}{$genus[0]} } @genus;
	}
	# search high confidence score:
	$SCO = \%SCOREBAK;
	foreach my $ge (@topGenus){
		@species = sort {$SCORE{S}{$ge}{$b}{score} <=> $SCORE{S}{$ge}{$a}{score}} keys %{$SCORE{S}{$ge}};
		if(keys %{$SCORE{S}{$ge}} && $SCORE{S}{$ge}{$species[0]}{score}){
			$SCO = \%SCORE; last;
		}
	}
	foreach my $ge (@topGenus){
		@species = sort {$$SCO{S}{$ge}{$b}{score} <=> $$SCO{S}{$ge}{$a}{score}} keys %{$$SCO{S}{$ge}};
		@topSpecies = grep {$$SCO{S}{$ge}{$_}{score} == $$SCO{S}{$ge}{$species[0]}{score} } @species;
		foreach my $sp (@topSpecies){
			next unless $$SCO{S}{$ge}{$sp}{score};
			#$sp_1st ||= $sp;
			if($$SCO{S}{$ge}{$sp}{score} > $$SCO{S}{$ge_1st}{$sp_1st}{score}){
				@sp_2nds = ($sp_1st) if $sp_1st; $sp_1st = $sp;
				@ge_2nds = ($ge_1st) if $ge_1st; $ge_1st = $ge;
			}elsif($$SCO{S}{$ge}{$sp}{score} == $$SCO{S}{$ge_1st}{$sp_1st}{score}){
				my $ident = ($$SCO{S}{$ge}{$sp}{mis} + 2 * $$SCO{S}{$ge}{$sp}{gap})/$$SCO{S}{$ge}{$sp}{len};
				my $id1st = ($$SCO{S}{$ge_1st}{$sp_1st}{mis} + 2 * $$SCO{S}{$ge_1st}{$sp_1st}{gap})/$$SCO{S}{$ge_1st}{$sp_1st}{len};
				# if(@sp_2nds && $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len} == 0 ){
				#   print STDERR "checkpoint\n";
				# }
				my $id2nd = (@sp_2nds && $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len}> 0)?(($$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{mis} + 2 * $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{gap})/$$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len}):1;
				if($ident < $id1st ){
					@sp_2nds = ($sp_1st) if $sp_1st; $sp_1st = $sp;
					@ge_2nds = ($ge_1st) if $ge_1st; $ge_1st = $ge;
				}elsif($ident < $id2nd){
					@sp_2nds = ($sp);
					@ge_2nds = ($ge);
				}elsif($ident == $id2nd){
					push @sp_2nds, $sp;
					push @ge_2nds, $ge;
				}
			}elsif($$SCO{S}{$ge}{$sp} > $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}){
				@sp_2nds = ($sp);
				@ge_2nds = ($ge);
			}elsif($$SCO{S}{$ge}{$sp} == $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}){
				my $ident = ($$SCO{S}{$ge}{$sp}{mis} + 2 * $$SCO{S}{$ge}{$sp}{gap})/$$SCO{S}{$ge}{$sp}{len};
				my $id2nd = (@sp_2nds)?(($$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{mis} + 2 * $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{gap})/$$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len}):1;
				if($ident < $id2nd){
					@sp_2nds = ($sp);
					@ge_2nds = ($ge);
				}elsif($ident == $id2nd){
					push @sp_2nds, $sp;
					push @ge_2nds, $ge;
				}
			}
		}
	}
	# check whether sp_1st is unique:
	my ($num2nds,$uniqAnno,$maxDiffLv2nds,$score1st,$score2nd,$id1st,$id2nd) = ((scalar @sp_2nds),"","",);
	$score1st = $$SCO{S}{$ge_1st}{$sp_1st}{score};
	$score2nd = (@sp_2nds)?$$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{score}:0;

	$id1st = ($$SCO{S}{$ge_1st}{$sp_1st}{len})?sprintf("%.2f",100-100*($$SCO{S}{$ge_1st}{$sp_1st}{mis} + 2 * $$SCO{S}{$ge_1st}{$sp_1st}{gap})/$$SCO{S}{$ge_1st}{$sp_1st}{len}):0;
	$id2nd = ($$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len})?sprintf("%.2f",100-100*(($$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{mis} + 2 * $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{gap})/$$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len})):0;
	if($score1st == $score2nd && $id1st == $id2nd){
		$maxDiffLv2nds = "sameSp."; $uniqAnno = "multi";
		for(my $i=0; $i<@sp_2nds; $i++){
			if($ge_2nds[$i] ne $ge_1st){
				$maxDiffLv2nds = "diffGe.";
				last;
			}elsif($sp_2nds[$i] ne $sp_1st){
				$maxDiffLv2nds = "diffSp.";
			}
		}
	}else{
		$maxDiffLv2nds = "diffGe.";  $uniqAnno = "unique";
		for(my $i=0; $i<@sp_2nds; $i++){
			if($sp_2nds[$i] eq $sp_1st){
				$maxDiffLv2nds = "sameSp.";
				last;
			}elsif($ge_2nds[$i] eq $ge_1st){
				$maxDiffLv2nds = "sameGe.";
			}
		}
	}
	#estimate the possibility of hybridize bead:
	my %HYB;
	foreach my $cNo (@cNos){
		my @support_sps = keys %{$SUPPORTCNO2TAX{$cNo}};
		if(not defined $SUPPORTCNO2TAX{$cNo}{$sp_1st}{C}){
			foreach my $sp (@support_sps){
				if($SUPPORTCNO2TAX{$cNo}{$sp}{G} ne $ge_1st && $SUPPORTCNO2TAX{$cNo}{$sp}{G} !~ /metagenome|unidentified|uncultured|unknown/){
					$HYB{hybCNo} ++; $HYB{hybLen} += $SUPPORTCNO2TAX{$cNo}{$sp}{L};
					last;
				}
			}
		}
		$HYB{CNoCount} ++
	}
	$HYB{hybCNo} ||= 0; $HYB{hybLen}||=0;
	$HYB{hybPct} = sprintf("%.2f",100* $HYB{hybCNo} / $HYB{CNoCount});
	#sum multi clips into a bead scale:
	my @binfs = ($CID,$BINFO{bLen},$BINFO{piece},0);
	my @res=("",0,0,0,0,1,0,1,0,0,0,$BINFO{bLen},0,"","","","","");
	my ($pUnit,$pTax) = ("SSU","");
	foreach my $unit("SSU","ITS","LSU"){
		my ($sLen,%REG,@ctypes,@desc) =(0,,,);
		if ($pUnit ne $unit){$res[13] .=","};
		foreach my $cNo (@cNos){
			my $tax = ($SCORE{S}{$ge_1st}{$sp_1st}{score})?$SCORE{S}{$ge_1st}{$sp_1st}{$unit}{$cNo}:$SCOREBAK{S}{$ge_1st}{$sp_1st}{$unit}{$cNo};
			next unless exists $$B{$cNo}{$tax}{$unit}{A}[0];
			$sLen ||= $$B{$cNo}{$tax}{$unit}{A}[13];
			$$B{$cNo}{$tax}{LOTU}=~/^(BI\d+)_k\d+_(\d+)_.+_(C[0-9\-]+)_/;
			my $ctype = "$2$3";
			$binfs[3] ++;                                                             #desc: hits of references
			$res[0]  ||=$sp_1st;                                                      #desc: taxonomy name
			$res[1]  += $$B{$cNo}{$tax}{$unit}{err};                                  #desc: identity (intermediate value)
			$res[2]  += $$B{$cNo}{$tax}{$unit}{A}[3];                                 #sum of mapped base (may contained overlaped region in multi fragments)
			$res[3]  += $$B{$cNo}{$tax}{$unit}{A}[4];                                 #desc: mismatch
			$res[4]  += $$B{$cNo}{$tax}{$unit}{A}[5];                                 #desc: gap
			$res[13] .= $$B{$cNo}{$tax}{$unit}{std};                                  #desc: qeury strand
			$res[10] += $$B{$cNo}{$tax}{$unit}{A}[11];                                #desc: bitScore
			$res[17] ||= $$B{$cNo}{$tax}{$unit}{A}[1];
			push @ctypes, "$ctype:$$B{$cNo}{$tax}{$unit}{A}[6]-$$B{$cNo}{$tax}{$unit}{A}[7]";     #desc: query Map desc
			#push @desc, (($pTax ne $tax)?"$$B{$cNo}{$tax}{taxN}:":"")."$$B{$cNo}{$tax}{$unit}{sf}-$$B{$cNo}{$tax}{$unit}{st}";     #desc: subject Map desc
			push @desc, (($pTax ne $tax)?"$tax:":"")."$$B{$cNo}{$tax}{$unit}{sf}-$$B{$cNo}{$tax}{$unit}{st}";     #desc: subject Map desc
			$REG{$$B{$cNo}{$tax}{$unit}{sf}} = $$B{$cNo}{$tax}{$unit}{st};
			$pTax = $tax;
		}
		my $qMapLen = &calCovRegion(\%REG); $res[6] += $qMapLen;
		$res[12] += $sLen; $res[8] += $sLen;
		$res[14] .= (@ctypes)?"$unit(".join(",",@ctypes).")":"";                    #desc of clip types covered in each subunits' refs
		$res[15] .= (@desc)?"$unit(".join(",",@desc).")":"";                        #desc of ref region be aligned
		$res[16] .= ($sLen>0)?sprintf("$unit(%d)",100*$qMapLen/$sLen):"";
		$pUnit = $unit;
	}
	$res[1] = $res[2]?sprintf("%.2f",100*(1-$res[1]/$res[2])):0;

	print OUT (join("\t",join("|",@binfs),@res)."\t".join("\t", $uniqAnno, $num2nds, $maxDiffLv2nds, $score1st, $score2nd, $id1st, $id2nd, join(",",((@sp_2nds)?@sp_2nds:""))));
	print OUT "\t$HYB{hybCNo}\t$HYB{hybLen}\t$HYB{hybPct}\n";

	return(%MHASH);
}

sub summaryBead{
  my ($BID,$B) = @_;
  my (%QP,%BINFO,$fi,%RE,%MAP,%SCORE,%SCOREBAK,%SUPPORTCNO2TAX);
  return () if $BID eq "";
  my @cNos = sort {$a<=>$b} keys %$B;
  foreach my $cNo (@cNos){
    my @s = sort {$$B{$cNo}{$b}{score} <=> $$B{$cNo}{$a}{score}} keys %{$$B{$cNo}};
    my ($s1,$s2) = ($s[0],$s[1]);
    $BINFO{piece} ++;
    $BINFO{bLen}  += $$B{$cNo}{$s1}{qlen};
    # sort regions on clips
    my (%TMPSCORE,%TMPSCOREBAK) = ();
    for(my $i=0;$i<@s;$i++){
      my $tax = $s[$i]; # the ref ID
      my $sp = $$B{$cNo}{$tax}{taxN};
      my $ge = $$B{$cNo}{$tax}{taxN1};
      my $spPenalty = ($sp =~ /metagenome|unidentified|uncultured|unknown| sp\.$/)?10:0;
      my $gePenalty = ($ge =~ /metagenome|unidentified|uncultured|unknown/)?10:0;
      # pick best/alternative score for genus and species(if multiple subspecies were found) level:
      foreach my $unit ("SSU","ITS","LSU"){
        #genus
        # if($cNo eq "BI00004826.1.-" && $sp eq "Bacillus sp. 17376" && $unit eq "LSU"){
        #   print STDERR "checkpoint\n";
        # }
        if($$B{$cNo}{$tax}{$unit}{A}[11] && $$B{$cNo}{$tax}{$unit}{A}[11] > $TMPSCORE{G}{$ge}{$unit}{score}){
          if($$B{$cNo}{$tax}{$unit}{A}[2] > 97){
            $TMPSCORE{G}{$ge}{$unit}{score} = $$B{$cNo}{$tax}{$unit}{A}[11];
            $TMPSCORE{G}{$ge}{$unit}{len} = $$B{$cNo}{$tax}{$unit}{A}[3];
            $TMPSCORE{G}{$ge}{$unit}{mis} = $$B{$cNo}{$tax}{$unit}{A}[4];
            $TMPSCORE{G}{$ge}{$unit}{gap} = $$B{$cNo}{$tax}{$unit}{A}[5];
            $TMPSCORE{G}{$ge}{$unit}{$cNo} = $tax; #needed in later summarise part
          }else{
            $TMPSCOREBAK{G}{$ge}{$unit}{score} = $$B{$cNo}{$tax}{$unit}{A}[11];
            $TMPSCOREBAK{G}{$ge}{$unit}{len} = $$B{$cNo}{$tax}{$unit}{A}[3];
            $TMPSCOREBAK{G}{$ge}{$unit}{mis} = $$B{$cNo}{$tax}{$unit}{A}[4];
            $TMPSCOREBAK{G}{$ge}{$unit}{gap} = $$B{$cNo}{$tax}{$unit}{A}[5];
            $TMPSCOREBAK{G}{$ge}{$unit}{$cNo} = $tax;
          }
        }
        #sp
        if($$B{$cNo}{$tax}{$unit}{A}[11] && $$B{$cNo}{$tax}{$unit}{A}[11] > $TMPSCORE{S}{$ge}{$sp}{$unit}{score}){
          if($$B{$cNo}{$tax}{$unit}{A}[2] > 99){
            $TMPSCORE{S}{$ge}{$sp}{$unit}{score} = $$B{$cNo}{$tax}{$unit}{A}[11];
            $TMPSCORE{S}{$ge}{$sp}{$unit}{$cNo} = $tax; #needed in later summarise part
          }else{
            $TMPSCOREBAK{S}{$ge}{$sp}{$unit}{score} = $$B{$cNo}{$tax}{$unit}{A}[11];
            $TMPSCOREBAK{S}{$ge}{$sp}{$unit}{$cNo} = $tax;
          }
        }
        # record supports of clips to this tax
        if($$B{$cNo}{$tax}{$unit}{A}[11] && $$B{$cNo}{$tax}{$unit}{A}[2] > 99){
          $SUPPORTCNO2TAX{$cNo}{$sp}{S} += $$B{$cNo}{$tax}{$unit}{A}[11];
          $SUPPORTCNO2TAX{$cNo}{$sp}{L} += $$B{$cNo}{$tax}{$unit}{A}[3];
          $SUPPORTCNO2TAX{$cNo}{$sp}{C} = 1;
          $SUPPORTCNO2TAX{$cNo}{$sp}{G} = $ge;
        }
        # For unique tax records:
        $SCORE{T}{$tax}{score} += $$B{$cNo}{$tax}{$unit}{A}[11];
      }
    }
    #sum total score for each tax from all ref database:
    foreach my $ge (sort keys %{$TMPSCORE{G}}){
      foreach my $unit ("SSU","ITS","LSU"){
        #sum bitscore
        $SCORE{G}{$ge}  += $TMPSCORE{G}{$ge}{$unit}{score};
      }
    }
    foreach my $ge (sort keys %{$TMPSCORE{S}}){
      foreach my $unit ("SSU","ITS","LSU"){
        foreach my $sp (sort keys %{$TMPSCORE{S}{$ge}}){
          # if($cNo eq "BI00004826.1.-" && $sp eq "Bacillus sp. 17376" && $unit eq "LSU"){
          #   print STDERR "checkpoint\n";
          # }
          $SCORE{S}{$ge}{$sp}{score}  += $TMPSCORE{S}{$ge}{$sp}{$unit}{score};
          # calculate more info for species level
          if(my $tax = $TMPSCORE{S}{$ge}{$sp}{$unit}{$cNo}){
            $SCORE{S}{$ge}{$sp}{len} += $$B{$cNo}{$tax}{$unit}{A}[3];
            $SCORE{S}{$ge}{$sp}{mis} += $$B{$cNo}{$tax}{$unit}{A}[4];
            $SCORE{S}{$ge}{$sp}{gap} += $$B{$cNo}{$tax}{$unit}{A}[5];
            $SCORE{S}{$ge}{$sp}{$unit}{$cNo} = $tax;
          }
        }
      }
    }
    #bak
    foreach my $ge (sort keys %{$TMPSCOREBAK{G}}){
      foreach my $unit ("SSU","ITS","LSU"){
        #sum bitscore
        $SCOREBAK{G}{$ge} += $TMPSCOREBAK{G}{$ge}{$unit}{SCOREBAK};
      }
    }
    foreach my $ge (sort keys %{$TMPSCOREBAK{S}}){
      foreach my $unit ("SSU","ITS","LSU"){
        foreach my $sp (sort keys %{$TMPSCOREBAK{S}{$ge}}){
          $SCOREBAK{S}{$ge}{$sp}{score}  += $TMPSCOREBAK{S}{$ge}{$sp}{$unit}{score};
          # calculate more info for species level
          if(my $tax = $TMPSCOREBAK{S}{$ge}{$sp}{$unit}{$cNo}){
            $SCOREBAK{S}{$ge}{$sp}{len} += $$B{$cNo}{$tax}{$unit}{A}[3];
            $SCOREBAK{S}{$ge}{$sp}{mis} += $$B{$cNo}{$tax}{$unit}{A}[4];
            $SCOREBAK{S}{$ge}{$sp}{gap} += $$B{$cNo}{$tax}{$unit}{A}[5];
            $SCOREBAK{S}{$ge}{$sp}{$unit}{$cNo} = $tax;
          }
        }
      }
    }
  }
  # summary
  ############################################################################
  # pick top tax among top species which belonging to the top genus:
  my (@genus,@topGenus,@species,@topSpecies,$SCO,$ge_1st,$sp_1st,@ge_2nds,@sp_2nds);
  @genus = sort {$SCORE{G}{$b} <=> $SCORE{G}{$a}} keys %{$SCORE{G}};
  if( @genus && $SCORE{G}{$genus[0]} ){
    @topGenus = grep {$SCORE{G}{$_} == $SCORE{G}{$genus[0]} } @genus;
  }else{
    @genus = sort (sort {$SCOREBAK{G}{$b} <=> $SCOREBAK{G}{$a}} keys %{$SCOREBAK{G}});
    @topGenus = grep {$SCOREBAK{G}{$_} == $SCOREBAK{G}{$genus[0]} } @genus;
  }
  # search high confidence score:
  $SCO = \%SCOREBAK;
  foreach my $ge (@topGenus){
    @species = sort {$SCORE{S}{$ge}{$b}{score} <=> $SCORE{S}{$ge}{$a}{score}} keys %{$SCORE{S}{$ge}};
    if(keys %{$SCORE{S}{$ge}} && $SCORE{S}{$ge}{$species[0]}{score}){
      $SCO = \%SCORE; last;
    }
  }
  foreach my $ge (@topGenus){
    @species = sort {$$SCO{S}{$ge}{$b}{score} <=> $$SCO{S}{$ge}{$a}{score}} keys %{$$SCO{S}{$ge}};
    @topSpecies = grep {$$SCO{S}{$ge}{$_}{score} == $$SCO{S}{$ge}{$species[0]}{score} } @species;
    foreach my $sp (@topSpecies){
      next unless $$SCO{S}{$ge}{$sp}{score};
      #$sp_1st ||= $sp;
      if($$SCO{S}{$ge}{$sp}{score} > $$SCO{S}{$ge_1st}{$sp_1st}{score}){
        @sp_2nds = ($sp_1st) if $sp_1st; $sp_1st = $sp;
        @ge_2nds = ($ge_1st) if $ge_1st; $ge_1st = $ge;
      }elsif($$SCO{S}{$ge}{$sp}{score} == $$SCO{S}{$ge_1st}{$sp_1st}{score}){
        my $ident = ($$SCO{S}{$ge}{$sp}{mis} + 2 * $$SCO{S}{$ge}{$sp}{gap})/$$SCO{S}{$ge}{$sp}{len};
        my $id1st = ($$SCO{S}{$ge_1st}{$sp_1st}{mis} + 2 * $$SCO{S}{$ge_1st}{$sp_1st}{gap})/$$SCO{S}{$ge_1st}{$sp_1st}{len};
        # if(@sp_2nds && $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len} == 0 ){
        #   print STDERR "checkpoint\n";
        # }
        my $id2nd = ($$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len}>0)?(($$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{mis} + 2 * $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{gap})/$$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len}):1;
        if($ident < $id1st ){
          @sp_2nds = ($sp_1st) if $sp_1st; $sp_1st = $sp;
          @ge_2nds = ($ge_1st) if $ge_1st; $ge_1st = $ge;
        }elsif($ident < $id2nd){
          @sp_2nds = ($sp);
          @ge_2nds = ($ge);
        }elsif($ident == $id2nd){
          push @sp_2nds, $sp;
          push @ge_2nds, $ge;
        }
      }elsif($$SCO{S}{$ge}{$sp} > $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}){
        @sp_2nds = ($sp);
        @ge_2nds = ($ge);
      }elsif($$SCO{S}{$ge}{$sp} == $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}){
        my $ident = ($$SCO{S}{$ge}{$sp}{mis} + 2 * $$SCO{S}{$ge}{$sp}{gap})/$$SCO{S}{$ge}{$sp}{len};
        my $id2nd = (@sp_2nds)?(($$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{mis} + 2 * $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{gap})/$$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len}):1;
        if($ident < $id2nd){
          @sp_2nds = ($sp);
          @ge_2nds = ($ge);
        }elsif($ident == $id2nd){
          push @sp_2nds, $sp;
          push @ge_2nds, $ge;
        }
      }
    }
  }
  # check whether sp_1st is unique:
  my ($num2nds,$uniqAnno,$maxDiffLv2nds,$score1st,$score2nd,$id1st,$id2nd) = ((scalar @sp_2nds),"","",);
  $score1st = $$SCO{S}{$ge_1st}{$sp_1st}{score};
  $score2nd = (@sp_2nds)?$$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{score}:0;
  $id1st = $$SCO{S}{$ge_1st}{$sp_1st}{len}?sprintf("%.2f",100-100*($$SCO{S}{$ge_1st}{$sp_1st}{mis} + 2 * $$SCO{S}{$ge_1st}{$sp_1st}{gap})/$$SCO{S}{$ge_1st}{$sp_1st}{len}):0;
  $id2nd = ($$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len})?sprintf("%.2f",100-100*(($$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{mis} + 2 * $$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{gap})/$$SCO{S}{$ge_2nds[0]}{$sp_2nds[0]}{len})):0;
  if($score1st == $score2nd && $id1st == $id2nd){
    $maxDiffLv2nds = "sameSp."; $uniqAnno = "multi";
    for(my $i=0; $i<@sp_2nds; $i++){
      if($ge_2nds[$i] ne $ge_1st){
        $maxDiffLv2nds = "diffGe.";
        last;
      }elsif($sp_2nds[$i] ne $sp_1st){
        $maxDiffLv2nds = "diffSp.";
      }
    }
  }else{
    $maxDiffLv2nds = "diffGe.";  $uniqAnno = "unique";
    for(my $i=0; $i<@sp_2nds; $i++){
      if($sp_2nds[$i] eq $sp_1st){
        $maxDiffLv2nds = "sameSp.";
        last;
      }elsif($ge_2nds[$i] eq $ge_1st){
        $maxDiffLv2nds = "sameGe.";
      }
    }
  }
  #estimate the possibility of hybridize bead:
  my %HYB;
  foreach my $cNo (@cNos){
    my @support_sps = keys %{$SUPPORTCNO2TAX{$cNo}};
    if(not defined $SUPPORTCNO2TAX{$cNo}{$sp_1st}{C}){
      foreach my $sp (@support_sps){
        if($SUPPORTCNO2TAX{$cNo}{$sp}{G} ne $ge_1st && $SUPPORTCNO2TAX{$cNo}{$sp}{G} !~ /metagenome|unidentified|uncultured|unknown/){
          $HYB{hybCNo} ++; $HYB{hybLen} += $SUPPORTCNO2TAX{$cNo}{$sp}{L};
          last;
        }
      }
    }
    $HYB{CNoCount} ++
  }
  $HYB{hybCNo} ||= 0; $HYB{hybLen}||=0;
  $HYB{hybPct} = sprintf("%.2f",100* $HYB{hybCNo} / $HYB{CNoCount});
  #sum multi clips into a bead scale:
  my @binfs = ($BID,$BINFO{bLen},$BINFO{piece},0);
  my @res=("",0,0,0,0,1,0,1,0,0,0,$BINFO{bLen},0,"","","","","");
  my ($pUnit,$pTax) = ("SSU","");
  foreach my $unit("SSU","ITS","LSU"){
    my ($sLen,%REG,@ctypes,@desc) =(0,,,);
    if ($pUnit ne $unit){$res[13] .=","};
    foreach my $cNo (@cNos){
      my $tax = ($SCORE{S}{$ge_1st}{$sp_1st}{score})?$SCORE{S}{$ge_1st}{$sp_1st}{$unit}{$cNo}:$SCOREBAK{S}{$ge_1st}{$sp_1st}{$unit}{$cNo};
      next unless exists $$B{$cNo}{$tax}{$unit}{A}[0];
      $sLen ||= $$B{$cNo}{$tax}{$unit}{A}[13];
      $$B{$cNo}{$tax}{LOTU}=~/^(BI\d+)_k\d+_(\d+)_.+_(C[0-9\-]+)_/;
      my $ctype = "$2$3";
      $binfs[3] ++;                                                             #desc: hits of references
      $res[0]  ||=$sp_1st;                                                      #desc: taxonomy name
      $res[1]  += $$B{$cNo}{$tax}{$unit}{err};                                  #desc: identity (intermediate value)
      $res[2]  += $$B{$cNo}{$tax}{$unit}{A}[3];                                 #sum of mapped base (may contained overlaped region in multi fragments)
      $res[3]  += $$B{$cNo}{$tax}{$unit}{A}[4];                                 #desc: mismatch
      $res[4]  += $$B{$cNo}{$tax}{$unit}{A}[5];                                 #desc: gap
      $res[13] .= $$B{$cNo}{$tax}{$unit}{std};                                  #desc: qeury strand
      $res[10] += $$B{$cNo}{$tax}{$unit}{A}[11];                                #desc: bitScore
      $res[17] ||= $$B{$cNo}{$tax}{$unit}{A}[1];
      push @ctypes, "$ctype:$$B{$cNo}{$tax}{$unit}{A}[6]-$$B{$cNo}{$tax}{$unit}{A}[7]";     #desc: query Map desc
      #push @desc, (($pTax ne $tax)?"$$B{$cNo}{$tax}{taxN}:":"")."$$B{$cNo}{$tax}{$unit}{sf}-$$B{$cNo}{$tax}{$unit}{st}";     #desc: subject Map desc
      push @desc, (($pTax ne $tax)?"$tax:":"")."$$B{$cNo}{$tax}{$unit}{sf}-$$B{$cNo}{$tax}{$unit}{st}";     #desc: subject Map desc
      $REG{$$B{$cNo}{$tax}{$unit}{sf}} = $$B{$cNo}{$tax}{$unit}{st};
       $pTax = $tax;
    }
    my $qMapLen = &calCovRegion(\%REG); $res[6] += $qMapLen;
    $res[12] += $sLen; $res[8] += $sLen;
    $res[14] .= (@ctypes)?"$unit(".join(",",@ctypes).")":"";                    #desc of clip types covered in each subunits' refs
    $res[15] .= (@desc)?"$unit(".join(",",@desc).")":"";                        #desc of ref region be aligned
    $res[16] .= ($sLen>0)?sprintf("$unit(%d)",100*$qMapLen/$sLen):"";
    $pUnit = $unit;
  }
  $res[1] = $res[2]?sprintf("%.2f",100*(1-$res[1]/$res[2])):0;

  print BEAD (join("\t",join("|",@binfs),@res)."\t".join("\t", $uniqAnno, $num2nds, $maxDiffLv2nds, $score1st, $score2nd, $id1st, $id2nd, join(",",((@sp_2nds)?@sp_2nds:""))));
  print BEAD "\t$HYB{hybCNo}\t$HYB{hybLen}\t$HYB{hybPct}\n";

}

sub calCovRegion{
  my $R = shift;
  my ($cov,$f0,$t0) = (0,,);
  my @f = sort {$a<=>$b} keys %$R;
  for(my $i=0;$i<@f;$i++){
    my $fi = $f[$i];
    my $ti  = $$R{$fi};
    $f0 ||= $fi;
    $t0 = ($ti>$t0)?$ti:$t0;
    if($i eq $#f|| $t0 < $f[$i+1]){
      $cov += $t0 - $f0 + 1;
      #reset fi
      $f0 = "";
    }
  }
  return($cov);
}



################################################################################
# For fungi, summary SSU, LSU and UNITE annotation for each otu
################################################################################
sub usage4otupair {
  my $msg = shift;
  print <<USAGE;
$msg
usage [v0.1] :
  $0 -m otupair -i <lotu fasta>,<taxon tree> -o output
    -i  LOTU fa , taxon tree, id map
    -o  output
    -v  verbose
    -h  show help info
USAGE
}


sub run_otupair {
  &verbose("[otupair] Mode start ... \n");
  my @files = split (",",$inf);
  &verbose("[otupair] needs 3 files. Pls check\n") & die $! if @files != 3;
  open INF, "<$files[0]" or die $!;
  open TAX, "<$files[1]" or die $!;
  open MAP, "<$files[2]" or die $!;
  open OUT, ">$out" or die $!;
  ####
  my (%RANK,%RRANK,%MAP);
  my ($rs,$ss,)=(10,0);
  ($RANK{all},$RANK{root},$RANK{major_clade}) = (0,1,99);
  ($RRANK{0},$RRANK{1},$RRANK{99}) = ("all","root","major_clade");
  foreach my $r ("domain","kingdom","phylum","class","order","family","genus","species","LOTU"){
    $ss = -1;
    foreach my $s ("super","","sub","infra"){
      $RANK{"$s$r"} = $rs + $ss;
      $RRANK{$rs+$ss} = "$s$r";
      $ss ++;
    }
    $rs += 10;
  }
  ####
  &verbose("  Reading taxid ... ");
  while(<MAP>){
    chomp;
    my @s = split /\t/;
    next if $s[1] eq "";
    $MAP{$s[1]}{$s[0]} ++;
    $MAP{$s[0]}{$s[1]} ++;
  }
  close MAP;
  ####
  my(%TREE,%LOTU,%LAST,%UC);
  &verbose("done\n  Reading tax and run vserach ... ");
  my $LB = $/; $/ = ">";
  while(<INF>){
    chomp; next unless $_;
    my @s = split /\n/;
    my $L1 = shift @s;
    my ($id,$tax);
    if($L1 =~/^(.*BI\d+)_k\d+_(\d+)_.*C([0-9]+|-)_/){
      $id = "$1.$2.$3";
      $tax = (keys %{$MAP{$id}})[0];
    }elsif($L1 =~ /^(\S+) (.*)$/){
      $id = $1;
      my @t = split /;/,$2;
      $tax = $t[-1];
    }else{
      $id = $L1;
      $tax = (keys %{$MAP{$id}})[0];
    }
    #$id =~ s/LOTU_//;
    if($id eq "O4BI00000584.1.-"){
      my $debug = 1;
    }
    $LOTU{$id} = join("\n",@s);
    next if exists $MAP{$id};
    $MAP{$tax}{$id} ++;
  }
  close INF; $/ = $LB;
  &verbose("done\n  Reading taxon:\n");
  chomp(my $ramtmp=("-e /dev/shm")?`mktemp -d /dev/shm/otupair.XXXXXXXXXX`:`mktemp -dp $out.XXX`);
  my $summary = 1;
  while(($_ = <TAX>) || $summary ){
    my(@s,@t,$rankScore,$findOTU,@lastRs,$d,$sameRank);
    @lastRs = sort {$a<=>$b} keys %{$LAST{RT}};
    my $sd = $#lastRs;
    if($_){
      chomp;
      if($. == 102673){
        $debug =1 ;
      }
      @s = split /\t/;
      my @st = split /;/,$s[0];
      next if $st[-1] =~/major_clade__/;
      foreach my $i (@st){
        push @t, $i unless $i =~/major_clade__/;
      }
      $rankScore = $RANK{$s[2]};
      #$findOTU = ($t[-1] =~ s/LOTU_//)?1:0;

      $findOTU = ($t[-1] =~ /species__/)?1:0;
      if($findOTU && $t[-2] =~ /genus__/){
        if($t[-1] =~ /species__Alexandrium ostenfeldii/){
          my $debug = 1;
        }
        (my $geName = $t[-2]) =~  s/^genus__//;
        #if($t[-1] =~/$geName (\S+) (.+)/ && $1 !~ /^sp(\.|)$/){
        if($t[-1] =~/$geName (\S+) (.*)/){
          my $spName = "species__$geName $1";
          #record a new species level
          unless($LAST{RT}{$rankScore} eq $spName){
            $LAST{R} = $rankScore;
            $LAST{T} = $spName;
            $LAST{RT}{$rankScore} = $spName;
          }
          #add a level for subspecies
          push @t, "sub$t[-1]";
          $t[-2] = $spName;
          $rankScore = $RANK{"subspecies"};
        }
        #elsif($t[-1] =~/$geName (sp|sp\.)$/){
        #   #not a valid species rank, cautious
        #   $TREE{$LAST{R}}{$LAST{T}}{$findOTU}{$s[-2]} = $rankScore;
        #   pop @t;
        #   $rankScore = $RANK{"species"};
        # }

      }
      $sameRank = "";
      for($d=0;$d<@lastRs;$d ++){
        if($sameRank eq "" && $LAST{RT}{$lastRs[$d]} ne $t[$d]){
          $sameRank =  $lastRs[$d-1];
          $sd = $d-1;
        }
      }
      $LAST{R} = $lastRs[-1];
      $LAST{T} = $t[-2];
    }else{
      $summary = 0;
      $sd = -1;
      #summary after end line
    }
    #summary previouse taxonomies
    for(my $i=$#lastRs;$i>$sd;$i--){
      my $r = $lastRs[$i];
      my $tax = $LAST{RT}{$r};
      if($tax =~ "kingdom__Fungi"){
        $debug = 1;
      }
      unless($tax =~ /uncultured|unknown|unidentified/){

        my @o = sort {$a<=>$b} keys %{$TREE{$r}{$tax}{1}};
        # prepare fasta belonging to this taxono
        open LOTU,"> $ramtmp/lotu.fa" or die $!;
        my $oCount = 0;
        my $firstSeq = "";
        foreach my $o (@o){
          if($o == 14501){
            my $debug = 1;
          }
          foreach my $s (sort keys %{$MAP{$o}}){
            next unless $LOTU{$s};
            print LOTU ">$s\n$LOTU{$s}\n";
            $oCount ++;
            $firstSeq ||= $s;
            last if $oCount > 99;
          }
        }
        my @o2 = sort {$a<=>$b} keys %{$TREE{$r}{$tax}{2}};
        if(@o2){
          foreach my $o (@o2){
            next unless $LOTU{$o};
            print LOTU ">$o\n$LOTU{$o}\n";
            $oCount ++;
            $firstSeq ||= $o;
            last if $oCount > 99;
          }
        }

        close LOTU;
        # run vsearch alignemnt:
        my ($minIdent,@idents,@minOtus) = (101,,);
        if($oCount > 1){
          if($tax =~ / sp/){
            my $debug = 1;
          }
          my $cmd = "vsearch --threads ".(($oCount <= 50)?1:($oCount <= 100)?8:24)." --allpairs_global $ramtmp/lotu.fa --acceptall --uc - 2> /dev/null ";
          &verbose(" Processing $r: $RRANK{$r} | $tax: $oCount LOTUs ... ");
          open UC,"$cmd|" or die $!;
          while(<UC>){
            chomp;
            my @u = split /\t/;
            next unless $u[0] eq "H";
            #next if $u[3] < 80;
            my @ut = sort {$a<=>$b} ($u[8],$u[9]);
            push @idents, $u[3];
            @minOtus = ($u[3]<$minIdent)?($ut[0],$ut[1]):@minOtus;
            $minIdent = ($u[3]<$minIdent)?$u[3]:$minIdent;
          }
          close UC;
          my @sorts = sort {$a<=>$b} @idents;
          my $snum = $#sorts + 1;
          my $min = $sorts[0];
          my $q01 = $sorts[sprintf("%d",$snum/100)];
          my $q05 = $sorts[sprintf("%d",$snum/20)];
          my $q25 = $sorts[sprintf("%d",$snum/4)];
          my $q50 = $sorts[sprintf("%d",$snum/2)];
          my $q75 = $sorts[sprintf("%d",$snum*3/4)];
          my $q95 = $sorts[sprintf("%d",$snum*19/20)];
          my $q99 = $sorts[sprintf("%d",$snum*99/100)];
          print OUT "$RRANK{$r}\t$tax\t$oCount\t$min\t$q01,$q05,$q25,$q50,$q75,$q95,$q99\n";
          &verbose("done\n");
        }
        #pick OTU to next level:

        my $nextR = $lastRs[$i-1];
        my $nextT = $LAST{RT}{$nextR};
        if(@minOtus){
          $TREE{$nextR}{$nextT}{2}{$minOtus[0]} = $nextT;
          $TREE{$nextR}{$nextT}{2}{$minOtus[1]} = $nextT;
        }elsif($firstSeq){
          $TREE{$nextR}{$nextT}{2}{$firstSeq} = $nextT;
        }
      }
      #remove recorded tax:
      delete $LAST{RT}{$r};
    }
    #summary current  taxonomies
    if($sd != $#t ){
      print SATDERR "$.. might missing a rank\n";
    }
    if(exists $LAST{RT}{$rankScore}){
      #duplicate rank found:
      $rankScore +=1;
    }
    $LAST{R} = $rankScore;
    $LAST{T} = $t[-1];
    $LAST{RT}{$rankScore} = $t[-1];
    if(@t - scalar keys %{$LAST{RT}} > 1 ){
      $debug =1;
    }
    $TREE{$LAST{R}}{$LAST{T}}{$findOTU}{$s[-2]} = $rankScore;

  }
  close TAX;
  close OUT;

  &verbose("done\n");
}


################################################################################
# OTU GC content calculator
################################################################################

sub run_otuGC {
  &verbose("[otuGC] Mode start ... \n");
  my @files = split (",",$inf);
  &verbose("[otuGC] needs 3 files. Pls check\n") & die $! if @files != 3;
  open INF, "<$files[0]" or die $!;
  open TAX, "<$files[1]" or die $!;
  open MAP, "<$files[2]" or die $!;
  open OUT, ">$out" or die $!;
  ####
  my (%RANK,%RRANK,%MAP);
  my ($rs,$ss,)=(10,0);
  ($RANK{all},$RANK{root},$RANK{major_clade}) = (0,1,99);
  ($RRANK{0},$RRANK{1},$RRANK{99}) = ("all","root","major_clade");
  foreach my $r ("domain","kingdom","phylum","class","order","family","genus","species","LOTU"){
    $ss = -1;
    foreach my $s ("super","","sub","infra"){
      $RANK{"$s$r"} = $rs + $ss;
      $RRANK{$rs+$ss} = "$s$r";
      $ss ++;
    }
    $rs += 10;
  }
  ####
  &verbose("  Reading taxid ... ");
  while(<MAP>){
    chomp;
    my @s = split /\t/;
    $MAP{$s[1]}{$s[0]} ++;
    $MAP{$s[0]}{$s[1]} ++;
  }
  close MAP;
  ####
  my(%TREE,%LOTU,%LAST,%UC);
  &verbose("done\n  Reading tax ... ");
  my $LB = $/; $/ = ">";
  while(<INF>){
    chomp; next unless $_;
    my @s = split /\n/;
    my $L1 = shift @s;
    my ($id,$tax);
    if($L1 =~ /^(\S+) (.*)$/){
      $id = $1;
      my @t = split /;/,$2;
      $tax = $t[-1];
    }else{
      $id = $L1;
      $tax = (keys %{$MAP{$id}})[0];
    }
    my $str = join("\n",@s);
    my $cN  = length($str);
    my $cGC = $str =~ tr/GC/GC/;
    $LOTU{$id} = $cGC / $cN;
    next if exists $MAP{$id};
    $MAP{$tax}{$id} ++;
  }
  close INF; $/ = $LB;
  &verbose("done\n  Reading taxon:\n");
  my $summary = 1;
  while(($_ = <TAX>) || $summary ){
    my(@s,@t,$rankScore,$findOTU,@lastRs,$d,$sameRank);
    @lastRs = sort {$a<=>$b} keys %{$LAST{RT}};
    my $sd = $#lastRs;
    if($_){
      chomp;
      @s = split /\t/;
      my @st = split /;/,$s[0];
      next if $st[-1] =~/major_clade__/;
      foreach my $i (@st){
        push @t, $i unless $i =~/major_clade__/;
      }
      if($RANK{$s[3]}){
        $rankScore = $RANK{$s[3]};
      }else{
        if($s[3] =~ /CLADE([0-9.]+)/){
          $rankScore = $1 * 100;
        }else{
          $rankScore = 101;
        }
      }

      $findOTU = ($t[-1] =~ /species__/)?1:0;
      if($findOTU && $t[-2] =~ /genus__/){
        if($t[-1] =~ /species__Alexandrium ostenfeldii/){
          my $debug = 1;
        }
        (my $geName = $t[-2]) =~  s/^genus__//;
        #if($t[-1] =~/$geName (\S+) (.+)/ && $1 !~ /^sp(\.|)$/){
        if($t[-1] =~/$geName (\S+) (.*)/){
          my $spName = "species__$geName $1";
          #record a new species level
          unless($LAST{RT}{$rankScore} eq $spName){
            $LAST{R} = $rankScore;
            $LAST{T} = $spName;
            $LAST{RT}{$rankScore} = $spName;
          }
          #add a level for subspecies
          push @t, "sub$t[-1]";
          $t[-2] = $spName;
          $rankScore = $RANK{"subspecies"};
        }
      }
      $sameRank = "";
      for($d=0;$d<@lastRs;$d ++){
        if($sameRank eq "" && $LAST{RT}{$lastRs[$d]} ne $t[$d]){
          $sameRank =  $lastRs[$d-1];
          $sd = $d-1;
        }
      }
      $LAST{R} = $lastRs[-1];
      $LAST{T} = $t[-2];
    }else{
      $summary = 0;
      $sd = -1;
      #summary after end line
    }
    #summary previouse taxonomies
    for(my $i=$#lastRs;$i>$sd;$i--){
      my $r = $lastRs[$i];
      my $tax = $LAST{RT}{$r};
      unless($tax =~ /uncultured|unknown|unidentified/){

        my @o = (sort {$a<=>$b} keys %{$TREE{$r}{$tax}{0}}, sort {$a<=>$b} keys %{$TREE{$r}{$tax}{1}});
        # prepare fasta belonging to this taxono
        my %GCSTAT;
        my $oCount = 0;
        my $firstSeq = "";
        foreach my $o (@o){
          foreach my $s (sort keys %{$MAP{$o}}){
            $GCSTAT{sum} += $LOTU{$s};
            $GCSTAT{num} ++;
            $oCount ++;
            $firstSeq ||= $s;
          }
        }

        if($TREE{$r}{$tax}{2}{num}){
          $GCSTAT{sum} += $TREE{$r}{$tax}{2}{sum};
          $GCSTAT{num} += $TREE{$r}{$tax}{2}{num};
          $oCount ++;
        }

        # run vsearch alignemnt:
        my ($minIdent,@idents,@minOtus) = (101,,);
        if($oCount > 0){
          my $meanGC = $GCSTAT{sum} / $GCSTAT{num};
          print OUT sprintf("%d\t%s\t%d\t%.4f\n",$LAST{RI}{$r},$tax,$GCSTAT{num},$meanGC);
          #&verbose(" Processing $r: $tax | $LAST{RC}{$r} : $oCount LOTUs \n");
        }
        my $nextR = $lastRs[$i-1];
        my $nextT = $LAST{RT}{$nextR};
        #pick OTU to next level:
        $TREE{$nextR}{$nextT}{2}{num} = $GCSTAT{num};
        $TREE{$nextR}{$nextT}{2}{sum} = $GCSTAT{sum};
      }
      #remove recorded tax:
      delete $LAST{RT}{$r};
    }
    #summary current  taxonomies
    if($sd != $#t ){
      print SATDERR "$.. might missing a rank\n";
    }
    if(exists $LAST{RT}{$rankScore}){
      #duplicate rank found:
      $rankScore +=1;
    }
    $LAST{R} = $rankScore;
    $LAST{T} = $t[-1];
    $LAST{RT}{$rankScore} = $t[-1];
    $LAST{RI}{$rankScore} = $s[1];
    $LAST{RC}{$rankScore} = $s[4];
    if(@t - scalar keys %{$LAST{RT}} > 1 ){
      $debug =1;
    }
    $TREE{$LAST{R}}{$LAST{T}}{$findOTU}{$s[1]} = $rankScore;

  }
  close TAX;
  close OUT;

  &verbose("done\n");
}
