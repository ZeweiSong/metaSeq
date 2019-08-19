#!/usr/bin/env perl
# Stat identity of annotation in blast6 format.
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
  $0 -r reference -a blast6 -o output
    -r  reference
    -a  alignment in blast6 format
    -o  output filename
    -v  verbose mode
    -h  show help info
USAGE
}

my ($ref,$anno,$out,$verbose,$help);
GetOptions(
  "r=s" => \$ref,
  "a=s" => \$anno,
  "o=s" => \$out,
  "v" => \$verbose,
  "h|help|?" => \$help,
);
&usage && exit if $help;

open REF,"<$ref" or die "Can not open $ref. $!\n";
open ANN,($anno)?"<$anno":"-" or die "Can not open $anno. $!\n";
open OUT,">$out" or die "Can not open $out. $!\n";

my (%ANN);
while(<REF>){
  chomp;
  my @L = split (/\t/,$_);
  my @a = split (/\|/,$L[0]);
  # identity annoation
  $ANN{$a[1]}{'attr'}{'taxonomy'} = $a[0];
  $ANN{$a[1]}{'attr'}{'ID2'} = $a[2];
  $ANN{$a[1]}{'attr'}{'type'} = $a[3];
  my @TAX=split (';',$a[4]);
  my $parent = "root";
  for(my $i=0;$i<@TAX;$i++){
    my @tax = split('__',$TAX[$i]);
    #%ANN{ID}{ATTR} = value;
    #k(indom) -> p(hylum) -> c(lass) -> o(rder)
    # -> f(amily) -> g(enus) -> s(pecies)
    $ANN{$a[1]}{$parent}{$tax[0]} = $tax[1];
    $parent = $tax[1];
  }
  # identity length
  $L[1] =~ s/ bp\.//;
  $ANN{$a[1]}{'attr'}{'length'} = $L[1];
  # Identity TIS regions
  for(my $i=2;$i<@L;$i++){
    my @reg = split (": ",$L[$i]);
    $ANN{$a[1]}{'region'}{$reg[0]} = $reg[1];
  }
}
close REF;

my (%SC,%CT); #Define the score hash
my @lv = ("k","p","c","o","f","g","s");

my ($OTU,@L,%WEIGHT);
#first record
#my @L = split (/\t/,<ANN>);
#$OTU{1} = \@L;
my $end = 1;
my $prv = "";
while(defined($_=<ANN>) || $end){
  chomp;

  #special treatment for last line
  if(defined $_){
    @L = split (/\t/,$_);
  }else{
    $end = 0;
    @L = ("","");
  }
  # read all hits of each OTU
  if($L[1] ne $prv && $. > 1){
    # Summary last OTU score
    for(my $r=0;$r<@$OTU;$r++){
      # pick one record of an OTU
      my @LI = @{$$OTU[$r]};
      my @T = split (/;/,$LI[0]);
      my ($size,$iden,$cov) = (($T[1]=~s/size=//),$LI[2],$LI[3]);
      my ($rID,$parent) = ($L[1],"root");
      for(my $i=0;$i<@lv;$i++){
        # socre for each taxonomy level
        my $tax = (defined $ANN{$rID}{$parent}{$lv[$i]})?$ANN{$rID}{$parent}{$lv[$i]}:"unidentified";
        my $weight = $WEIGHT{$i}{$tax}/$WEIGHT{$i}{'sumAll'};
        $SC{$parent}{$tax} += $iden*$cov * $weight;
        $CT{$parent}{$tax} += $weight;
        $parent = ($tax eq "unidentified")?"$parent\_$tax":$tax;
      }
      # turn to next record
    }
    # init stat of each level
   ($OTU, %WEIGHT) = ();
  }
  # Treat new record
  push @$OTU, \@L;
  my ($rID,$parent) = ($L[1],"root");
  for(my $i=0;$i<@lv;$i++){
    my $tax = (defined $ANN{$rID}{$parent}{$lv[$i]})?$ANN{$rID}{$parent}{$lv[$i]}:"unidentified";
    $WEIGHT{$i}{$tax} ++;
    $WEIGHT{$i}{'sumAll'} ++;
    $parent = ($tax eq "unidentified")?"$parent\_$tax":$tax;
  }
  #prepare for next line
  $prv = $L[1];
}
close ANN;

foreach my $k (sort {$SC{'root'}{$b} <=> $SC{'root'}{$a}} keys %{$SC{'root'}} ){
  #print OUT "$CT{'root'}{$k}\t$SC{'root'}{$k}\tk__$k\n";
  printf OUT ("%d\t%d\t%s\n",$CT{'root'}{$k},$SC{'root'}{$k},"k__$k") if $CT{'root'}{$k} > 0;
  foreach my $p (sort {$SC{$k}{$b} <=> $SC{$k}{$a}} keys %{$SC{$k}} ){
    (my $sp = $p ) =~ s/$k\_//;
    printf OUT ("%d\t%d\t%s\n",$CT{$k}{$p},$SC{$k}{$p},"k__$k|p__$sp") if $CT{$k}{$p} >0 ;
    foreach my $c (sort {$SC{$p}{$b} <=> $SC{$p}{$a}} keys %{$SC{$p}} ){
      (my $sc = $c ) =~ s/$p\_//;
      printf OUT ("%d\t%d\t%s\n",$CT{$p}{$c},$SC{$p}{$c},"k__$k|p__$p|c_$sc") if $CT{$p}{$c} >0 ;
      foreach my $o (sort {$SC{$c}{$b} <=> $SC{$c}{$a}} keys %{$SC{$c}} ){
        (my $so = $o ) =~ s/$c\_//;
        printf OUT ("%d\t%d\t%s\n",$CT{$c}{$o},$SC{$c}{$o},"k__$k|p__$p|c_$c|o__$so") if $CT{$c}{$o} >0 ;
        foreach my $f (sort {$SC{$o}{$b} <=> $SC{$o}{$a}} keys %{$SC{$o}} ){
          (my $sf = $f ) =~ s/$o\_//;
          printf OUT ("%d\t%d\t%s\n",$CT{$o}{$f},$SC{$o}{$f},"k__$k|p__$p|c_$c|o__$o|f__$sf") if $CT{$o}{$f} >0 ;
          foreach my $g (sort {$SC{$f}{$b} <=> $SC{$f}{$a}} keys %{$SC{$f}} ){
            (my $sg = $g ) =~ s/$f\_//;
            printf OUT ("%d\t%d\t%s\n",$CT{$f}{$g},$SC{$f}{$g},"k__$k|p__$p|c_$c|o__$o|f__$f|g__$sg") if $CT{$f}{$g} >0 ;
            foreach my $s (sort {$SC{$g}{$b} <=> $SC{$g}{$a}} keys %{$SC{$g}} ){
              (my $ss = $s ) =~ s/$g\_//;
              printf OUT ("%d\t%d\t%s\n",$CT{$g}{$s},$SC{$g}{$s},"k__$k|p__$p|c_$c|o__$o|f__$f|g__$g|s__$ss") if $CT{$g}{$s} >0 ;
            }
          }
        }
      }
    }
  }
}

close OUT;
exit;
