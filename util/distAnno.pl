#!/usr/bin/perl
use strict;

unless (@ARGV){
  print "Usage: perl distAnno.pl <dist> <bead anno> > <output>\n";
  exit;
}

my($in,$anno) = @ARGV;
my(%HS,%HS2);
open I,"<$anno";
while(<I>){
  next unless $_ =~ /^\d/;
  my @a=split;
  $HS{$a[0]} = $a[1];
}
close I;
open I, "<$in";
while(<I>){
  my @a=split;
  next if $a[0]==$a[1];
  my @b=sort ($HS{$a[0]},$HS{$a[1]});
  $HS2{$b[0]}{$b[1]}{sprintf("%.4f",$a[2])} ++;
}
close I;
print "tax1\ttax2\tdistance\tcount\n";
foreach my $t1 (sort keys %HS2){
  foreach my $t2 (sort keys %{$HS2{$t1}}){
    foreach my $dst (sort keys %{$HS2{$t1}{$t2}}){
      printf("%s\t%s\t%f\t%d\n",$t1,$t2,$dst,$HS2{$t1}{$t2}{$dst});
    }
  }
}

exit;
