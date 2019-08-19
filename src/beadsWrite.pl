#!/usr/bin/env perl
use strict;

my ($in,$list,$cut,$sfx,$out) = @ARGV;
my (%FILES,%LIST,%STAT);
my $FHcount=0;
my $startTimeStamp = time();
my @timeS;
`mkdir -p $out`;

open LST,"<$list" or die $!;
my $fns=0;
while(<LST>){
	chomp;
	my @str = split;
	if($str[1]>=$cut){
		$LIST{$str[0]} = $str[1];
		$fns++;
	}else{
		last;
	}
}
my $f10 = (int $fns/1000) *100;
close LST;

open IN, "pigz -dc $in|" or die $!;
#open INF, ">$out.info" or die $!;
open F0, "> $out/0000.$sfx.fa" or die $!;
open FX, "> $out/xxxx.$sfx.fa" or die $!;

my(@p);
while(<IN>){
  # Read one sequence info
  chomp(my $id = $_);
  chomp(my $seq= <IN>);
  <IN>;
  chomp(my $qua= <IN>);
  # Detect barcode
  $id =~ /\/(\d+)_(\d+)_(\d+)\//;
  my @bcode= ($1,$2,$3);
  my @bc = ($bcode[0],"$bcode[0]_$bcode[1]_$bcode[2]");
  if($bcode[0] eq "0000"||$bcode[1] eq "0000"||$bcode[2] eq "0000"){
     print F0 "$id\n$seq\n";
  }elsif(not defined $LIST{$bc[1]}){
  	 print FX "$id\n$seq\n";
  }elsif(not defined $FILES{$bc[0]}{$bc[1]}{'count'}){
		close FS if defined $FILES{$p[0]}{$p[1]}{'count'};
		my $ct = 1;
		$FILES{$bc[0]}{$bc[1]}{'count'} = 1;
		`mkdir -p $out/$bc[0]`;
		open FS, "> $out/$bc[0]/$bc[1].$sfx.fa" or die $!;
		print FS "$id\n$seq\n";
		$STAT{'barcodes'} ++;
		@p = @bc;
		$FHcount ++;
  }else{
    $FILES{$bc[0]}{$bc[1]}{'count'} ++;
		print FS "$id\n$seq\n";
  }
}
close IN;
close F0;
close FX;
close FS;

@timeS =  &getTime(time,$startTimeStamp);
print STDERR sprintf("[%02d:%02d:%02d] %d beads fa written. Reporting...\n",
  $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0],$STAT{'barcodes'});

foreach my $b1 (sort keys %FILES){
  foreach my $b123 (sort keys %{$FILES{$b1}}){
      print STDOUT "$b123\t$FILES{$b1}{$b123}{'count'}\n";
  }
}
#close INF;


@timeS =  &getTime(time,$startTimeStamp);
print STDERR sprintf("[%02d:%02d:%02d] All done.\n",
  $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0]);

exit;

sub getTime {
  my $now = shift;
  my $start = shift;
  my ($diff,$d,$h,$m,$s);
  $diff = $now - $start;
  $h = $diff%(3600*24);
  $d = ($diff -$h)/(3600*24);
  $m = $h%3600;
  $h = ($diff -$m)/3600;
  $s = $m%60;
  $m = ($m - $s)/60;
  return($s,$m,$h,$d);
}
