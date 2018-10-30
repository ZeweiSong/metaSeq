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
		$nfs++;
	}else{
		last;
	}
}
my $f10 = (int $fns/1000) *100;
close LST;

open IN, "pigz -dc $in|" or die $!;
open INF, ">$out.info" or die $!;
open F0, "> $out/0000.$sfx.fa" or die $!;
#open FX, "> $out/xxxx.$sfx.fa" or die $!;
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
  	 print STDOUT "$id\n$seq\n"; 
  }elsif(not defined $FILES{$bc[0]}{$bc[1]}{'count'}){
    my $ct = 1;
    $FILES{$bc[0]}{$bc[1]}{'count'} = $ct;
    $FILES{$bc[0]}{$bc[1]}{$ct} = "$id\n$seq\n";

    $STAT{'barcodes'} ++;
  }else{
    $FILES{$bc[0]}{$bc[1]}{'count'} ++;
    if($FILES{$bc[0]}{$bc[1]}{'count'} < $cut){
      my $ct = $FILES{$bc[0]}{$bc[1]}{'count'};
      $FILES{$bc[0]}{$bc[1]}{$ct} = "$id\n$seq\n";
    }elsif($FILES{$bc[0]}{$bc[1]}{'count'} == $cut){
      `mkdir -p $out/$bc[0]`;
      open $FILES{$bc[0]}{$bc[1]}{'FILEHANDLE'}, "> $out/$bc[0]/$bc[1].$sfx.fa" or die $!;
      $FHcount ++;

      # Write the pre-member(s)
      for(my $ct==1;$ct<$cut;$ct++){
        $FILES{$bc[0]}{$bc[1]}{'FILEHANDLE'}->print($FILES{$bc[0]}{$bc[1]}{$ct});
        delete $FILES{$bc[0]}{$bc[1]}{$ct};
      }
      $FILES{$bc[0]}{$bc[1]}{'FILEHANDLE'}->print("$id\n$seq\n");

      #Report
      if($FHcount==1 ||$FHcount%$f10==0){
        @timeS =  &getTime(time,$startTimeStamp);
        print STDERR sprintf("[%02d:%02d:%02d] Filehandle opened: %5d\n",$timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0],$FHcount);
      }

      $STAT{'multiHits'} ++;
    }else{
      $FILES{$bc[0]}{$bc[1]}{'FILEHANDLE'}->print("$id\n$seq\n");
    }
  }
}
close IN;

@timeS =  &getTime(time,$startTimeStamp);
print STDERR sprintf("[%02d:%02d:%02d] %d of %d barcodes found multiple (>=%d) reads. Writing rest reads into single file...\n",
  $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0],$STAT{'multiHits'},$STAT{'barcodes'},$cut);

foreach my $b1 (sort keys %FILES){
  foreach my $b123 (sort keys %{$FILES{$b1}}){
    if($FILES{$b1}{$b123}{'count'}>=$cut){
      print INF "$b123\t$FILES{$b1}{$b123}{'count'}\n";
      close $FILES{$b1}{$b123}{'FILEHANDLE'};
    }else{
      print F0 $FILES{$b1}{$b123}{'the1st'};
    }
  }
}
close INF;
close F0;
close FX;

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
