#!/usr/bin/env perl
use strict;
use Getopt::Long;

sub usage {
    print <<USAGE;
usage: perl beadsWrite.pl  --r1 <single fastq> -L <list> -c <cutoff> -s <sufix> -o <output>
       perl beadsWrite.pl  --r1 <read 1 fastq> --r2 <read 2 fastq> -b <barcode> -f <format> -o <output>

default options:
    --r1        single fq file
    -L        barcode frequency file list
    -c        [int] cutoff
    -s        [1|2] sufix
    -p        process
    -o        output dir
extract specific barcodes:
    --r2        read2 fq file
    -b        barcode id
    -f        format [fq|fa|fq+fa]

USAGE
  exit;
};

&usage unless @ARGV >= 10;
my ($in1,$in2,$list,$cut,$sfx,$bcode,$proc,$fmt,$out);
GetOptions(
    "r1:s"        => \$in1,
    "r2:s"        => \$in2,
    "L|list:s"    => \$list,
    "c|cutoff:s"  => \$cut,
    "s|sufix:s"   => \$sfx,
    "b|barcode:s" => \$bcode,
    "p|proc:s"    => \$proc,
    "F|format:s"  => \$fmt,
    "o|outdir:s"  => \$out,
);
$fmt ||= "fq";
$proc||= 1;

my (%FILES,%LIST,%STAT);
my $FHcount=0;
my $startTimeStamp = time();
my @timeS;
`mkdir -p $out`;

if($list){
  &writeAll($in1,$list,$cut,$sfx,$proc,$fmt,$out);
}elsif($bcode){
  &specifcBarcode($in1,$in2,$bcode,$fmt,$out);
}

exit;

###############################
#
###############################
sub specifcBarcode{
  my($in1,$in2,$bcode,$fmt,$out) = @_;
  @timeS =  &getTime(time,$startTimeStamp);
  print STDERR sprintf("[%02d:%02d:%02d] specifcBarcode mode. Start.\n", $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0]);
  my @bcs = split(",",$bcode);
  while(@bcs){
    my $bc = shift @bcs;
    $bc =~ /(\d+)_(\d+)_(\d+)/;
    my @bc = ($1,$bc);
    $LIST{$bc[1]} = $bc[1];
  }
  my $bnum = keys %LIST;
  open R1,"pigz -p $proc -dc $in1|" or die $!;
  open R2,"pigz -p $proc -dc $in2|" or die $!;

  my @p;
  while(<R1>){
    # Read one sequence info
    my $id = $_;

    # Detect barcode
    my @bc = &getIdBcode($id);

    if (defined $LIST{$bc[1]}){
      if($bc[1] ne $p[1]){
        if($FHcount >0 ){
          close O1;
          close O2;
        }
        if($fmt eq "fq"){
          open O1,"|pigz -p $proc > $out/$bc[0]/$bc[1].1.fq.gz" or die $!;
          open O2,"|pigz -p $proc > $out/$bc[0]/$bc[1].2.fq.gz" or die $!;
        }else{
          open O1,"|pigz -p $proc > $out/$bc[0]/$bc[1].1.fa.gz" or die $!;
          open O2,"|pigz -p $proc > $out/$bc[0]/$bc[1].2.fa.gz" or die $!;
        }
        $FHcount ++ ;
        @timeS =  &getTime(time,$startTimeStamp);
        print STDERR sprintf("[%02d:%02d:%02d] Find #%02s# %s. Start writing.\n",
          $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0],$FHcount,$bc[1]);
        @p = @bc;
      }
      if($fmt eq "fq"){
        print O1 $id.<R1>.<R1>.<R1>;
        print O2 <R2>.<R2>.<R2>.<R2>;
      }else{
        print O1 $id.<R1>;
        print O2 <R2>.<R2>;
        <R1>;<R1>;
        <R2>;<R2>;
      }
      $STAT{'barcodes'} ++;
    }else{
      last if $FHcount eq $bnum;
      <R1>;<R1>;<R1>;
      <R2>;<R2>;<R2>;<R2>;
    }
  }
  close R1;
  close R2;
  close O1;
  close O2;
  @timeS =  &getTime(time,$startTimeStamp);
  print STDERR sprintf("[%02d:%02d:%02d] %d beads with %d reads in total written. Done.\n",
    $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0],$FHcount,$STAT{'barcodes'});

}

sub writeAll {
  my ($in1,$list,$cut,$sfx,$proc,$fmt,$out) = @_;
  @timeS =  &getTime(time,$startTimeStamp);
  print STDERR sprintf("[%02d:%02d:%02d] WriteAll mode. Start.\n", $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0]);
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
  my $f10 = (int $fns/10);
  $f10 = ($f10==0)?1:$f10;
  close LST;
  open IN, "pigz -p $proc -dc $in1|" or die $!;
  if($fmt =~ /fq/){
    open Q0, "|pigz -p $proc > $out/0000.$sfx.fq.gz" or die $!;
    open QX, "|pigz -p $proc > $out/xxxx.$sfx.fq.gz" or die $!;
  }
  if($fmt =~ /fa/){
    open A0, "> $out/0000.$sfx.fa" or die $!;
    open AX, "> $out/xxxx.$sfx.fa" or die $!;
  }

  my(@p);
  while(<IN>){
    # Read one sequence info
    my $id = $_;
    my $seq = <IN>;
    my $qua = <IN>.<IN>;
    # Detect barcode
    my @bc = &getIdBcode($id);

    if (defined $LIST{$bc[1]}){
      if($bc[1] ne $p[1]){
        if($FHcount > 0 ){
          if($fmt =~ /fq/){close QT};
          if($fmt =~ /fa/){close AT};
        }
        `mkdir -p $out/$bc[0]`;
        if($fmt =~ /fq/){open QT,"|pigz -p $proc > $out/$bc[0]/$bc[1].$sfx.fq.gz" or die $!;}
        if($fmt =~ /fa/){open AT,"> $out/$bc[0]/$bc[1].$sfx.fa" or die $!;}
        $FHcount ++ ;
        if($FHcount % $f10 == 0){
          @timeS =  &getTime(time,$startTimeStamp);
          print STDERR sprintf("[%02d:%02d:%02d] Find #%04s / %04s beads. writing...\n",
            $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0],$FHcount,$fns);
        }
        @p = @bc;
      }
      if($fmt =~ /fq/){print QT $id.$seq.$qua};
      if($fmt =~ /fq/){print AT $id.$seq};
      $STAT{$bc[1]} ++;
    }elsif($bc[1] =~ /0000/){
      if($fmt =~ /fq/){print Q0 $id.$seq.$qua;}
      if($fmt =~ /fa/){print A0 $id.$seq;}
    }else{
      if($fmt =~ /fq/){print QX $id.$seq.$qua;}
      if($fmt =~ /fa/){print AX $id.$seq;}
    }
  }

  close IN;
  if($fmt =~ /fq/){
    close QT;
    close Q0;
    close QX;
  }
  if($fmt =~ /fa/){
    close AT;
    close A0;
    close AX;
  }

  @timeS =  &getTime(time,$startTimeStamp);
  print STDERR sprintf("[%02d:%02d:%02d] %d beads written. Reporting...\n",
    $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0],$FHcount);

  foreach my $b1 (sort keys %STAT){
    print STDOUT "$b1\t$STAT{$b1}\n";
  }
  #close INF;

  @timeS =  &getTime(time,$startTimeStamp);
  print STDERR sprintf("[%02d:%02d:%02d] All done.\n",
    $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0]);

}


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

sub getIdBcode {
  my $id = shift;
  $id =~ /\/(\d+)_(\d+)_(\d+)\//;
  my @bcode= ($1,$2,$3);
  my @bc = ($bcode[0],"$bcode[0]_$bcode[1]_$bcode[2]");
  return(@bc);
}
